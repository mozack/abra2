package abra;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.intel.gkl.compression.IntelDeflaterFactory;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;

public class SortedSAMWriter {
	
	private static final int TEMP_COMPRESSION_LEVEL = 1;
	private static final int FINAL_COMPRESSION_LEVEL = 5;
	
	public static final int GENOMIC_RANGE_TO_CACHE = 1000000;
	private static final int UNMAPPED_INDEX = 0;
	private static final int ASYNC_READ_CACHE_SIZE = 1000000;
	
	private Map<String, Integer> chromIdx = new HashMap<String, Integer>();
	private SAMFileWriterFactory writerFactory = new SAMFileWriterFactory();
	
	private SAMFileWriter writers[][];
	private String tempDir;
	private String[] outputFiles;
	private SAMFileHeader[] samHeaders;
	private boolean isKeepTmp;
	
	public SortedSAMWriter(String[] outputFiles, String tempDir, SAMFileHeader[] samHeaders, boolean isKeepTmp) {
	
		this.samHeaders = samHeaders;
		this.outputFiles = outputFiles;
		this.tempDir = tempDir;
		this.isKeepTmp = isKeepTmp;
		
		// TODO: Compare sequence headers across samples
		List<SAMSequenceRecord> sequences = samHeaders[0].getSequenceDictionary().getSequences();
		int idx = 1;
		
		for (SAMSequenceRecord seq : sequences) {
			chromIdx.put(seq.getSequenceName(), idx++);
		}
		
		// Unmapped reads assigned to slot 0
		chromIdx.put("*", UNMAPPED_INDEX);
		
		writerFactory.setUseAsyncIo(false);
		IntelDeflaterFactory intelDeflater = new IntelDeflaterFactory();
		writerFactory.setDeflaterFactory(intelDeflater);
		
		Logger.info("Using intel deflator: " + intelDeflater.usingIntelDeflater());
		
		writers = new SAMFileWriter[outputFiles.length][];
		
		for (int i=0; i<writers.length; i++) {
			writers[i] = new SAMFileWriter[sequences.size()+1];
		}
	}
	
	private void deleteOnExit(File file) {
		if (!isKeepTmp) {
			file.deleteOnExit();
		}
	}
	
	private String getTempFilename(int sampleIdx, int chrom) {
		return String.format("%s/%d.%d.bam", tempDir, sampleIdx, chrom);
	}
	
	public void addAlignment(int sampleIdx, SAMRecord read) {
		int chrom  = chromIdx.get(read.getReferenceName());
		writers[sampleIdx][chrom].addAlignment(read);
	}
	
	public void initChromosome(String chromosome) {
		for (int i=0; i<outputFiles.length; i++) {
			initChromosome(i, chromosome);
		}
	}
	
	private void initChromosome(int sampleIdx, String chromosome) {
		Logger.debug("Writer init: %d, %s", sampleIdx, chromosome);
		int chrom  = chromIdx.get(chromosome);
		writers[sampleIdx][chrom] = writerFactory.makeBAMWriter(samHeaders[sampleIdx], false, new File(getTempFilename(sampleIdx, chrom)), TEMP_COMPRESSION_LEVEL);
	}
	
	private void finishChromosome(int sampleIdx, String chromosome) {
		Logger.debug("Writer finish: %d, %s", sampleIdx, chromosome);
		int chrom  = chromIdx.get(chromosome);
		writers[sampleIdx][chrom].close();
	}
	
	public void finishChromosome(String chromosome) {
		for (int i=0; i<outputFiles.length; i++) {
			finishChromosome(i, chromosome);
		}
	}
	
	public void outputFinal() throws IOException {
		for (int i=0; i<outputFiles.length; i++) {
			Logger.info("Finishing: " + outputFiles[i]);
			writerFactory.setUseAsyncIo(true);
			writerFactory.setAsyncOutputBufferSize(ASYNC_READ_CACHE_SIZE);
			writerFactory.setCompressionLevel(FINAL_COMPRESSION_LEVEL);
			samHeaders[i].setSortOrder(SortOrder.coordinate);
			SAMFileWriter output = writerFactory.makeBAMWriter(samHeaders[i], true, new File(outputFiles[i]), FINAL_COMPRESSION_LEVEL);
			
			for (int chrom=1; chrom<chromIdx.size(); chrom++) {
				processChromosome(output, i, chrom);
			}
			
			processUnmapped(output, i);
			
			output.close();
		}
	}
	
	private void processChromosome(SAMFileWriter output, int sampleIdx, int chrom) throws IOException {
		
		Logger.debug("Final processing for: %d, %d", sampleIdx, chrom);
		
		List<SAMRecord> reads = new ArrayList<SAMRecord>();
		String filename = getTempFilename(sampleIdx, chrom);
		
		File file = new File(filename);
		
		if (file.exists()) {
			deleteOnExit(file);
			
			SamReader reader = SAMRecordUtils.getSamReader(filename);
	
			for (SAMRecord read : reader) {
				reads.add(read);
				
				if (read.getAlignmentStart() - reads.get(0).getAlignmentStart() > GENOMIC_RANGE_TO_CACHE*2) {
					Collections.sort(reads, new SAMCoordinateComparator());
					
					int start = reads.get(0).getAlignmentStart();
					int i = 0;
					while (i < reads.size() && reads.get(i).getAlignmentStart() - start < GENOMIC_RANGE_TO_CACHE) {
						output.addAlignment(reads.get(i));
						i += 1;
					}
					Logger.trace("Reads output: %d", i);
					reads.subList(0, i).clear();
				}
			}
			
			// Output any remaining reads
			Collections.sort(reads, new SAMCoordinateComparator());
			for (SAMRecord read : reads) {
				output.addAlignment(read);
			}
			
			reader.close();
		}
	}
	
	private void processUnmapped(SAMFileWriter output, int sampleIdx) throws IOException {
		String filename = getTempFilename(sampleIdx, UNMAPPED_INDEX);
		File file = new File(filename);
		
		if (file.exists()) {
		
			SamReader reader = SAMRecordUtils.getSamReader(filename);
	
			for (SAMRecord read : reader) {
				output.addAlignment(read);
			}
			
			reader.close();
		}
	}
	
	static class SAMCoordinateComparator implements Comparator<SAMRecord> {

		@Override
		public int compare(SAMRecord o1, SAMRecord o2) {
			return o1.getAlignmentStart()-o2.getAlignmentStart();
		}
	}
	
	public static void main(String[] args) throws IOException {
		String in = args[0];
		String out = args[1];
		String tempDir = args[2];
		
//		String in = "/home/lmose/dev/abra2_dev/sort/0.1.bam";
//		String out = "/home/lmose/dev/abra2_dev/sort/output.bam";
//		String tempDir = "/home/lmose/dev/abra2_dev/sort";
		
		Logger.setLevel("trace");
		
		
		
		SamReader reader = SAMRecordUtils.getSamReader(in);
		
//		SortedSAMWriter writer = new SortedSAMWriter(new String[] { "/home/lmose/dev/abra2_dev/sort/output.bam" }, "/home/lmose/dev/abra2_dev/sort", new SAMFileHeader[] { reader.getFileHeader() });
		
		SortedSAMWriter writer = new SortedSAMWriter(new String[] { out }, tempDir, new SAMFileHeader[] { reader.getFileHeader() }, false);

		/*
		Set<String> chromosomes = new HashSet<String>();
		
		for (SAMRecord read : reader) {
			if (!chromosomes.contains(read.getReferenceName())) {
				writer.initChromosome(0, read.getReferenceName());
				chromosomes.add(read.getReferenceName());
			}
			
			writer.addAlignment(0, read);
		}
		
		for (String chromosome : chromosomes) {
			writer.finishChromsome(0, chromosome);
		}
		*/
		
		writer.outputFinal();
	}
}

package abra;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.ValidationStringency;

public class SortedSAMWriter {
	
	private static final int TEMP_COMPRESSION_LEVEL = 1;
	
	private static final int GENOMIC_RANGE_TO_CACHE = 1000000;
	private static final int UNMAPPED_INDEX = 0;
	
	private Map<String, Integer> chromIdx = new HashMap<String, Integer>();
	private SAMFileWriterFactory writerFactory = new SAMFileWriterFactory();
	private SAMFileWriter writers[][];
	private String tempDir;
	private String[] outputFiles;
	private SAMFileHeader[] samHeaders;
	
	public SortedSAMWriter(String[] outputFiles, String tempDir, SAMFileHeader[] samHeaders) {
		this.samHeaders = samHeaders;
		this.outputFiles = outputFiles;
		this.tempDir = tempDir;
		
		// TODO: Compare sequence headers across samples
		List<SAMSequenceRecord> sequences = samHeaders[0].getSequenceDictionary().getSequences();
		int idx = 1;
		
		for (SAMSequenceRecord seq : sequences) {
			chromIdx.put(seq.getSequenceName(), idx++);
		}
		
		// Unmapped reads assigned to slot 0
		chromIdx.put("*", UNMAPPED_INDEX);
		
		writerFactory.setUseAsyncIo(false);
		
		writers = new SAMFileWriter[outputFiles.length][];
		
		for (int i=0; i<writers.length; i++) {
			writers[i] = new SAMFileWriter[sequences.size()+1];
		}
	}
	
	private String getTempFilename(int sampleIdx, int chrom) {
		return String.format("%s/%d.%d.bam", tempDir, sampleIdx, chrom);
	}
	
	public void addAlignment(int sampleIdx, SAMRecord read) {
		int chrom  = chromIdx.get(read.getReferenceName());
		writers[sampleIdx][chrom].addAlignment(read);
	}
	
	public void initChromosome(int sampleIdx, String chromosome) {
		int chrom  = chromIdx.get(chromosome);
		writers[sampleIdx][chrom] = writerFactory.makeBAMWriter(samHeaders[sampleIdx], false, new File(getTempFilename(sampleIdx, chrom)), TEMP_COMPRESSION_LEVEL);
	}
	
	public void finishChromsome(int sampleIdx, String chromosome) {
		int chrom  = chromIdx.get(chromosome);
		writers[sampleIdx][chrom].close();
	}
	
	public void outputFinal() {
		for (int i=0; i<outputFiles.length; i++) {
			writerFactory.setUseAsyncIo(true);
			SAMFileWriter output = writerFactory.makeBAMWriter(samHeaders[i], false, new File(outputFiles[i]));

			for (int chrom=1; chrom<chromIdx.size(); chrom++) {
				processChromosome(output, i, chrom);
			}
			
			processUnmapped(output, i);
			
			output.close();
		}
	}
	
	private void processChromosome(SAMFileWriter output, int sampleIdx, int chrom) {
		
		List<SAMRecord> reads = new ArrayList<SAMRecord>();
		
		SAMFileReader reader = new SAMFileReader(new File(getTempFilename(sampleIdx, chrom)));
		reader.setValidationStringency(ValidationStringency.SILENT);

		for (SAMRecord read : reader) {
			reads.add(read);
			
			if (read.getAlignmentStart() - reads.get(0).getAlignmentStart() > GENOMIC_RANGE_TO_CACHE*2) {
				Collections.sort(reads, new SAMCoordinateComparator());
				reads.sort(new SAMCoordinateComparator());
				
				int start = reads.get(0).getAlignmentStart();
				int i = 0;
				while (i < reads.size() && reads.get(i).getAlignmentStart() - start < GENOMIC_RANGE_TO_CACHE) {
					output.addAlignment(reads.get(i));
				}
				reads.subList(0, i).clear();
			}
		}
		
		reader.close();
	}
	
	private void processUnmapped(SAMFileWriter output, int sampleIdx) {
		SAMFileReader reader = new SAMFileReader(new File(getTempFilename(sampleIdx, UNMAPPED_INDEX)));
		reader.setValidationStringency(ValidationStringency.SILENT);

		for (SAMRecord read : reader) {
			output.addAlignment(read);
		}
		
		reader.close();
	}
	
	static class SAMCoordinateComparator implements Comparator<SAMRecord> {

		@Override
		public int compare(SAMRecord o1, SAMRecord o2) {
			return o1.getAlignmentStart()-o2.getAlignmentStart();
		}
	}
	
	public static void main(String[] args) {
		String in = args[0];
		String out = args[1];
		String tempDir = args[2];
		
		SAMFileReader reader = new SAMFileReader(new File(in));
		reader.setValidationStringency(ValidationStringency.SILENT);
		
		SortedSAMWriter writer = new SortedSAMWriter(new String[] { out }, tempDir, new SAMFileHeader[] { reader.getFileHeader() });
		
		Set<String> chromosomes = new HashSet<String>();
		
		for (SAMRecord read : reader) {
			if (!chromosomes.contains(read.getReferenceName())) {
				writer.initChromosome(0, read.getReferenceName());
			}
			
			writer.addAlignment(0, read);
		}
		
		for (String chromosome : chromosomes) {
			writer.finishChromsome(0, chromosome);
		}
	}
}

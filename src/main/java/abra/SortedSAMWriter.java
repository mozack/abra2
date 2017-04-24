package abra;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import com.intel.gkl.compression.IntelDeflaterFactory;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;


public class SortedSAMWriter {
	
	private static final int TEMP_COMPRESSION_LEVEL = 1;
	private static final int FINAL_COMPRESSION_LEVEL = 5;
	
	public static final int GENOMIC_RANGE_TO_CACHE = 1000000;
	private static int UNMAPPED_INDEX = -1;
	private static final int ASYNC_READ_CACHE_SIZE = 1000000;
	
	private SAMFileWriterFactory writerFactory = new SAMFileWriterFactory();
	
	private SAMFileWriter writers[][];
	private String tempDir;
	private String[] outputFiles;
	private SAMFileHeader[] samHeaders;
	private boolean isKeepTmp;
	private ChromosomeChunker chromosomeChunker;
	
	private Set<Integer> chunksReady = new HashSet<Integer>();
	
	public SortedSAMWriter(String[] outputFiles, String tempDir, SAMFileHeader[] samHeaders,
			boolean isKeepTmp, ChromosomeChunker chromosomeChunker) {
	
		this.samHeaders = samHeaders;
		this.outputFiles = outputFiles;
		this.tempDir = tempDir;
		this.isKeepTmp = isKeepTmp;
		this.chromosomeChunker = chromosomeChunker;
		
		UNMAPPED_INDEX = chromosomeChunker.getChunks().size();
		
		writerFactory.setUseAsyncIo(false);
		IntelDeflaterFactory intelDeflater = new IntelDeflaterFactory();
		writerFactory.setDeflaterFactory(intelDeflater);
		
		Logger.info("Using intel deflator: " + intelDeflater.usingIntelDeflater());
		
		writers = new SAMFileWriter[outputFiles.length][];
		
		for (int i=0; i<writers.length; i++) {
			//writers[i] = new SAMFileWriter[sequences.size()+1];
			writers[i] = new SAMFileWriter[chromosomeChunker.getChunks().size()+1];
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
	
	public void addAlignment(int sampleIdx, SAMRecord read, int chromosomeChunkIdx) {
		Feature chunk = this.chromosomeChunker.getChunks().get(chromosomeChunkIdx);
		
		// Only output reads with start pos within specified chromosomeChunk
		// Avoids reads being written in 2 different chunks
		if (read.getAlignmentStart() >= chunk.getStart() && read.getAlignmentStart() <= chunk.getEnd()) {
			writers[sampleIdx][chromosomeChunkIdx].addAlignment(read);
		}
	}
	
	public void initChromosomeChunk(int chromosomeChunkIdx) {
		for (int i=0; i<outputFiles.length; i++) {
			initChromosomeChunk(i, chromosomeChunkIdx);
		}
	}
	
	private void initChromosomeChunk(int sampleIdx, int chromosomeChunkIdx) {
		Logger.debug("Writer init: %d, %d", sampleIdx, chromosomeChunkIdx);
		//int chrom  = chromIdx.get(chromosome);
		samHeaders[sampleIdx].setSortOrder(SortOrder.unsorted);
		writers[sampleIdx][chromosomeChunkIdx] = writerFactory.makeBAMWriter(samHeaders[sampleIdx], false, new File(getTempFilename(sampleIdx, chromosomeChunkIdx)), TEMP_COMPRESSION_LEVEL);
	}
	
	private void finishChromosomeChunk(int sampleIdx, int chromosomeChunkIdx) throws IOException {
		Logger.debug("Writer finish: %d, %d", sampleIdx, chromosomeChunkIdx);
		//int chrom  = chromIdx.get(chromosome);
		writers[sampleIdx][chromosomeChunkIdx].close();
	}
	
	public void finishChromosomeChunk(int chromosomeChunkIdx) throws IOException {
		for (int i=0; i<outputFiles.length; i++) {
			finishChromosomeChunk(i, chromosomeChunkIdx);
		}
		
		chunksReady.add(chromosomeChunkIdx);
//		outputFinalizedChromosomes();
	}
	
	public void outputFinal(int sampleIdx, String inputBam) throws IOException {
		
		Logger.info("Finishing: " + outputFiles[sampleIdx]);

		writerFactory.setUseAsyncIo(true);
		writerFactory.setAsyncOutputBufferSize(ASYNC_READ_CACHE_SIZE);
		writerFactory.setCompressionLevel(FINAL_COMPRESSION_LEVEL);
		samHeaders[sampleIdx].setSortOrder(SortOrder.coordinate);
		SAMFileWriter output = writerFactory.makeBAMWriter(samHeaders[sampleIdx], true, new File(outputFiles[sampleIdx]), FINAL_COMPRESSION_LEVEL);
		
		for (String chromosome : chromosomeChunker.getChromosomes()) {
			processChromosome(output, sampleIdx, chromosome);
		}
		
		processUnmapped(output, inputBam);
		
		output.close();
	}
	
	private void processChromosome(SAMFileWriter output, int sampleIdx, String chromosome) throws IOException {
		
		Logger.debug("Final processing for: %d, %s", sampleIdx, chromosome);
		
		List<SAMRecord> reads = new ArrayList<SAMRecord>();
		List<Integer> chunks = chromosomeChunker.getChunkGroups().get(chromosome);
		
		for (int chunk : chunks) {
			Logger.debug("Outputting chunk: %d", chunk);
			String filename = getTempFilename(sampleIdx, chunk);
			
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
						Logger.debug("Reads output: %d", i);
						reads.subList(0, i).clear();
					}
				}
				
				reader.close();
			}
		}
		
		// Output any remaining reads for the current chromosome
		Collections.sort(reads, new SAMCoordinateComparator());
		int i = 0;
		for (SAMRecord read : reads) {
			output.addAlignment(read);
			i += 1;
		}
		Logger.debug("Final reads output: %d", i);
	}
		
	private void processUnmapped(SAMFileWriter output, String inputBam) throws IOException {
		
		Logger.debug("Processing unmapped reads...");
		
		SamReader reader = SAMRecordUtils.getSamReader(inputBam);

		// This should give us only read pairs with both ends unmapped
		Iterator<SAMRecord> iter = reader.queryUnmapped();
		
		while (iter.hasNext()) {
			SAMRecord read = iter.next();
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
	
	/*
	public static void main(String[] args) throws IOException {
		String bam = "/home/lmose/dev/abra2_dev/sort_issue3/test2.bam";
		
		SAMFileWriterFactory writerFactory = new SAMFileWriterFactory();
		writerFactory.setUseAsyncIo(false);
//			writerFactory.setAsyncOutputBufferSize(ASYNC_READ_CACHE_SIZE);
//		writerFactory.setCompressionLevel(FINAL_COMPRESSION_LEVEL);
//		samHeaders[i].setSortOrder(SortOrder.coordinate);
		

		SamReaderFactory factory = SamReaderFactory.make()
        .validationStringency(ValidationStringency.SILENT)
        .disable(Option.EAGERLY_DECODE)
        .samRecordFactory(DefaultSAMRecordFactory.getInstance());

		SamReader reader = factory.open(new File(bam));
		
		List<SAMRecord> reads = new ArrayList<SAMRecord>();
		for (SAMRecord read : reader) {
			reads.add(read);
		}

		long start = System.currentTimeMillis();
		writerFactory.setCompressionLevel(5);
		
		SAMFileWriter out = writerFactory.makeBAMWriter(reader.getFileHeader(), true, new File("/home/lmose/dev/abra2_dev/sort_issue3/test2.out.bam"),5);
		
		for (SAMRecord read : reads) {
			out.addAlignment(read);;
		}
		
		out.close();
		
		long end = System.currentTimeMillis();
		
		System.out.println("Elapsed: " + (end-start));
	}
	*/
	
//	public static void main(String[] args) throws IOException {
////		String in = args[0];
////		String out = args[1];
////		String tempDir = args[2];
//		
////		String in = "/home/lmose/dev/abra2_dev/sort/0.1.bam";
////		String out = "/home/lmose/dev/abra2_dev/sort/output.bam";
////		String tempDir = "/home/lmose/dev/abra2_dev/sort";
//		
//		Logger.setLevel("trace");
//		
//		String ref = "/home/lmose/dev/reference/hg38/chr1.fa";
//		
//		CompareToReference2 c2r = new CompareToReference2();
//		c2r.init(ref);
//		
//		ChromosomeChunker cc = new ChromosomeChunker(c2r);
//		
//		cc.init();
//		
//		SamReader reader = SAMRecordUtils.getSamReader("/home/lmose/dev/abra2_dev/sort_issue3/0.5.bam");
//		
////		SortedSAMWriter writer = new SortedSAMWriter(new String[] { "/home/lmose/dev/abra2_dev/sort/output.bam" }, "/home/lmose/dev/abra2_dev/sort", new SAMFileHeader[] { reader.getFileHeader() });
//		
//		SortedSAMWriter writer = new SortedSAMWriter(new String[] { "/home/lmose/dev/abra2_dev/sort_issue3/final.bam" }, "/home/lmose/dev/abra2_dev/sort_issue3", new SAMFileHeader[] { reader.getFileHeader() }, true, cc);
//		
//		reader.close();
//
//		/*
//		Set<String> chromosomes = new HashSet<String>();
//		
//		for (SAMRecord read : reader) {
//			if (!chromosomes.contains(read.getReferenceName())) {
//				writer.initChromosome(0, read.getReferenceName());
//				chromosomes.add(read.getReferenceName());
//			}
//			
//			writer.addAlignment(0, read);
//		}
//		
//		for (String chromosome : chromosomes) {
//			writer.finishChromsome(0, chromosome);
//		}
//		*/
//		
//		long start = System.currentTimeMillis();
//		writer.outputFinal();
//		long stop = System.currentTimeMillis();
//	}
}

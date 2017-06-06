package abra;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
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
	
	private static final int ASYNC_READ_CACHE_SIZE = 100000;
	
	private SAMFileWriterFactory writerFactory = new SAMFileWriterFactory();
	
	private SAMFileWriter writers[][];
	private String tempDir;
	private String[] outputFiles;
	private SAMFileHeader[] samHeaders;
	private boolean isKeepTmp;
	private ChromosomeChunker chromosomeChunker;
	private int finalCompressionLevel;
	private boolean shouldSort;
	private int genomicRangeToCache;
	
	private Set<Integer> chunksReady = new HashSet<Integer>();
	
	private ReverseComplementor rc = new ReverseComplementor();
	
	public SortedSAMWriter(String[] outputFiles, String tempDir, SAMFileHeader[] samHeaders,
			boolean isKeepTmp, ChromosomeChunker chromosomeChunker, int finalCompressionLevel,
			boolean shouldSort, int genomicRangeToCache) {
	
		this.samHeaders = samHeaders;
		this.outputFiles = outputFiles;
		this.tempDir = tempDir;
		this.isKeepTmp = isKeepTmp;
		this.chromosomeChunker = chromosomeChunker;
		this.finalCompressionLevel = finalCompressionLevel;
		this.shouldSort = shouldSort;
		this.genomicRangeToCache = genomicRangeToCache;
		
		writerFactory.setUseAsyncIo(false);
		IntelDeflaterFactory intelDeflater = new IntelDeflaterFactory();
		writerFactory.setDeflaterFactory(intelDeflater);
		
		Logger.info("Using intel deflator: " + intelDeflater.usingIntelDeflater());
		
		writers = new SAMFileWriter[outputFiles.length][];
		
		for (int i=0; i<writers.length; i++) {
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
	
	public void addAlignment(int sampleIdx, SAMRecordWrapper samRecord, int chromosomeChunkIdx) {
		Feature chunk = this.chromosomeChunker.getChunks().get(chromosomeChunkIdx);
		
		SAMRecord read = samRecord.getSamRecord();
		
		// Only output reads with start pos within specified chromosomeChunk
		// Avoids reads being written in 2 different chunks
		if (read.getAlignmentStart() >= chunk.getStart() && read.getAlignmentStart() <= chunk.getEnd()) {
			
			if (samRecord.isUnalignedRc() && read.getReadUnmappedFlag()) {
				// This read was reverse complemented, but not updated.
				// Change it back to its original state.
				read.setReadString(rc.reverseComplement(read.getReadString()));
				read.setBaseQualityString(rc.reverse(read.getBaseQualityString()));
				read.setReadNegativeStrandFlag(!read.getReadNegativeStrandFlag());
			}
			
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
		samHeaders[sampleIdx].setSortOrder(SortOrder.unsorted);
		writers[sampleIdx][chromosomeChunkIdx] = writerFactory.makeBAMWriter(samHeaders[sampleIdx], false, new File(getTempFilename(sampleIdx, chromosomeChunkIdx)), TEMP_COMPRESSION_LEVEL);
	}
	
	private void finishChromosomeChunk(int sampleIdx, int chromosomeChunkIdx) throws IOException {
		Logger.debug("Writer finish: %d, %d", sampleIdx, chromosomeChunkIdx);
		writers[sampleIdx][chromosomeChunkIdx].close();
	}
	
	public void finishChromosomeChunk(int chromosomeChunkIdx) throws IOException {
		for (int i=0; i<outputFiles.length; i++) {
			finishChromosomeChunk(i, chromosomeChunkIdx);
		}
		
		chunksReady.add(chromosomeChunkIdx);
	}
	
	public void outputFinal(int sampleIdx, String inputBam) throws IOException {
		
		Logger.info("Finishing: " + outputFiles[sampleIdx]);

		writerFactory.setUseAsyncIo(true);
		writerFactory.setAsyncOutputBufferSize(ASYNC_READ_CACHE_SIZE);
		writerFactory.setCompressionLevel(finalCompressionLevel);
		if (shouldSort) {
			samHeaders[sampleIdx].setSortOrder(SortOrder.coordinate);
		} else {
			samHeaders[sampleIdx].setSortOrder(SortOrder.unsorted);
		}
		SAMFileWriter output = writerFactory.makeBAMWriter(samHeaders[sampleIdx], true, new File(outputFiles[sampleIdx]), finalCompressionLevel);
		
		for (String chromosome : chromosomeChunker.getChromosomes()) {
			processChromosome(output, sampleIdx, chromosome);
		}
		
		processUnmapped(output, inputBam);
		
		output.close();
	}
	
	private void setMateInfo(SAMRecord read, Map<MateKey, SAMRecord> mates) {
		SAMRecord mate = mates.get(getMateKey(read));
		if (mate != null) {
			// Only update mate info if a read has been modified
			if (read.getAttribute("YO") != null || mate.getAttribute("YO") != null) {
				read.setMateAlignmentStart(mate.getAlignmentStart());
				read.setMateUnmappedFlag(mate.getReadUnmappedFlag());
				read.setMateNegativeStrandFlag(mate.getReadNegativeStrandFlag());
				 
				int start = read.getAlignmentStart() < mate.getAlignmentStart() ? read.getAlignmentStart() : mate.getAlignmentStart();
				int stop  = read.getAlignmentEnd() > mate.getAlignmentEnd() ? read.getAlignmentEnd() : mate.getAlignmentEnd();
				
				int insert = stop-start+1;
				
				if (read.getAlignmentStart() > mate.getAlignmentStart()) {
					insert *= -1;
				} else if (read.getAlignmentStart() == mate.getAlignmentStart() && mate.getFirstOfPairFlag()) {
					insert *= -1;
				}
				
				read.setInferredInsertSize(insert);
			}
		}
	}
	
	private void processChromosome(SAMFileWriter output, int sampleIdx, String chromosome) throws IOException {
		
		Logger.debug("Final processing for: %d, %s", sampleIdx, chromosome);
		
		List<SAMRecord> reads = new ArrayList<SAMRecord>();
		List<Integer> chunks = chromosomeChunker.getChunkGroups().get(chromosome);
		
		// Cache read by name and original position for pairing with mate
		Map<MateKey, SAMRecord> mates = new HashMap<MateKey, SAMRecord>();
		
		for (int chunk : chunks) {
			Logger.debug("Outputting chunk: %d", chunk);
			String filename = getTempFilename(sampleIdx, chunk);
			
			File file = new File(filename);
			
			if (file.exists()) {
				deleteOnExit(file);
				
				SamReader reader = SAMRecordUtils.getSamReader(filename);
		
				for (SAMRecord read : reader) {
					if (shouldSort) {
						reads.add(read);
	
						MateKey mateKey = getOriginalReadInfo(read);
						SAMRecord existingMate = mates.get(mateKey);
						if (existingMate == null || (existingMate.getFlags() & 0xA00) != 0) {
							// Cache read info giving priority to primary alignments
							mates.put(mateKey, read);
						}
						
						if (read.getAlignmentStart() - reads.get(0).getAlignmentStart() > genomicRangeToCache*3) {
							Collections.sort(reads, new SAMCoordinateComparator());
							
							int start = reads.get(0).getAlignmentStart();
							int i = 0;
							int last = -1;
							while (i < reads.size() && reads.get(i).getAlignmentStart() - start < genomicRangeToCache) {
								SAMRecord currRead = reads.get(i);
								setMateInfo(currRead, mates);
								output.addAlignment(currRead);
								last = currRead.getAlignmentStart();
								i += 1;
							}
							
							// Clear out of scope mate keys
							List<MateKey> toRemove = new ArrayList<MateKey>();
							for (MateKey key : mates.keySet()) {
								if (key.start - start < -genomicRangeToCache) {
									toRemove.add(key);
								}
							}
							
							for (MateKey key : toRemove) {
								mates.remove(key);
							}
							
							Logger.debug("%s - Reads output: %d @ %d - %d, curr: %d", chromosome, i, start, last, read.getAlignmentStart());
							reads.subList(0, i).clear();
						}
					} else {
						output.addAlignment(read);
					}
				}
				
				reader.close();
			}
		}
		
		// Output any remaining reads for the current chromosome
		Collections.sort(reads, new SAMCoordinateComparator());
		int i = 0;
		for (SAMRecord read : reads) {
			setMateInfo(read, mates);
			output.addAlignment(read);
			i += 1;
		}
		Logger.debug("Final reads output: %d", i);
	}
	
	public MateKey getMateKey(SAMRecord read) {
		
		// If mate is mapped, use read flag.
		// If mate is not mapped, use opposite of this read's RC flag
		boolean isMateRevOrientation = read.getMateUnmappedFlag() ? !read.getReadNegativeStrandFlag() : read.getMateNegativeStrandFlag();
		int matePos = read.getMateUnmappedFlag() ? -1 : read.getMateAlignmentStart();
		
		int mateNum = read.getFirstOfPairFlag() ? 2 : 1;
		
		return new MateKey(read.getReadName(), matePos,
				read.getMateUnmappedFlag(), isMateRevOrientation, mateNum, read.getAlignmentStart());
	}
	
	public MateKey getOriginalReadInfo(SAMRecord read) {
		int pos = read.getAlignmentStart();
		boolean isUnmapped = read.getReadUnmappedFlag();
		boolean isRc = read.getReadNegativeStrandFlag();
		
		String yo = read.getStringAttribute("YO");
		
		if (yo != null) {
			if (yo.equals("N/A")) {
				// Original alignment was unmapped
				isUnmapped = true;
//				isRc = false;
				// Orientation is forced to be opposite of mate during realignment
				// regardless of the original alignment.
				isRc = !read.getMateNegativeStrandFlag();
			} else {
				String[] fields = yo.split(":");
				pos = Integer.parseInt(fields[1]);
				isUnmapped = false;
				isRc = fields[2].equals("-") ? true : false;
			}
		}
		
		int readNum = read.getFirstOfPairFlag() ? 1 : 2;
		
		return new MateKey(read.getReadName(), pos, isUnmapped, isRc, readNum, read.getAlignmentStart());
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
	
	static class MateKey {
		String readId;
		int pos;
		boolean isUnmapped;
		boolean isRc;
		int readNum;  // 1st or 2nd read
		int start; // Used for cache clearing.  Not part of identity
		
		MateKey(String readId, int pos, boolean isUnmapped, boolean isRc, int readNum, int start) {
			this.readId = readId;
			if (isUnmapped) {
				this.pos = -1;
			} else {
				this.pos = pos;
			}
			this.isUnmapped = isUnmapped;
			this.isRc = isRc;
			this.readNum = readNum;
			this.start = start;
		}

		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result + (isRc ? 1231 : 1237);
			result = prime * result + (isUnmapped ? 1231 : 1237);
			result = prime * result + pos;
			result = prime * result + ((readId == null) ? 0 : readId.hashCode());
			result = prime * result + readNum;
			return result;
		}

		@Override
		public boolean equals(Object obj) {
			if (this == obj)
				return true;
			if (obj == null)
				return false;
			if (getClass() != obj.getClass())
				return false;
			MateKey other = (MateKey) obj;
			if (isRc != other.isRc)
				return false;
			if (isUnmapped != other.isUnmapped)
				return false;
			if (pos != other.pos)
				return false;
			if (readId == null) {
				if (other.readId != null)
					return false;
			} else if (!readId.equals(other.readId))
				return false;
			if (readNum != other.readNum)
				return false;
			return true;
		}

		@Override
		public String toString() {
			return "MateKey [readId=" + readId + ", pos=" + pos + ", isUnmapped=" + isUnmapped + ", isRc=" + isRc
					+ ", readNum=" + readNum + ", start=" + start + "]";
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
	
	public static void main(String[] args) throws IOException {
		Logger.LEVEL = Logger.Level.TRACE;
		String ref = "/home/lmose/dev/reference/hg38/chr1.fa";
		
		CompareToReference2 c2r = new CompareToReference2();
		c2r.init(ref);
		
		ChromosomeChunker cc = new ChromosomeChunker(c2r);
		
		cc.init();
		
		SamReader reader = SAMRecordUtils.getSamReader("/home/lmose/dev/abra2_dev/sort_issue3/0.5.bam");
		SAMFileHeader header = reader.getFileHeader();
		header.setSortOrder(SortOrder.coordinate);
		
		SortedSAMWriter writer = new SortedSAMWriter(new String[] { "/home/lmose/dev/abra2_dev/sort_issue4/final.bam" }, "/home/lmose/dev/abra2_dev/sort_issue4", new SAMFileHeader[] { reader.getFileHeader() }, true, cc,
				1,true,1000);

		SAMFileWriterFactory writerFactory = new SAMFileWriterFactory();
		SAMFileWriter out = writerFactory.makeBAMWriter(reader.getFileHeader(), true, new File("/home/lmose/dev/abra2_dev/sort_issue4/test.bam"),1);
		
		writer.processChromosome(out, 1, "chr1");
		
		out.close();
	}
	
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

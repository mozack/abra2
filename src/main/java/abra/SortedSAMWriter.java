package abra;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
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
	private boolean shouldUnsetDuplicates;
	private boolean shouldCreateIndex;
	
	private int maxRecordsInRam;
	
	private Set<Integer> chunksReady = new HashSet<Integer>();
	
	private ReverseComplementor rc = new ReverseComplementor();
	
	public SortedSAMWriter(String[] outputFiles, String tempDir, SAMFileHeader[] samHeaders,
			boolean isKeepTmp, ChromosomeChunker chromosomeChunker, int finalCompressionLevel,
			boolean shouldSort, int genomicRangeToCache, boolean shouldUnsetDuplicates,
			boolean shouldCreateIndex, boolean shouldUseGkl, int maxReadsInRam) {
	
		this.samHeaders = samHeaders;
		this.outputFiles = outputFiles;
		this.tempDir = tempDir;
		this.isKeepTmp = isKeepTmp;
		this.chromosomeChunker = chromosomeChunker;
		this.finalCompressionLevel = finalCompressionLevel;
		this.shouldSort = shouldSort;
		this.genomicRangeToCache = genomicRangeToCache;
		this.shouldUnsetDuplicates = shouldUnsetDuplicates;
		this.shouldCreateIndex = shouldCreateIndex;

		this.maxRecordsInRam = maxReadsInRam;

		if (shouldUseGkl) {
			writerFactory.setUseAsyncIo(false);
			IntelDeflaterFactory intelDeflater = new IntelDeflaterFactory();
			writerFactory.setDeflaterFactory(intelDeflater);
			
			Logger.info("Using intel deflator: " + intelDeflater.usingIntelDeflater());
		} else {
			Logger.info("Intel deflater disabled");
		}
		
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
		
		// Only output reads with original start pos within specified chromosomeChunk
		// Avoids reads being written in 2 different chunks
		
		int origAlignmentStart = read.getAlignmentStart();
		String yo = read.getStringAttribute("YO");
		if (yo != null) {
			String[] fields = yo.split(":");
			origAlignmentStart = Integer.parseInt(fields[1]);
		}
		
		if (origAlignmentStart >= chunk.getStart() && origAlignmentStart <= chunk.getEnd()) {
			
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
		
		SAMRecord[] readsByNameArray = null;
		SAMRecord[] readsByCoordArray = null;

		if (shouldSort) {
			
			// Initialize internal read buffers used by SortingCollection2
			// These are initialized once per thread and re-used each time
			// a SortingCollection2 is initialized.
			readsByNameArray = new SAMRecord[maxRecordsInRam];
			readsByCoordArray = new SAMRecord[maxRecordsInRam];
			
			// Only allow buffering if sorting
			writerFactory.setUseAsyncIo(true);
			writerFactory.setAsyncOutputBufferSize(ASYNC_READ_CACHE_SIZE);
			writerFactory.setCreateIndex(shouldCreateIndex);
		} else {
			writerFactory.setUseAsyncIo(false);
		}
		
		writerFactory.setCompressionLevel(finalCompressionLevel);
		if (shouldSort) {
			samHeaders[sampleIdx].setSortOrder(SortOrder.coordinate);
		} else {
			samHeaders[sampleIdx].setSortOrder(SortOrder.unsorted);
		}
		SAMFileWriter output = writerFactory.makeBAMWriter(samHeaders[sampleIdx], true, new File(outputFiles[sampleIdx]), finalCompressionLevel);
		
		for (String chromosome : chromosomeChunker.getChromosomes()) {
			processChromosome(output, sampleIdx, chromosome, readsByNameArray, readsByCoordArray);
		}
		
		processUnmapped(output, inputBam);
		
		output.close();
	}
	
	private void setMateInfo(SAMRecord read, Map<MateKey, SAMRecord> mates) {
		
		if (read.getReadPairedFlag()) {
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
	}
	
	private SortingSAMRecordCollection updateReadMatesAndSortByCoordinate(SortingSAMRecordCollection readsByName, SAMFileHeader header, SAMRecord[] readsByCoordArray) {
		
		SortingSAMRecordCollection readsByCoord = SortingSAMRecordCollection.newSortByCoordinateInstance(readsByCoordArray, header, maxRecordsInRam, tempDir); 
		
		Iterator<SAMRecord> iter = readsByName.iterator();
		
		// Cache read by name and original position for pairing with mate
		Map<MateKey, SAMRecord> mates = new HashMap<MateKey, SAMRecord>();
		
		String currReadName = "";
		List<SAMRecord> currReads = new ArrayList<SAMRecord>();
		
		while (iter.hasNext()) {
			SAMRecord read = iter.next();
			
			if (currReadName.isEmpty()) {
				currReadName = read.getReadName();
			}
			else if (!read.getReadName().equals(currReadName)) {
				// Go through all reads with same read name and update mates
				for (SAMRecord currRead : currReads) {
					setMateInfo(currRead, mates);
					
					// Final set of output reads are sorted by coordinate
					readsByCoord.add(currRead);
				}
				
				mates.clear();
				currReads.clear();
				currReadName = read.getReadName();
			}
			
			// Cache reads with same read name
			currReads.add(read);
			
			// Cache read by mate info
			if (read.getReadPairedFlag() && read.getSupplementaryAlignmentFlag() != true && (read.getFlags() & 0x100) == 0) {
				MateKey mateKey = getOriginalReadInfo(read);
				SAMRecord existingMate = mates.get(mateKey);
				
	//			if (existingMate == null || existingMate.getSupplementaryAlignmentFlag() == true || (existingMate.getFlags() & 0x100) != 0) {
				if (existingMate == null) {
					// Cache read info giving priority to primary alignments
					mates.put(mateKey, read);
				}
			}
		}
		
		// Update stragglers
		for (SAMRecord currRead : currReads) {
			setMateInfo(currRead, mates);
			
			// Final set of output reads are sorted by coordinate
			readsByCoord.add(currRead);
		}
		
		return readsByCoord;
	}
	
	private void processChromosome(SAMFileWriter output, int sampleIdx, String chromosome,
			SAMRecord[] readsByNameArray, SAMRecord[] readsByCoordArray) throws IOException {
		
		Logger.debug("Final processing for: %d, %s", sampleIdx, chromosome);
		
		SAMFileHeader sortByCoordHeader = output.getFileHeader();
		sortByCoordHeader.setSortOrder(SortOrder.coordinate);
		
		SAMFileHeader sortByNameHeader = output.getFileHeader();
		sortByNameHeader.setSortOrder(SortOrder.queryname);
		
		List<Integer> chunks = chromosomeChunker.getChunkGroups().get(chromosome);
		
		SortingSAMRecordCollection readsByName = SortingSAMRecordCollection.newSortByNameInstance(readsByNameArray, sortByNameHeader, maxRecordsInRam, tempDir);
		
		for (int chunk : chunks) {
			Logger.debug("Outputting chunk: %d", chunk);
			String filename = getTempFilename(sampleIdx, chunk);
			
			File file = new File(filename);
			
			if (file.exists()) {
				deleteOnExit(file);
				
				SamReader reader = SAMRecordUtils.getSamReader(filename);

				int firstReadPos = -1;
		
				for (SAMRecord read : reader) {
					if (shouldUnsetDuplicates) {
						read.setDuplicateReadFlag(false);
					}
					
					if (shouldSort) {
						readsByName.add(read);
						
						if (firstReadPos < 0) {
							firstReadPos = read.getAlignmentStart();
						}
						
						// Require 4x maxDist reads to be cached before flushing the first 1x maxDist reads
						if (firstReadPos >= 0 && read.getAlignmentStart() - firstReadPos > genomicRangeToCache*4) {
							
							SortingSAMRecordCollection readsByCoord = updateReadMatesAndSortByCoordinate(readsByName, sortByCoordHeader, readsByCoordArray);
							
							// Reset readsByName collection
							readsByName.cleanup();
							readsByName = SortingSAMRecordCollection.newSortByNameInstance(readsByNameArray, sortByNameHeader, maxRecordsInRam, tempDir);
							
							int start = firstReadPos;
							int i = 0;
							int last = -1;
							firstReadPos = -1;
							
							Iterator<SAMRecord> iter = readsByCoord.iterator();
							while (iter.hasNext()) {
								SAMRecord sortedRead = iter.next();
								
								if (sortedRead.getAlignmentStart() - start < genomicRangeToCache) {
									output.addAlignment(sortedRead);
									last = sortedRead.getAlignmentStart();
									i += 1;
								} else {
									if (firstReadPos == -1) {
										firstReadPos = sortedRead.getAlignmentStart();
									}
									readsByName.add(sortedRead);
								}
							}
							
							readsByCoord.cleanup();
							
							Logger.debug("%s - Reads output: %d @ %d - %d, curr: %d", chromosome, i, start, last, read.getAlignmentStart());
						}
					} else {
						output.addAlignment(read);
					}
				}
				
				reader.close();
			}
		}
		
		// Output any remaining reads for the current chromosome
		if (shouldSort) {
			SortingSAMRecordCollection readsByCoord = updateReadMatesAndSortByCoordinate(readsByName, output.getFileHeader(), readsByCoordArray);
			readsByName.cleanup();
			
			int i = 0;
			Iterator<SAMRecord> iter = readsByCoord.iterator();
			while (iter.hasNext()) {
				SAMRecord sortedRead = iter.next();
				output.addAlignment(sortedRead);
				i += 0;
			}
			Logger.debug("Final reads output: %d", i);
		}
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
			if (yo.startsWith("N/A")) {
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
//		Logger.LEVEL = Logger.Level.TRACE;
		String ref = "/home/lmose/dev/reference/hg38b/hg38.fa";
		
		CompareToReference2 c2r = new CompareToReference2();
		c2r.init(ref);
		
		ChromosomeChunker cc = new ChromosomeChunker(c2r);
		
		cc.init();
		
//		SamReader reader = SAMRecordUtils.getSamReader("/home/lmose/dev/abra2_dev/sort_issue3/0.5.bam");
//		SamReader reader = SAMRecordUtils.getSamReader("/home/lmose/dev/abra2_dev/sort_issue4/1.6.bam");
		SamReader reader = SAMRecordUtils.getSamReader("/home/lmose/dev/abra2_dev/mate_fix/0.87.bam");

		SAMFileHeader header = reader.getFileHeader();
		header.setSortOrder(SortOrder.coordinate);
		
		int maxRecordsInRam = 1000000;
		SAMRecord[] readsByNameArray = new SAMRecord[maxRecordsInRam];
		SAMRecord[] readsByCoordArray = new SAMRecord[maxRecordsInRam];
		
		SortedSAMWriter writer = new SortedSAMWriter(new String[] { "/home/lmose/dev/abra2_dev/mate_fix" }, "/home/lmose/dev/abra2_dev/mate_fix", new SAMFileHeader[] { reader.getFileHeader() }, true, cc,
				1,true,1000,false, false, false, maxRecordsInRam);

		SAMFileWriterFactory writerFactory = new SAMFileWriterFactory();
		SAMFileWriter out = writerFactory.makeBAMWriter(reader.getFileHeader(), true, new File("/home/lmose/dev/abra2_dev/mate_fix/test.bam"),1);
		
		long start = System.currentTimeMillis();
		
		writer.processChromosome(out, 0, "chr12", readsByNameArray, readsByCoordArray);
		
		out.close();
		
		long stop = System.currentTimeMillis();
		
		System.out.println("Elapsed msecs: " + (stop-start));
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

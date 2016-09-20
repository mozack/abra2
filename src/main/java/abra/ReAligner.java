/* Copyright 2013 University of North Carolina at Chapel Hill.  All rights reserved. */
package abra;

import static abra.Logger.log;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.net.InetAddress;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import abra.ReadEvaluator.Alignment;
import abra.SSWAligner.SSWAlignerResult;
import abra.SimpleMapper.Orientation;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.ValidationStringency;

/**
 * ABRA's main entry point
 * 
 * @author Lisle E. Mose (lmose at unc dot edu)
 */
public class ReAligner {

	private static final int DEFAULT_MAX_UNALIGNED_READS = 1000000;
	
	public static final int MAX_REGION_LENGTH = 400;
	private static final int MIN_REGION_REMAINDER = 200;
	public static final int REGION_OVERLAP = 200;
	
	// Minimum sequence length recommended for use with bwa mem
	private static final int MIN_CONTIG_LENGTH = 70;
	
	// Cannot be larger than buffer in assembler.c
	private static final int MAX_KMER_SIZE = 199;
	
	private SAMFileHeader[] samHeaders;
	
	private List<Feature> regions;

	private String regionsBed;

	private String tempDir;
	
	private String unalignedRegionSam;

	private String reference;
	
	private String bwaIndex;

	private AssemblerSettings assemblerSettings;
	
	private int numThreads;
	
	private int maxUnalignedReads = DEFAULT_MAX_UNALIGNED_READS;
	
	private boolean shouldReprocessUnaligned = true;
	
	private String localRepeatFile;
	
	private String[] inputSams;
	private SAMFileWriter[] writers;
	
	private int readLength = -1;
	private int maxMapq = -1;
	private int minInsertLength = Integer.MAX_VALUE;
	private int maxInsertLength = -1;
	
	private boolean isPairedEnd = false;
	
	private BufferedWriter contigWriter;
	
	private CompareToReference2 c2r;
	
	private ThreadManager threadManager;
	
	private int minMappingQuality;
	
	private boolean isDebug;
	
	// If true, the input target file specifies kmer values
	private boolean hasPresetKmers = false;
	
	public static final int COMPRESSION_LEVEL = 1;
	
	public void reAlign(String[] inputFiles, String[] outputFiles) throws Exception {
		
		this.inputSams = inputFiles;
		
		logStartupInfo(outputFiles);
				
		init();
		
		c2r = new CompareToReference2();
		c2r.init(this.reference);

		log("Reading Input SAM Header and identifying read length");
		getSamHeaderAndReadLength();
		
		log("Read length: " + readLength);
		
		log("Loading target regions");
		loadRegions();
		
		Clock clock = new Clock("Assembly");
		clock.start();
		
		String contigFasta = tempDir + "/" + "all_contigs.fasta";
		contigWriter = new BufferedWriter(new FileWriter(contigFasta, false));
		
		SAMFileWriterFactory writerFactory = new SAMFileWriterFactory();
		
//		writerFactory.setUseAsyncIo(true);
//		writerFactory.setAsyncOutputBufferSize(500000);
		
		writerFactory.setUseAsyncIo(false);
		
		writers = new SAMFileWriter[inputSams.length];
		
		for (int i=0; i<inputSams.length; i++) {
			// init BAM writer
			writers[i] = writerFactory.makeBAMWriter(
					samHeaders[i], false, new File(outputFiles[i]), COMPRESSION_LEVEL);
		}

		int currRegionIdx = -1;
		
		MultiSamReader reader = new MultiSamReader(this.inputSams, this.minMappingQuality, this.isPairedEnd);
		
		List<List<SAMRecordWrapper>> currReads = new ArrayList<List<SAMRecordWrapper>>();
		for (int i=0; i<this.inputSams.length; i++) {
			currReads.add(new ArrayList<SAMRecordWrapper>());
		}
		
		Map<Feature, Map<SimpleMapper, SSWAlignerResult>> regionContigs = new HashMap<Feature, Map<SimpleMapper, SSWAlignerResult>>();
		int readCount = 0;
		
		// TODO: Consider parallelizing by chromosome
		for (SAMRecordWrapper record : reader) {
			int regionIdx = Feature.findFirstOverlappingRegion(reader.getSAMFileHeader(), record.getSamRecord(), regions, currRegionIdx);
			
			System.err.println("currRegionIdx: " + currRegionIdx);
			System.err.println("regionIdx: " + regionIdx);
			
			if ((regionIdx == -1 && currRegionIdx >= 0) ||
				(regionIdx > currRegionIdx)) {
				
				if (currRegionIdx >= 0) {
					
					// We've moved beyond the current region
					// Assemble reads
					Feature region = regions.get(currRegionIdx);
					System.err.println("Processing region: " + region);
					Map<SimpleMapper, SSWAlignerResult> mappedContigs = processRegion(region, currReads);
					System.err.println("Region: " + region + " assembled: " + mappedContigs.keySet().size() + " contigs");
					regionContigs.put(region, mappedContigs);
				}
				
				currRegionIdx = regionIdx;
			}
			
			if (regionIdx >= 0) {
				// Read is in a target region
				if (currRegionIdx == regionIdx) {
					// Read is in the currently tracked region.  Cache it for processing at end of region
					currReads.get(record.getSampleIdx()).add(record);
				} else {
					throw new IllegalStateException("Out of order reads / regions for sample: " + record.getSampleIdx() + 
							", read: " + record.getSamRecord().getSAMString() + ", region: " + regions.get(regionIdx) + 
							", curr_region: " + regions.get(currRegionIdx));
					
				}
			}
			
			// Todo - make constant or parameterize
			int MAX_READ_RANGE = 1000 + this.readLength;
			
			// Check for out of scope reads every 2500 reads (TODO: is 2500 the best number?)
			if (readCount % 2500 == 0) {
				// Remap / output / clear out of scope reads
				List<List<SAMRecordWrapper>> readsToRemap = new ArrayList<List<SAMRecordWrapper>>();
				
				// Initialize per sample lists
				for (List<SAMRecordWrapper> origSample : currReads) {
					List<SAMRecordWrapper> sampleReadsToRemap = new ArrayList<SAMRecordWrapper>();
					readsToRemap.add(sampleReadsToRemap);
					
					Iterator<SAMRecordWrapper> iter = origSample.iterator();
					while (iter.hasNext()) {
						SAMRecordWrapper read = iter.next();
						if (record.getSamRecord().getAlignmentStart() - read.getSamRecord().getAlignmentStart() > MAX_READ_RANGE) {
							sampleReadsToRemap.add(read);
							iter.remove();
						}
					}					
				}

				// Remap out of scope reads
				remapReads(regionContigs, readsToRemap);
				
				// Remove out of scope region assemblies
				List<Feature> regionsToRemove = new ArrayList<Feature>();
				for (Feature region : regionContigs.keySet()) {
					if (region.getStart() - getFirstStartPos(currReads) > MAX_READ_RANGE) {
						regionsToRemove.add(region);
					}
				}
				
				for (Feature region : regionsToRemove) {
					System.err.println("Removing contigs for region: " + region);
					regionContigs.remove(region);
				}
			}
			
			readCount += 1;
		}
		
		// Process last region
		if (currRegionIdx >= 0) {
			
			// We've moved beyond the current region
			// Assemble reads
			Feature region = regions.get(currRegionIdx);
			System.err.println("Processing region: " + region);
			Map<SimpleMapper, SSWAlignerResult> mappedContigs = processRegion(region, currReads);
			System.err.println("Region: " + region + " assembled: " + mappedContigs.keySet().size() + " contigs");
			regionContigs.put(region, mappedContigs);
		}
		
		System.err.println("Remapping reads");
		// Remap remaining reads
		remapReads(regionContigs, currReads);
		currReads.clear();
		regionContigs.clear();
		
		reader.close();
		
		/*
		log("Iterating over regions");
		
		int count = 0;
		for (Feature region : regions) {
			count += 1;
			spawnRegionThread(region, null);
			if ((count % 1000) == 0) {
				System.err.println("Processing region: " + count + " of " + regions.size());
			}
		}
		*/
		
		log("Waiting for all threads to complete");
		threadManager.waitForAllThreadsToComplete();
		
		contigWriter.close();
		
		clock.stopAndPrint();		
		
		for (SAMFileWriter writer : this.writers) {
			writer.close();
		}
		
		System.err.println("Done.");
	}
	
	private int getFirstStartPos(List<List<SAMRecordWrapper>> readsList) {
		int minPos = Integer.MAX_VALUE;
		for (List<SAMRecordWrapper> reads : readsList) {
			if (reads.get(0).getSamRecord().getAlignmentStart() < minPos) {
				minPos = reads.get(0).getSamRecord().getAlignmentStart(); 
			}
		}
		
		return minPos;
	}
	
	private void logStartupInfo(String[] outputFiles) {
		
		int ctr = 0;
		for (String input : inputSams) {
			System.err.println("input" + ctr + ": " + input);
		}

		ctr = 0;
		for (String output : outputFiles) {
			System.err.println("output" + ctr + ": " + output);
		}
		
		System.err.println("regions: " + regionsBed);
		System.err.println("reference: " + reference);
		System.err.println("bwa index: " + bwaIndex);
		System.err.println("working dir: " + tempDir);
		System.err.println("num threads: " + numThreads);
		System.err.println("max unaligned reads: " + maxUnalignedReads);
		System.err.println(assemblerSettings.getDescription());
		System.err.println("paired end: " + isPairedEnd);
		
		String javaVersion = System.getProperty("java.version");
		System.err.println("Java version: " + javaVersion);
		if (javaVersion.startsWith("1.6") || javaVersion.startsWith("1.5") || javaVersion.startsWith("1.4")) {
			throw new RuntimeException("Please upgrade to Java 7 or later to run ABRA.");
		}
		
		try {
			InetAddress localhost = java.net.InetAddress.getLocalHost();
			String hostname = localhost.getHostName();
			System.err.println("hostname: " + hostname);
		} catch (Throwable t) {
			System.err.println("Error getting hostname: " + t.getMessage());
		}
	}
		
	private void spawnRegionThread(Feature region, String inputSam) throws InterruptedException {
		ReAlignerRunnable thread = new ReAlignerRunnable(threadManager, this, region);
		threadManager.spawnThread(thread);
	}
	
	private synchronized void appendContigs(String contigs) throws IOException {
		contigWriter.write(contigs);
	}
	
	private void remapRead(ReadEvaluator readEvaluator, SAMRecord read, int origEditDist) {
		Alignment alignment = readEvaluator.getImprovedAlignment(origEditDist, read.getReadString(), read);
		if (alignment != null) {
			
			System.err.println("Updating: " + read);
			
			int readPos = alignment.pos;
			
			// Set contig alignment info for all reads that map to contigs (even if read is unchanged)
			String ya = alignment.chromosome + ":" + alignment.contigPos + ":" + alignment.contigCigar;
			read.setAttribute("YA", ya);
			
			// If the read has actually moved, updated it
			if (read.getReadUnmappedFlag() || read.getAlignmentStart() != readPos || !read.getCigarString().equals(alignment.cigar)) {

				// Original alignment info
				String yo = "N/A";
				if (!read.getReadUnmappedFlag()) {
					String origOrientation = read.getReadNegativeStrandFlag() ? "-" : "+";
					yo = read.getReferenceName() + ":" + read.getAlignmentStart() + ":" + origOrientation + ":" + read.getCigarString();
				} else {
					read.setReadUnmappedFlag(false);
				}
				read.setAttribute("YO", yo);

				// Update alignment position and cigar and orientation
				read.setAlignmentStart(alignment.pos);
				read.setCigarString(alignment.cigar);
				
				// If this is true, the read was already reverse complemented in the original alignment
				if (read.getReadNegativeStrandFlag()) {
					read.setReadNegativeStrandFlag(alignment.orientation == Orientation.FORWARD ? true : false);
				} else {
					read.setReadNegativeStrandFlag(alignment.orientation == Orientation.FORWARD ? false : true);
				}
				
				// Number of mismatches to contig
				read.setAttribute("YM", alignment.numMismatches);

				// Original edit distance
				read.setAttribute("YX",  origEditDist);
				
				// Updated edit distance
				read.setAttribute("NM", SAMRecordUtils.getEditDistance(read, c2r));
				
				//TODO: Compute mapq intelligently???
				read.setMappingQuality(Math.min(read.getMappingQuality(), 60));
			}
		}
	}
	
	private void remapReads(Map<Feature, Map<SimpleMapper, SSWAlignerResult>> mappedContigs, List<List<SAMRecordWrapper>> readsList) throws Exception {
		ReadEvaluator readEvaluator = new ReadEvaluator(mappedContigs);
		
		int sampleIdx = 0;
		// For each sample.
		for (List<SAMRecordWrapper> reads : readsList) {
			
			// For each read.
			for (SAMRecordWrapper readWrapper : reads) {
				SAMRecord read = readWrapper.getSamRecord();
				// TODO: Use NM tag if available (need to handle soft clipping though!)
				int origEditDist = SAMRecordUtils.getEditDistance(read, c2r);
//				int origEditDist = c2r.numMismatches(read);
				
				if (origEditDist > 0) {
					System.err.println("Remapping: " + read.getSAMString());
					remapRead(readEvaluator, read, origEditDist);
				}
			}

			// Output all reads for this sample - synchronize on the current BAM
			synchronized(this.writers[sampleIdx]) {
				for (SAMRecordWrapper read : reads) {
					this.writers[sampleIdx].addAlignment(read.getSamRecord());
				}
			}
			
			sampleIdx += 1;
		}
	}
	
	private List<List<SAMRecordWrapper>> subsetReads(Feature region, List<List<SAMRecordWrapper>> readsList) {
		List<List<SAMRecordWrapper>> subset = new ArrayList<List<SAMRecordWrapper>>();
		
		// Initialize per sample lists
		for (List<SAMRecordWrapper> origSample : readsList) {
			List<SAMRecordWrapper> subsetSample = new ArrayList<SAMRecordWrapper>();
			subset.add(subsetSample);
			
			for (SAMRecordWrapper read : origSample) {
				if (region.overlapsRead(read.getSamRecord())) {
					subsetSample.add(read);
				}
			}
		}
		
		return subset;
	}
	
	public Map<SimpleMapper, SSWAlignerResult> processRegion(Feature region, List<List<SAMRecordWrapper>> reads) throws Exception {
		if (isDebug) {
			log("Processing region: " + region.getDescriptor());
		}
		
		Map<SimpleMapper, SSWAlignerResult> mappedContigs = new HashMap<SimpleMapper, SSWAlignerResult>();
		
		List<List<SAMRecordWrapper>> readsList = subsetReads(region, reads);
		
		
		try {
			String contigsFasta = tempDir + "/" + region.getDescriptor() + "_contigs.fasta";
			
			List<String> bams = new ArrayList<String>(Arrays.asList(this.inputSams));
			
			// Assemble contigs
			if (region.getKmer() > this.readLength-15) {
				System.err.println("Skipping assembly of region: " + region.getDescriptor() + " - " + region.getKmer());
			} else {
				NativeAssembler assem = (NativeAssembler) newAssembler(region);
				List<Feature> regions = new ArrayList<Feature>();
				regions.add(region); 
//				List<List<SAMRecordWrapper>> readsList = ReadLoader.getReads(bams, regions.get(0), this);
				String contigs = assem.assembleContigs(bams, contigsFasta, tempDir, regions, region.getDescriptor(), true, this, c2r, readsList);
				
				if (!contigs.equals("<ERROR>") && !contigs.equals("<REPEAT>") && !contigs.isEmpty()) {
					
					appendContigs(contigs);
					
					// Get reference sequence matching current region (TODO: Add extra padding to allow for indels in partially overlapping reads ?)
					System.out.println("Getting reference for: " + region.getSeqname() + ":" + (region.getStart()-this.readLength) + ", length: " + (region.getLength()+this.readLength*2));
					
					int refStart = (int) region.getStart() - this.readLength;
					String refSeq = c2r.getSequence(region.getSeqname(), refStart, (int) region.getLength()+this.readLength*2);
					
					SSWAligner ssw = new SSWAligner(refSeq, region.getSeqname(), refStart);
					
					// Map contigs to reference
					String[] contigSequences = contigs.split("\n");
					for (String contig : contigSequences) {
						if (!contig.startsWith(">")) {
							SSWAlignerResult sswResult = ssw.align(contig);
							// TODO: In multi-region processing, check to ensure identical contigs have identical mappings
							mappedContigs.put(new SimpleMapper(contig), sswResult);
						}
					}
					
					// remap reads
//					remapReads(mappedContigs, readsList, region, refStart);
				}
			}
		}
		catch (Exception e) {
			e.printStackTrace();
			throw e;
		}
		
		return mappedContigs;
	}
	
	static List<Feature> getRegions(String regionsBed, int readLength, boolean hasPresetKmers) throws IOException {
		RegionLoader loader = new RegionLoader();
		List<Feature> regions = loader.load(regionsBed, hasPresetKmers);
		if (regions.size() > 0 && (regions.get(0).getKmer() == 0)) {
			regions = RegionLoader.collapseRegions(regions, readLength);
			regions = splitRegions(regions);
		}

		return regions;
	}
		
	private void loadRegions() throws IOException {
		this.regions = getRegions(regionsBed, readLength, hasPresetKmers);
		
		System.err.println("Num regions: " + regions.size());
		if (isDebug) {
			for (Feature region : regions) {
				System.err.println(region.getSeqname() + "\t" + region.getStart() + "\t" + region.getEnd() + "\t" + region.getKmer());
			}
		}
	}

	public void setRegionsBed(String bedFile) {
		this.regionsBed = bedFile;
	}

	private void getSamHeaderAndReadLength() {
		
		log("Identifying header and determining read length");
		
		this.samHeaders = new SAMFileHeader[this.inputSams.length];
		
		for (int i=0; i<this.inputSams.length; i++) {
		
			SAMFileReader reader = new SAMFileReader(new File(inputSams[i]));
			try {
				reader.setValidationStringency(ValidationStringency.SILENT);
		
				samHeaders[i] = reader.getFileHeader();
				samHeaders[i].setSortOrder(SAMFileHeader.SortOrder.unsorted);
				
				Iterator<SAMRecord> iter = reader.iterator();
				
				int cnt = 0;
				while ((iter.hasNext()) && (cnt < 1000000)) {
					SAMRecord read = iter.next();
					this.readLength = Math.max(this.readLength, read.getReadLength());
					this.maxMapq = Math.max(this.maxMapq, read.getMappingQuality());
					
					// Assumes aligner sets proper pair flag correctly
					if ((isPairedEnd) && (read.getReadPairedFlag()) && (read.getProperPairFlag())) {
						this.minInsertLength = Math.min(this.minInsertLength, Math.abs(read.getInferredInsertSize()));
						this.maxInsertLength = Math.max(this.maxInsertLength, Math.abs(read.getInferredInsertSize()));
					}
				}
				
				// Allow some fudge in insert length
				minInsertLength = Math.max(minInsertLength - 2*readLength, 0);
				maxInsertLength = maxInsertLength + 2*readLength;
				
			} finally {
				reader.close();
			}
		}
		
		System.err.println("Min insert length: " + minInsertLength);
		System.err.println("Max insert length: " + maxInsertLength);
				
		log("Max read length is: " + readLength);
		if (assemblerSettings.getMinContigLength() < 1) {
			assemblerSettings.setMinContigLength(Math.max(readLength+1, MIN_CONTIG_LENGTH));
		}
		log("Min contig length: " + assemblerSettings.getMinContigLength());
	}
			
	static class Pair<T, Y> {
		private T t;
		private Y y;
		public Pair(T t, Y y) {
			this.t = t;
			this.y = y;
		}
		
		public T getFirst() {
			return t;
		}
		
		public Y getSecond() {
			return y;
		}
	}
	
	static List<Feature> splitRegions(List<Feature> regions, 
			int maxRegionLength, int minRegionRemainder, int regionOverlap) {
		
		List<Feature> splitRegions = new ArrayList<Feature>();
		
		for (Feature region : regions) {
			if (region.getLength() <= maxRegionLength + minRegionRemainder) {
				splitRegions.add(region);
			} else {
				splitRegions.addAll(splitWithOverlap(region, maxRegionLength, minRegionRemainder, regionOverlap));
			}
		}
		
		return splitRegions;
	}
	
	/**
	 *  If any of the input list of features is greater than maxSize, split them into multiple features. 
	 */
	public static List<Feature> splitRegions(List<Feature> regions) {
		return splitRegions(regions, MAX_REGION_LENGTH, MIN_REGION_REMAINDER, REGION_OVERLAP);
	}
	
	public static List<Feature> splitWithOverlap(Feature region) {
		return splitWithOverlap(region, MAX_REGION_LENGTH, MIN_REGION_REMAINDER, REGION_OVERLAP);
	}
	
	static List<Feature> splitWithOverlap(Feature region, int maxRegionLength,
			int minRegionRemainder, int regionOverlap) {
		List<Feature> regions = new ArrayList<Feature>();
		
		long pos = region.getStart();
		long end = pos-1;
		
		while (end < region.getEnd()) {
			long start = pos;
			end = pos + maxRegionLength;
			long marker = end;
			
			// If we're at or near the end of the region, stop at region end.
			if (end > (region.getEnd() - minRegionRemainder)) {
				end = region.getEnd();
			}
			
			pos = marker - regionOverlap;
			
			regions.add(new Feature(region.getSeqname(), start, end));
		}
		
		return regions;
	}
	
	int[] getKmers(Feature region) {
		int[] kmerSizes = null;
		
		int kmerSize = region.getKmer();
		
		if (kmerSize > 0) {
			kmerSizes = toKmerArray(kmerSize, readLength);
		} else {
			kmerSizes = assemblerSettings.getKmerSize();
		}
		
		return kmerSizes;
	}
	
	int[] toKmerArray(int kmerSize, int readLength) {
		int[] kmerSizes = null;
		
		int maxKmerSize = this.readLength-15; 
		if (maxKmerSize > MAX_KMER_SIZE) {
			maxKmerSize = MAX_KMER_SIZE;
		}
		List<Integer> kmers = new ArrayList<Integer>();
		
		while (kmerSize < maxKmerSize) {
			kmers.add(kmerSize);
			kmerSize += 2;
		}
		
		kmerSizes = new int[kmers.size()];
		
		int i=0;
		for (int kmer : kmers) {
			kmerSizes[i++] = kmer;
		}

		return kmerSizes;
	}
		
	private NativeAssembler newAssembler(Feature region) {
		NativeAssembler assem = new NativeAssembler();

		assem.setTruncateOutputOnRepeat(true);
		assem.setMaxContigs(assemblerSettings
				.getMaxPotentialContigs());

		assem.setMaxPathsFromRoot(100000);
		assem.setReadLength(readLength);
		//assem.setKmer(assemblerSettings.getKmerSize());
		assem.setKmer(getKmers(region));
		assem.setMinKmerFrequency(assemblerSettings.getMinNodeFrequncy());
		assem.setMinEdgeRatio(assemblerSettings.getMinEdgeRatio());
		assem.setMinBaseQuality(assemblerSettings.getMinBaseQuality());
		assem.setMaxNodes(assemblerSettings.getMaxNodes());
		assem.setMinReadCandidateFraction(assemblerSettings.getMinReadCandidateFraction());
		assem.setMaxAverageDepth(assemblerSettings.getMaxAverageDepth());
		assem.setAverageDepthCeiling(assemblerSettings.getAverageDepthCeiling());
		assem.setDebug(assemblerSettings.isDebug());

		return assem;
	}

	private void init() throws IOException {
		
		File workingDir = new File(tempDir);
		if (workingDir.exists()) {
			if (!workingDir.delete()) {
				throw new IllegalStateException("Unable to delete: " + tempDir);
			}
		}

		if (!workingDir.mkdir()) {
			throw new IllegalStateException("Unable to create: " + tempDir);
		}
		
		File unalignedTempDir = new File(tempDir + "/unaligned");
		
		if (!unalignedTempDir.mkdir()) {
			throw new IllegalStateException("Unable to create: " + tempDir + "/unaligned");
		}
		
		new NativeLibraryLoader().load(tempDir);
		
		threadManager = new ThreadManager(numThreads);
	}
	
	public void setReference(String reference) {
		this.reference = reference;
	}
	
	public void setBwaIndex(String bwaIndex) {
		this.bwaIndex = bwaIndex;
	}

	public void setTempDir(String temp) {
		this.tempDir = temp;
	}

	public void setAssemblerSettings(AssemblerSettings settings) {
		this.assemblerSettings = settings;
	}
	
	public void setNumThreads(int numThreads) {
		this.numThreads = numThreads;
	}
		
	public void setShouldReprocessUnaligned(boolean shouldReprocessUnaligned) {
		this.shouldReprocessUnaligned = shouldReprocessUnaligned;
	}
	
	public void setMaxUnalignedReads(int maxUnalignedReads) {
		this.maxUnalignedReads = maxUnalignedReads;
	}
	
	public CompareToReference2 getC2r() {
		return this.c2r;
	}
	
	public int getMinMappingQuality() {
		return this.minMappingQuality;
	}
	
	public int getMaxInsertLength() {
		return this.maxInsertLength;
	}
	
	public int getMinInsertLength() {
		return this.minInsertLength;
	}
	
	public void setMaxInsertLength(int maxInsertLen) {
		this.maxInsertLength = maxInsertLen;
	}
	
	public void setMinInsertLength(int minInsertLen) {
		this.minInsertLength = minInsertLen;
	}
	
	boolean isFiltered(SAMRecord read) {
		return SAMRecordUtils.isFiltered(isPairedEnd, read);
	}

	public static void run(String[] args) throws Exception {
		
		System.err.println("Starting 0.97 ...");
		
		ReAlignerOptions options = new ReAlignerOptions();
		options.parseOptions(args);

		if (options.isValid()) {

			AssemblerSettings assemblerSettings = new AssemblerSettings();

			assemblerSettings.setKmerSize(options.getKmerSizes());
			assemblerSettings.setMinContigLength(options.getMinContigLength());
			assemblerSettings.setMinNodeFrequncy(options.getMinNodeFrequency());
			assemblerSettings.setMaxPotentialContigs(options
					.getMaxPotentialContigs());
			assemblerSettings.setMinUnalignedNodeFrequency(options.getMinUnalignedNodeFrequency());
			assemblerSettings.setMinBaseQuality(options.getMinBaseQuality());
			assemblerSettings.setMinReadCandidateFraction(options.getMinReadCandidateFraction());
			assemblerSettings.setMaxAverageDepth(options.getMaxAverageRegionDepth());
			assemblerSettings.setAverageDepthCeiling(options.getAverageDepthCeiling());
			assemblerSettings.setMinEdgeRatio(options.getMinEdgeRatio());
			assemblerSettings.setDebug(options.isDebug());
			assemblerSettings.setMaxNodes(options.getMaxNodes());

			ReAligner realigner = new ReAligner();
			realigner.setReference(options.getReference());
			realigner.setBwaIndex(options.getBwaIndex());
			realigner.setRegionsBed(options.getTargetRegionFile());
			realigner.setTempDir(options.getWorkingDir());
			realigner.setAssemblerSettings(assemblerSettings);
			realigner.setNumThreads(options.getNumThreads());
			realigner.isPairedEnd = options.isPairedEnd();
			realigner.minMappingQuality = options.getMinimumMappingQuality();
			realigner.hasPresetKmers = options.hasPresetKmers();
			realigner.isDebug = options.isDebug();

			long s = System.currentTimeMillis();
			
			realigner.reAlign(options.getInputFiles(), options.getOutputFiles());

			long e = System.currentTimeMillis();

			System.err.println("Elapsed seconds: " + (e - s) / 1000);
		} else {
			System.exit(-1);
		}
	}
	
	public static void main(String[] args) throws Exception {
//		String inp = "--in /home/lmose/dev/ayc/opt/mem/test_tumor.bam --kmer 43 --mc-mapq 25 --mcl 101 --mcr -1.0 --mnf 2 --umnf 2 --mpc 50000 --out /home/lmose/dev/ayc/opt/mem/test_tumor.abra.bam --ref /home/lmose/reference/test/test.fa --targets /home/lmose/dev/ayc/opt/mem/test.gtf --threads 2 --working /home/lmose/dev/ayc/opt/mem/work1 --mur 50000000 --no-unalign --mbq 20 --rcf .02";
		String inp = "--in /home/lmose/dev/ayc/opt/mem/test_tumor.bam --kmer 43 --out /home/lmose/dev/ayc/opt/mem/test_tumor.abra3.bam --ref /home/lmose/reference/test/test.fa --targets /home/lmose/dev/ayc/opt/mem/test2.bed --threads 2 --working /home/lmose/dev/ayc/opt/mem/work3";
		run(inp.split("\\s+"));
	}
}


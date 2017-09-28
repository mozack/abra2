/* Copyright 2013 University of North Carolina at Chapel Hill.  All rights reserved. */
package abra;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.net.InetAddress;
import java.net.URL;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.attribute.PosixFilePermission;
import java.nio.file.attribute.PosixFilePermissions;
import java.security.CodeSource;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Properties;
import java.util.Set;
import java.util.TreeSet;
import java.util.UUID;

import abra.JunctionUtils.JunctionComparator;
import abra.JunctionUtils.TooManyJunctionPermutationsException;
import abra.ReadEvaluator.Alignment;
import abra.ContigAligner.ContigAlignerResult;
import abra.SimpleMapper.Orientation;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMProgramRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.TextCigarCodec;

/**
 * ABRA's main entry point
 * 
 * @author Lisle E. Mose (lmose at unc dot edu)
 */
public class ReAligner {
	
	public static int MAX_REGION_LENGTH = 400;
	private static int MIN_REGION_REMAINDER = 200;
	public static int REGION_OVERLAP = 200;
	
	// Minimum sequence length recommended for use with bwa mem
	private static final int MIN_CONTIG_LENGTH = 70;
	
	// Cannot be larger than buffer in assembler.c
	private static final int MAX_KMER_SIZE = 199;
	
	// These must match constants in C code (Aligner matrix dimensions)
	private static final int MAX_CONTIG_LEN = 2000-1;
	private static final int MAX_REF_REGION_LEN = 5000-1;
	
	private SAMFileHeader[] samHeaders;
	
	private List<Feature> regions;

	private String regionsBed;

	private String reference;

	private AssemblerSettings assemblerSettings;
	
	private int numThreads;
	
	private String[] inputSams;
	
	private int readLength = -1;
	private int maxMapq = -1;
	private int minInsertLength = Integer.MAX_VALUE;
	private int maxInsertLength = -1;
	
	private boolean isPairedEnd = false;
	
	private BufferedWriter contigWriter = null;
	
	public static CompareToReference2 c2r;
	
	private ThreadManager threadManager;
	
	private int minMappingQuality;
	private double maxMismatchRate;
	
	private boolean isDebug;
	private boolean isSkipAssembly;
	private boolean isSkipUnmappedTrigger;
	private boolean useSoftClippedReads;
	private boolean useObservedIndels;
	private boolean useConsensusSeq;
	private boolean isKeepTmp;
	private boolean shouldSort;
	private String tmpDir;
	private int finalCompressionLevel;
	private int maxRealignDist;
	private int maxAssembledContigs;
	
	// If true, the input target file specifies kmer values
	private boolean hasPresetKmers = false;
	
	private String contigFile = null;
	
	// RNA specific
	private String junctionFile;
	private String gtfJunctionFile;
	private Set<Feature> junctions = new HashSet<Feature>();
	
	private ReverseComplementor rc = new ReverseComplementor();
	
	private String version = "unknown";
	private String cl = "unknown";
	
	private int[] swScoring;
	private int[] softClipParams;
	
	private int maxCachedReads = 0;
	private int maxReadsInRegion;
	
	private int minAnchorLen;
	private int maxAnchorMismatches;
	
	private SortedSAMWriter writer;
	
	private ChromosomeChunker chromosomeChunker;
	
	private String chromosomesToSkipRegex;
	private ChromosomeRegex chromosomeSkipRegex;
	
	private boolean shouldUnsetDuplicates;
	private String inputVcf;
	
	private Map<String, List<Variant>> knownVariants;;
	
	private boolean shouldCreateIndex;
	private boolean shouldUseGkl;
	
	private int ambiguousMapq;
	private double maxReadNoise;
	
	private int maxReadsInRamForSort;
	
	public void reAlign(String[] inputFiles, String[] outputFiles) throws Exception {
		
		this.inputSams = inputFiles;
		
		logStartupInfo(outputFiles);
				
		String tempDir = init();
		
		c2r = new CompareToReference2();
		c2r.init(this.reference);
		
		chromosomeChunker = new ChromosomeChunker(c2r);
		chromosomeChunker.init();

		Logger.info("Reading Input SAM Header and identifying read length");
		getSamHeaderAndReadLength();
		
		Logger.info("Read length: " + readLength);
		
		Logger.info("Loading target regions");
		loadRegions();
		loadJunctions();
		
		Clock clock = new Clock("Realignment");
		clock.start();
		
		if (contigFile != null) {
			contigWriter = new BufferedWriter(new FileWriter(contigFile, false));
		}
				
		for (int i=0; i<inputSams.length; i++) {
			
			SAMProgramRecord pg = new SAMProgramRecord("ABRA2");
			pg.setProgramVersion(this.version);
			pg.setCommandLine(cl);
			samHeaders[i].addProgramRecord(pg);			
		}
		
		writer = new SortedSAMWriter(outputFiles, tempDir.toString(), samHeaders, isKeepTmp, chromosomeChunker,
				finalCompressionLevel, shouldSort, maxRealignDist, shouldUnsetDuplicates, shouldCreateIndex, shouldUseGkl, maxReadsInRamForSort);

		// Spawn thread for each chromosome
		// TODO: Validate identical sequence dictionary for each input file
		
		for (int i=0; i<this.chromosomeChunker.getChunks().size(); i++) {
			spawnChromosomeThread(i);
		}
		
		Logger.info("Waiting for processing threads to complete");
		threadManager.waitForAllThreadsToComplete();
		
		if (contigWriter != null) {
			contigWriter.close();
		}
		
		clock.stopAndPrint();
		
		clock = new Clock("Sort and cleanup");
		clock.start();
		
		// Cut num threads in half to allow for async writer thread
		threadManager = new ThreadManager(Math.max(numThreads / 2, 1));
		
		for (int i=0; i<outputFiles.length; i++) {
			SortedSAMWriterRunnable thread = new SortedSAMWriterRunnable(threadManager, writer, i, inputSams[i]);
			threadManager.spawnThread(thread);
		}
		
		Logger.info("Waiting for writer threads to complete");
		threadManager.waitForAllThreadsToComplete();
		
		clock.stopAndPrint();
		
		Logger.info("Done.");
	}
	
	void processChromosomeChunk(int chromosomeChunkIdx) throws Exception {
		
		Feature chromosomeChunk = chromosomeChunker.getChunks().get(chromosomeChunkIdx);
		String chromosome = chromosomeChunk.getSeqname();
		
		Logger.info("Processing chromosome chunk: " + chromosomeChunk);
		Clock clock = new Clock("Chromosome: " + chromosomeChunk);
		clock.start();
		
		writer.initChromosomeChunk(chromosomeChunkIdx);
		
		MultiSamReader reader = new MultiSamReader(this.inputSams, this.minMappingQuality, this.isPairedEnd, chromosomeChunk);
		
		List<List<SAMRecordWrapper>> currReads = new ArrayList<List<SAMRecordWrapper>>();
		for (int i=0; i<this.inputSams.length; i++) {
			currReads.add(new ArrayList<SAMRecordWrapper>());
		}
		
		List<List<SAMRecordWrapper>> outOfRegionReads = new ArrayList<List<SAMRecordWrapper>>();
		for (int i=0; i<this.inputSams.length; i++) {
			outOfRegionReads.add(new ArrayList<SAMRecordWrapper>());
		}
		
		Map<Feature, Map<SimpleMapper, ContigAlignerResult>> regionContigs = new HashMap<Feature, Map<SimpleMapper, ContigAlignerResult>>();
		int readCount = 0;
		
		// Identify regions overlapping the current chromosome chunk
		List<Feature> chromosomeRegions = new ArrayList<Feature>();
		for (Feature region : regions) {			
			if (region.getSeqname().equals(chromosome)) {
				if (region.getStart() > chromosomeChunk.getStart()-MAX_REGION_LENGTH && region.getEnd() < chromosomeChunk.getEnd()+MAX_REGION_LENGTH) {
					chromosomeRegions.add(region);
				}
			}
		}
		
		List<Feature> chromosomeJunctions = new ArrayList<Feature>();
		for (Feature junction : junctions) {
			if (junction.getSeqname().equals(chromosome)) {
				chromosomeJunctions.add(junction);
			}
		}
		
		List<Variant> knownVariants = getKnownVariants(chromosome);
		
		Map<Feature, List<Feature>> regionJunctions = JunctionUtils.getRegionJunctions(chromosomeRegions, chromosomeJunctions, readLength, MAX_REGION_LENGTH);
		
		Map<Feature, List<Variant>> regionVariants = Variant.groupByRegion(chromosomeRegions, knownVariants);
		
		Set<Integer> regionsToProcess = new TreeSet<Integer>();
	
		int searchStartRegionIdx = 0;
		
		for (SAMRecordWrapper record : reader) {
			
			// If this is an unmapped read anchored by its mate, check rc flag
			SAMRecord read1 = record.getSamRecord();
			if (read1.getReadUnmappedFlag() && !read1.getMateUnmappedFlag()) {
				
				if (!read1.getReadNegativeStrandFlag() && !read1.getMateNegativeStrandFlag()) {
					// Both ends in forward orientation.  Reverse the unmapped read
					read1.setReadString(rc.reverseComplement(read1.getReadString()));
					read1.setBaseQualityString(rc.reverse(read1.getBaseQualityString()));
					read1.setReadNegativeStrandFlag(true);
					record.setUnalignedRc(true);
				} else if (read1.getReadNegativeStrandFlag() && read1.getMateNegativeStrandFlag()) {
					// Both ends in reverse orientation.  Reverse the unmapped read
					read1.setReadString(rc.reverseComplement(read1.getReadString()));
					read1.setBaseQualityString(rc.reverse(read1.getBaseQualityString()));
					read1.setReadNegativeStrandFlag(false);
					record.setUnalignedRc(true);
				}				
			}
			
			List<Integer> overlappingRegions = new ArrayList<Integer>();
			
			if (chromosomeRegions.size() > 0) {
				// Advance currRegion to current locus @ record start - 1000000
				// TODO: Use move dist instead of 1000000 ?
				Feature searchStartRegion = chromosomeRegions.get(searchStartRegionIdx);
				while (searchStartRegion != null && searchStartRegionIdx < chromosomeRegions.size()-1 && searchStartRegion.getEnd() < record.getAdjustedAlignmentStart() - 1000000) {
					searchStartRegionIdx += 1;
					searchStartRegion = chromosomeRegions.get(searchStartRegionIdx);
				}
				
				overlappingRegions = Feature.findAllOverlappingRegions(reader.getSAMFileHeader(), record, chromosomeRegions, searchStartRegionIdx);
			}
			
//			int regionIdx = Feature.findFirstOverlappingRegion(reader.getSAMFileHeader(), record, chromosomeRegions, currRegionIdx);
						
			// Identify next region that is a candidate for processing
			// Note: Splicing can cause reads to go in and out of a region
//			if (regionIdx >= 0) {
			if (!overlappingRegions.isEmpty()) {
				regionsToProcess.addAll(overlappingRegions);
			}
			
			// Cache read for processing at end of region
			currReads.get(record.getSampleIdx()).add(record);
			
			Iterator<Integer> regionIter = regionsToProcess.iterator();
			while (regionIter.hasNext()) {
				int regionToProcessIdx = regionIter.next();
				// If start position for current read is beyond current region, trigger assembly
				Feature currRegion = chromosomeRegions.get(regionToProcessIdx);
				if (record.getAdjustedAlignmentStart() > currRegion.getEnd() + this.readLength*2) {
					Logger.debug("Processing region: %s", currRegion);
					Map<SimpleMapper, ContigAlignerResult> mappedContigs = processRegion(currRegion, currReads, regionJunctions.get(currRegion), regionVariants.get(currRegion));
					Logger.debug("Region: %s assembled: %d contigs", currRegion, mappedContigs.keySet().size());
					regionContigs.put(currRegion, mappedContigs);
					// Remove curr region from list of regions to process
					regionIter.remove();
				}
			}
			
			/*
			
			// TODO: Consider dropping this...  Reads are out of scope when we've moved beyond them via standard processing?
			if (overlappingRegions.isEmpty()) {
				
				// Process out of region read and output if ready.
				List<SAMRecordWrapper> outOfRegionReadsForSample = outOfRegionReads.get(record.getSampleIdx());
				outOfRegionReadsForSample.add(record);
				
				if (outOfRegionReads.get(record.getSampleIdx()).size() > 2500) {

					for (SAMRecordWrapper outOfRegionRead : outOfRegionReadsForSample) {
						this.writer.addAlignment(record.getSampleIdx(), outOfRegionRead.getSamRecord());
					}
					
					outOfRegionReadsForSample.clear();
				}
			}
			*/
			
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
						
						// record == most recent read.  read = cached read
						if (record.getSamRecord().getAlignmentStart() - read.getSamRecord().getAlignmentStart() > MAX_READ_RANGE) {

							// Only output reads with start pos within current chromosomeChunk
							if (read.getSamRecord().getAlignmentStart() >= chromosomeChunk.getStart() &&
								read.getSamRecord().getAlignmentStart() <= chromosomeChunk.getEnd()) {
							
								sampleReadsToRemap.add(read);
							}
							
							iter.remove();
						}
					}					
				}

				// Remap out of scope reads
				long start = System.currentTimeMillis();
				int totalReads = remapReads(regionContigs, readsToRemap, chromosomeChunkIdx);
				long stop = System.currentTimeMillis();
//				Logger.debug("REMAP_READS_MSECS:\t%d\t%d\t%s:%d", (stop-start), totalReads, record.getSamRecord().getReferenceName(), record.getSamRecord().getAlignmentStart());
				
				// Remove out of scope region assemblies
				List<Feature> regionsToRemove = new ArrayList<Feature>();
				for (Feature region : regionContigs.keySet()) {
					if (getFirstStartPos(currReads)-region.getStart() > MAX_READ_RANGE) {
						regionsToRemove.add(region);
					}
				}
				
				for (Feature region : regionsToRemove) {
					Logger.debug("Removing contigs for region: %s", region.toString());
					regionContigs.remove(region);
				}

				String logPrefix = record.getSamRecord().getReferenceName() + ":" + record.getSamRecord().getAlignmentStart() + " : ";
				
				if (regionContigs.size() > 10) {
					Logger.debug("%s\tregionContigs size: %d", logPrefix, regionContigs.size());
				}
				
				
				//TODO: Revisit this.  Is it still necessary?
				int currReadsCount = 0;
				int idx = 0;
				boolean shouldClear = false;
				for (List<SAMRecordWrapper> reads : currReads) {
					currReadsCount += reads.size();
					
					if (reads.size() >= this.maxCachedReads) {
						shouldClear = true;
						Logger.warn(logPrefix + " Too many reads for sample: " + idx + " num_reads: " + reads.size() + ", clearing.");
					}
					
					idx += 1;
				}
				
				if (shouldClear) {
					for (int i=0; i<currReads.size(); i++) {
						List<SAMRecordWrapper> reads = currReads.get(i);

						for (SAMRecordWrapper read : reads) {
							this.writer.addAlignment(i, read, chromosomeChunkIdx);
						}
						
						reads.clear();
					}
				}

				if (currReadsCount > 250000) {
					Logger.info(logPrefix + "\tCurr reads size: " + currReadsCount);
				}
				
				int outOfRegionCount = 0;
				for (List<SAMRecordWrapper> reads : outOfRegionReads) {
					outOfRegionCount += reads.size();
				}

				if (outOfRegionCount > 10000) {
					Logger.info(logPrefix + "\tOut of region reads size: " + outOfRegionCount);
				}
			}
			
			readCount += 1;
		}
		
		// Attempt to process last region if applicable
		Iterator<Integer> regionIter = regionsToProcess.iterator();
		while (regionIter.hasNext()) {
			int regionToProcessIdx = regionIter.next();
			
			// We've moved beyond the current region
			// Assemble reads
			Feature region = chromosomeRegions.get(regionToProcessIdx);
			Logger.debug("Processing region: %s", region);
			Map<SimpleMapper, ContigAlignerResult> mappedContigs = processRegion(region, currReads, regionJunctions.get(region), regionVariants.get(region));
			Logger.debug("Region: %s assembled: %d contigs", region, mappedContigs.keySet().size());
			regionContigs.put(region, mappedContigs);
		}
		
		// Remap remaining reads
		remapReads(regionContigs, currReads, chromosomeChunkIdx);
		currReads.clear();
		regionContigs.clear();
		
		// Output remaining out of region reads
		for (int i=0; i<outOfRegionReads.size(); i++) {
			List<SAMRecordWrapper> outOfRegionReadsForSample = outOfRegionReads.get(i);
			for (SAMRecordWrapper outOfRegionRead : outOfRegionReadsForSample) {
				this.writer.addAlignment(i, outOfRegionRead, chromosomeChunkIdx);
			}
			
			outOfRegionReadsForSample.clear();
		}
		
		reader.close();
		
		writer.finishChromosomeChunk(chromosomeChunkIdx);
		
		clock.stopAndPrint();
	}
	
	private int getFirstStartPos(List<List<SAMRecordWrapper>> readsList) {
		int minPos = Integer.MAX_VALUE;
		for (List<SAMRecordWrapper> reads : readsList) {
			if (reads.size() > 0 && reads.get(0).getSamRecord().getAlignmentStart() < minPos) {
				minPos = reads.get(0).getSamRecord().getAlignmentStart(); 
			}
		}
		
		return minPos;
	}
	
	private void logStartupInfo(String[] outputFiles) throws IOException {
		
		Logger.info("ABRA version: " + this.version);
		
		int ctr = 0;
		for (String input : inputSams) {
			Logger.info("input" + ctr + ": " + input);
		}

		ctr = 0;
		for (String output : outputFiles) {
			Logger.info("output" + ctr + ": " + output);
		}
		
		Logger.info("regions: " + regionsBed);
		Logger.info("reference: " + reference);
		Logger.info("num threads: " + numThreads);
		Logger.info(assemblerSettings.getDescription());
		Logger.info("paired end: " + isPairedEnd);
		Logger.info("isSkipAssembly: " + isSkipAssembly);
		Logger.info("useSoftClippedReads: " + useSoftClippedReads);
		Logger.info("SW scoring: " + Arrays.toString(swScoring));
		Logger.info("Soft clip params: " + Arrays.toString(softClipParams));
		
		String javaVersion = System.getProperty("java.version");
		Logger.info("Java version: " + javaVersion);
		if (javaVersion.startsWith("1.6") || javaVersion.startsWith("1.5") || javaVersion.startsWith("1.4")) {
			throw new RuntimeException("Please upgrade to Java 7 or later to run ABRA.");
		}
		
		try {
			InetAddress localhost = java.net.InetAddress.getLocalHost();
			String hostname = localhost.getHostName();
			Logger.info("hostname: " + hostname);
		} catch (Throwable t) {
			Logger.error("Error getting hostname: " + t.getMessage());
		}
	}
		
	private void spawnChromosomeThread(int chromosomeChunkIdx) throws InterruptedException {
		ReAlignerRunnable thread = new ReAlignerRunnable(threadManager, this, chromosomeChunkIdx);
		Logger.debug("Queuing thread for chromosome: " + chromosomeChunkIdx);
		threadManager.spawnThread(thread);
	}
	
	private synchronized void appendContigs(String contigs) throws IOException {
		if (contigWriter != null) {
			contigWriter.write(contigs);
		}
	}
	
	private void remapRead(ReadEvaluator readEvaluator, SAMRecord read, int origEditDist) {
		
		Alignment alignment = readEvaluator.getImprovedAlignment(origEditDist, read, c2r);
		if (alignment != null) {
			
			if (alignment == Alignment.AMBIGUOUS) {
				// Read maps equally well to reference and multiple differing contigs.  Flag with mapq of 1.
				if (ambiguousMapq >= 0 && read.getMappingQuality() > ambiguousMapq) {
					read.setMappingQuality(ambiguousMapq);
				}
			} else if (origEditDist == alignment.numMismatches && (read.getAlignmentStart() != alignment.pos || !read.getCigarString().equals(alignment.cigar))) {
				// Read maps ambiguously.  Downgrade mapping quality
				if (ambiguousMapq >= 0 && read.getMappingQuality() > ambiguousMapq) {
					read.setMappingQuality(ambiguousMapq);
				}
			} else if (Math.abs(read.getAlignmentStart() - alignment.pos) > maxRealignDist) {
				Logger.trace("Not moving read: " + read.getReadName() + " from: " + read.getAlignmentStart() + " to: " + alignment.pos);
			} else if (origEditDist > alignment.numMismatches) {
				
				SAMRecord orig = read.deepCopy();
				
				int readPos = alignment.pos;
				
				// Set contig alignment info for all reads that map to contigs (even if read is unchanged)
				String ya = alignment.chromosome + ":" + alignment.contigPos + ":" + alignment.contigCigar;
				
				// If no change to alignment, just record the YA tag
				if (!read.getReadUnmappedFlag() && read.getAlignmentStart() == readPos && read.getCigarString().equals(alignment.cigar)) {
					read.setAttribute("YA", ya);
				}
				
				// If the read has actually moved to an improved alignment, update
				if (read.getReadUnmappedFlag() || read.getAlignmentStart() != readPos || !read.getCigarString().equals(alignment.cigar)) {
					
					read.setAttribute("YA", ya);
	
					// Original alignment info
					String yo;
					if (!read.getReadUnmappedFlag()) {
						String origOrientation = read.getReadNegativeStrandFlag() ? "-" : "+";
						yo = read.getReferenceName() + ":" + read.getAlignmentStart() + ":" + origOrientation + ":" + read.getCigarString();
					} else {
						read.setReadUnmappedFlag(false);
						read.setMappingQuality(this.maxMapq);
						// Original alignment start info used in sort phase
						yo = "N/A:" + read.getAlignmentStart();
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
					
					String distTag = "NM";
					if (orig.getAttribute("NM") == null && orig.getAttribute("nM") != null) {
						distTag = "nM";
					}
					
					// Updated edit distance
					read.setAttribute(distTag, SAMRecordUtils.getEditDistance(read, c2r, false));
					
					//TODO: Compute mapq intelligently???
					read.setMappingQuality(Math.min(read.getMappingQuality()+10, this.maxMapq));
					
					// Check for realignments that are too noisy
					int numMismatches = c2r.numHighQualityMismatches(read, 20, false); 
					int numIndels = SAMRecordUtils.getNumIndels(read);
					
					int noise = numMismatches + (numIndels * 2);
					
					int maxNoise = (int) (SAMRecordUtils.getUnclippedLength(read) * maxReadNoise);
					
					if (noise > maxNoise || 
						read.getAlignmentEnd() >= c2r.getChromosomeLength(read.getReferenceName())-1) {
						// Read is too noisy or read maps off end of chromosome, revert
						read.setAlignmentStart(orig.getAlignmentStart());
						read.setCigar(orig.getCigar());
						read.setReadNegativeStrandFlag(orig.getReadNegativeStrandFlag());
						read.setAttribute("YA", null);
						read.setAttribute("YO", null);
						read.setAttribute("YM", null);
						read.setAttribute("YX", null);
						read.setAttribute(distTag, orig.getAttribute(distTag));
						read.setMappingQuality(orig.getMappingQuality());
						read.setReadUnmappedFlag(orig.getReadUnmappedFlag());
					}
				}
			}
		}
	}
	
	private int remapReads(Map<Feature, Map<SimpleMapper, ContigAlignerResult>> mappedContigs,
			List<List<SAMRecordWrapper>> readsList, int chromosomeChunkIdx) throws Exception {
		
		ReadEvaluator readEvaluator = new ReadEvaluator(mappedContigs);
		
		int sampleIdx = 0;
		int totalReads = 0;
		// For each sample.
		for (List<SAMRecordWrapper> reads : readsList) {
			
			// For each read.
			for (SAMRecordWrapper readWrapper : reads) {
				totalReads += 1;
				SAMRecord read = readWrapper.getSamRecord();
								
				if (read.getMappingQuality() >= this.minMappingQuality || read.getReadUnmappedFlag()) {
					
					// Don't remap reads with distant mate
					// Always allow single end to pass this check
					if (!read.getReadPairedFlag() ||
						(Math.abs(read.getAlignmentStart() - read.getMateAlignmentStart()) < maxRealignDist &&
								read.getReferenceName().equals(read.getMateReferenceName()))) {
					
						// TODO: Use NM tag if available (need to handle soft clipping though!)
						int origEditDist = SAMRecordUtils.getEditDistance(read, c2r, true);
		//				int origEditDist = c2r.numMismatches(read);
											
						remapRead(readEvaluator, read, origEditDist);
					}
				}
			}

			// Output all reads for this sample
			for (SAMRecordWrapper read : reads) {
				this.writer.addAlignment(sampleIdx, read, chromosomeChunkIdx);
			}
			
			sampleIdx += 1;
		}
		
		return totalReads;
	}
	
	private List<List<SAMRecordWrapper>> subsetReads(Feature region, List<List<SAMRecordWrapper>> readsList) {
		List<List<SAMRecordWrapper>> subset = new ArrayList<List<SAMRecordWrapper>>();
		
		// Initialize per sample lists
		for (List<SAMRecordWrapper> origSample : readsList) {
			
			// Track read pair ends
			Map<String, SAMRecordWrapper> firstReads = new HashMap<String, SAMRecordWrapper>();
			Map<String, SAMRecordWrapper> secondReads = new HashMap<String, SAMRecordWrapper>();
			
			List<SAMRecordWrapper> subsetSample = new ArrayList<SAMRecordWrapper>();
			subset.add(subsetSample);
			
			for (SAMRecordWrapper read : origSample) {
				if (region.overlapsRead(read.getSamRecord())) {
					subsetSample.add(read);
					
					if (read.getSamRecord().getReadPairedFlag()) {
						if (read.getSamRecord().getFirstOfPairFlag() && SAMRecordUtils.isPrimary(read.getSamRecord())) {
							firstReads.put(read.getSamRecord().getReadName() + "_" + read.getSamRecord().getAlignmentStart(), read);
						} else if (read.getSamRecord().getSecondOfPairFlag() && SAMRecordUtils.isPrimary(read.getSamRecord())) {
							secondReads.put(read.getSamRecord().getReadName() + "_" + read.getSamRecord().getAlignmentStart(), read);
						}
					}
				}
			}
			
			for (SAMRecordWrapper read : subsetSample) {
				if (SAMRecordUtils.hasPossibleAdapterReadThrough(read.getSamRecord(), firstReads, secondReads)) {
					read.setShouldAssemble(false);
				}
				
				// If reads overlap, attempt to generate merged sequence
				if (!read.hasMergedSeq()) {
					SAMRecordUtils.mergeReadPair(read, firstReads, secondReads);
				}
			}
		}
		
		return subset;
	}
	
	private List<Feature> getExtraJunctions(ContigAlignerResult result, List<Feature> junctions, List<Feature> junctions2) {
		// Look for junctions with adjacent deletions.
		// Treat these as putative novel junctions.
		Set<Feature> junctionSet = new HashSet<Feature>(junctions);
		junctionSet.addAll(junctions2);
		List<Feature> extraJunctions = new ArrayList<Feature>();
		
		Cigar cigar = TextCigarCodec.decode(result.getCigar());
		
		boolean isInGap = false;
		int gapStart = -1;
		int gapLength = -1;
		int numElems = -1;
		int refOffset = 0;
		
		for (CigarElement elem : cigar.getCigarElements()) {
			if (elem.getOperator() == CigarOperator.D || elem.getOperator() == CigarOperator.N) {
				if (!isInGap) {
					isInGap = true;
					gapStart = refOffset;
					gapLength = elem.getLength();
					numElems = 1;
				} else {
					gapLength += elem.getLength();
					numElems += 1;
				}
			} else {
				if (isInGap) {
					
					if (numElems > 1) {
						long start = result.getGenomicPos() + gapStart;
						long end = start + gapLength;
						Feature junc = new Feature(result.getChromosome(), start, end-1);
						if (!junctionSet.contains(junc)) {
							Logger.info("Extra junction idenfified: %s", junc);
							extraJunctions.add(junc);
						}
					}
					
					isInGap = false;
					gapStart = -1;
					gapLength = -1;
					numElems = -1;
				}
			}
			
			if (elem.getOperator() == CigarOperator.M || 
				elem.getOperator() == CigarOperator.D || 
				elem.getOperator() == CigarOperator.N ||
				elem.getOperator() == CigarOperator.X ||
				elem.getOperator() == CigarOperator.EQ) {
				
				refOffset += elem.getLength();
			}
		}
		
		return extraJunctions;
	}
	
	private List<Feature> getExonSkippingJunctions(ContigAlignerResult result, List<Feature> junctions) {

		// Handles special case where an exon skipping junction causes large gaps with a tiny (~1)
		// number of bases mapped somewhere in the gap
		Set<Feature> junctionSet = new HashSet<Feature>(junctions);
		List<Feature> extraJunctions = new ArrayList<Feature>();
		
		Cigar cigar = TextCigarCodec.decode(result.getCigar());
		
		boolean isInGap = false;
		int gapStart = -1;
		int gapLength = -1;
		int numElems = -1;
		int refOffset = 0;
		
		int maxBasesInMiddle = 5;
		int middleBases = -1;
		
		CigarOperator prev = CigarOperator.M;
		
		for (CigarElement elem : cigar.getCigarElements()) {
			if (elem.getOperator() == CigarOperator.D || elem.getOperator() == CigarOperator.N) {
				if (!isInGap) {
					isInGap = true;
					gapStart = refOffset;
					gapLength = elem.getLength();
					numElems = 1;
					middleBases = 0;
				} else {
					gapLength += elem.getLength();
					numElems += 1;
				}
			} else {
				if (isInGap) {
					
					if (elem.getLength() + middleBases < maxBasesInMiddle) {
						middleBases += elem.getLength();
					} else {
						
						if (numElems > 1 && middleBases > 0 && (prev == CigarOperator.D || prev == CigarOperator.N)) {
							long start = result.getGenomicPos() + gapStart;
							long end = start + gapLength;
							
							// Find junction start / end points closest to gap start / end point
							long closestStart = junctions.get(0).getStart();
							long closestEnd = junctions.get(0).getEnd();
							for (Feature junction : junctions) {
								if (Math.abs(junction.getStart()-start) < Math.abs(closestStart-start)) {
									closestStart = junction.getStart();
								}

								if (Math.abs(junction.getEnd()-end) < Math.abs(closestEnd-end)) {
									closestEnd = junction.getEnd();
								}
							}
							
							if (closestEnd > closestStart+1) {
								Feature junc = new Feature(result.getChromosome(), closestStart, closestEnd);
								if (!junctionSet.contains(junc)) {
									Logger.info("Potential exon skipping junction idenfified: %s", junc);
									extraJunctions.add(junc);
								}
							}
						}
						
						isInGap = false;
						gapStart = -1;
						gapLength = -1;
						numElems = -1;
						middleBases = -1;
					}
				}
			}
			
			if (elem.getOperator() == CigarOperator.M || 
				elem.getOperator() == CigarOperator.D || 
				elem.getOperator() == CigarOperator.N ||
				elem.getOperator() == CigarOperator.X ||
				elem.getOperator() == CigarOperator.EQ) {
				
				refOffset += elem.getLength();
			}
			
			prev = elem.getOperator();
		}
		
		return extraJunctions;
	}

	
	private ContigAlignerResult alignContig(Feature region, String contig, ContigAligner ssw, List<ContigAligner> sswJunctions, List<Feature> allJunctions,
			int chromosomeLength) {
		
		ContigAlignerResult bestResult = null;
		
		if (contig.length() > MAX_CONTIG_LEN) {
			Logger.warn(String.format("In Region: %s, contig too long: [%s]", region, contig));
		} else {

			int bestScore = -1;
			
			ContigAlignerResult sswResult;
			for (ContigAligner sswJunc : sswJunctions) {
				sswResult = sswJunc.align(contig);
				if (sswResult != null && sswResult.getScore() > bestScore) {
					bestScore = sswResult.getScore();
					bestResult = sswResult;
				}
			}
			
			sswResult = ssw.align(contig);
			if (sswResult != null && sswResult.getScore() > bestScore) {
				bestScore = sswResult.getScore();
				bestResult = sswResult;
			}
	
			if (bestResult != null && bestResult != ContigAlignerResult.INDEL_NEAR_END) {
				
				if (!allJunctions.isEmpty()) {
					// Check for additional potential exon skipping junctions masked by a base or 2 interrupting the gap
					// Using annotated exons here
					List<Feature> extraJunctions = getExonSkippingJunctions(bestResult, allJunctions);
					if (!extraJunctions.isEmpty()) {
						
						List<Feature> combinedJunctions = new ArrayList<Feature>(allJunctions);
						combinedJunctions.addAll(extraJunctions);
						
						Collections.sort(combinedJunctions, new JunctionComparator());
						
						List<List<Feature>> junctionPermutations = new ArrayList<List<Feature>>();
						try {
							junctionPermutations = JunctionUtils.combineJunctions(region, combinedJunctions, new HashSet<Feature>(extraJunctions), MAX_REGION_LENGTH, this.readLength);
						} catch (TooManyJunctionPermutationsException e) {
							Logger.warn("TOO_MANY_POTENTIAL_JUNCTION_PERMUTATIONS: " + region.getDescriptor());
						}
						
						for (List<Feature> permutation : junctionPermutations) {
							boolean hasExtra = false;
							for (Feature junc : permutation) {
								if (extraJunctions.contains(junc)) {
									hasExtra = true;
									break;
								}
							}
							
							if (hasExtra) {
								ContigAligner aligner = getContigAlignerForJunctionPermutation(permutation, region, chromosomeLength);
								if (aligner != null) {
									sswResult = aligner.align(contig);
									if (sswResult != null && sswResult.getScore() > bestScore) {
										bestScore = sswResult.getScore();
										bestResult = sswResult;
									}
								}
							}
						}
					}
					
					// Check for deletion adjacent to intron (i.e. skipped exon or unannotated splice)
					// Not relying on annotated exons here.
					extraJunctions = getExtraJunctions(bestResult, allJunctions, extraJunctions);
					if (!extraJunctions.isEmpty()) {
						
						List<Feature> combinedJunctions = new ArrayList<Feature>(allJunctions);
						combinedJunctions.addAll(extraJunctions);
						
						Collections.sort(combinedJunctions, new JunctionComparator());
						
						List<List<Feature>> junctionPermutations = new ArrayList<List<Feature>>();
						try {
							junctionPermutations = JunctionUtils.combineJunctions(region, combinedJunctions, new HashSet<Feature>(extraJunctions), MAX_REGION_LENGTH, this.readLength);
						} catch (TooManyJunctionPermutationsException e) {
							Logger.warn("TOO_MANY_POTENTIAL_JUNCTION_PERMUTATIONS: " + region.getDescriptor());
						}
						
						for (List<Feature> permutation : junctionPermutations) {
							boolean hasExtra = false;
							for (Feature junc : permutation) {
								if (extraJunctions.contains(junc)) {
									hasExtra = true;
									break;
								}
							}
							
							if (hasExtra) {
								ContigAligner aligner = getContigAlignerForJunctionPermutation(permutation, region, chromosomeLength);
								if (aligner != null) {
									sswResult = aligner.align(contig);
									if (sswResult != null && sswResult.getScore() > bestScore) {
										bestScore = sswResult.getScore();
										bestResult = sswResult;
									}
								}
							}
						}
					}
				}
				
				Logger.debug("BEST_SSW: %d : %s : %d: %d : %s",
						bestResult.getGenomicPos(), bestResult.getCigar(), bestResult.getRefPos(), bestResult.getScore(), bestResult.getSequence());
				
			} else {
				Logger.debug("NO_SSW: %s", contig);
			}
			
			//TODO: Check for tie scores with different final alignment
		}
		
		return bestResult;
	}
	
	private boolean assemble(List<ContigAlignerResult> results, Feature region, 
			String refSeq, List<String> bams, List<List<SAMRecordWrapper>> readsList, ContigAligner contigAligner,
			List<ContigAligner> junctionAligners, int mnf, int mbq, double mer, List<Feature> junctions,
			int chromosomeLength, int maxNumContigs) throws IOException {
		
		boolean shouldRetry = false;
		
		NativeAssembler assem = (NativeAssembler) newAssembler(region);
		List<Feature> regions = new ArrayList<Feature>();
		regions.add(region); 
		StringBuffer readBuffer = new StringBuffer();
		String contigs = assem.assembleContigs(bams, regions, region.getDescriptor(), true, this, c2r, readsList, mnf, mbq, mer, readBuffer);
		
		if (!contigs.equals("<ERROR>") && !contigs.equals("<REPEAT>") && !contigs.isEmpty()) {
			
			if (contigWriter != null) {
				appendContigs(contigs);
			}
			
			List<ScoredContig> scoredContigs = ScoredContig.convertAndFilter(contigs, maxNumContigs, readBuffer);
			
			// Map contigs to reference
			for (ScoredContig contig : scoredContigs) {
				// Filter contigs that match the reference
				if (!refSeq.contains(contig.getContig())) {
					
					ContigAlignerResult sswResult = alignContig(region, contig.getContig(), contigAligner, junctionAligners, junctions, chromosomeLength);
					
					if (sswResult == ContigAlignerResult.INDEL_NEAR_END) {
						shouldRetry = true;
					} else if (sswResult != null) {
						// TODO: In multi-region processing, check to ensure identical contigs have identical mappings
						results.add(sswResult);
					}
				}
			}
		} 
		
		return shouldRetry;
	}
	
	private ContigAligner getContigAlignerForJunctionPermutation(List<Feature> junctionPerm, Feature region, int chromosomeLength) {
		
		ContigAligner contigAligner = null;
		
		// List of junction positions within localized reference
		List<Integer> junctionPos = new ArrayList<Integer>();
		// List of junction lengths within localized reference
		List<Integer> junctionLengths = new ArrayList<Integer>();
		
		StringBuffer juncSeq = new StringBuffer();
		
		int refStart = Math.max((int) junctionPerm.get(0).getStart() - (int) region.getLength() - this.readLength*2, 1);
		String leftSeq = c2r.getSequence(region.getSeqname(), refStart, (int) junctionPerm.get(0).getStart() - refStart);
		juncSeq.append(leftSeq);
		junctionPos.add(leftSeq.length());
		junctionLengths.add((int) junctionPerm.get(0).getLength()+1);
		
		boolean isJunctionGapTooBig = false;
		
		for (int i=1; i<junctionPerm.size(); i++) {
			int midStart = (int) junctionPerm.get(i-1).getEnd()+1;
			String middleSeq = c2r.getSequence(region.getSeqname(), midStart, (int) junctionPerm.get(i).getStart() - midStart);
			if (middleSeq.length() > region.getLength()*2) {
				isJunctionGapTooBig = true;
				break;
			}
			juncSeq.append(middleSeq);
			junctionPos.add(juncSeq.length());
			junctionLengths.add((int) junctionPerm.get(i).getLength()+1);
		}
		
		// TODO: Tighten this up...
		if (!isJunctionGapTooBig && juncSeq.length() < region.getLength()*10) {
			
			// Sequence on right of last junction
			// Junction stop is exclusive, so add 1 to starting position (junction end + 1)
			Feature lastJunction = junctionPerm.get(junctionPerm.size()-1);
			int rightStart = (int) lastJunction.getEnd()+1;
			int rightStop = Math.min((int) lastJunction.getEnd() + (int) region.getLength() + this.readLength*2, chromosomeLength-1);
			
			if (rightStop-rightStart > 0) {
				String rightSeq = c2r.getSequence(region.getSeqname(), rightStart, rightStop-rightStart);
				juncSeq.append(rightSeq);
				// Junction pos and length should already be added
				if (juncSeq.length() > MAX_REF_REGION_LEN) {
					// Make sure we don't blow up the hardcoded size C matrix
					Logger.warn("Junction Ref Seq to long: " + juncSeq.toString());
					
				} else {
					contigAligner = new ContigAligner(juncSeq.toString(), region.getSeqname(), refStart, this.readLength, minAnchorLen, maxAnchorMismatches, junctionPos, junctionLengths);
				}
			}
		}

		return contigAligner;
	}
	
	public Map<SimpleMapper, ContigAlignerResult> processRegion(Feature region, List<List<SAMRecordWrapper>> reads, List<Feature> junctions, List<Variant> knownVariants) throws Exception {
		
		long start = System.currentTimeMillis();
		if (isDebug) {
			Logger.info("Processing region: " + region.getDescriptor());
		}
		
		if (region.getLength() > 10000) {
			throw new IllegalArgumentException("Region too big: [" + region + "]");
		}
		
		Map<SimpleMapper, ContigAlignerResult> mappedContigs = new HashMap<SimpleMapper, ContigAlignerResult>();
		
		List<List<SAMRecordWrapper>> readsList = subsetReads(region, reads);
		
		boolean isRegionOk = true;
		for (List<SAMRecordWrapper> sampleReads : readsList) {
			
			//TODO: Don't allow these reads to remap to neighboring regions.
			if (maxReadsInRegion < 0 || sampleReads.size() > this.maxReadsInRegion) {
				Logger.info("Too many reads in %s: %d", region, sampleReads.size());
				isRegionOk = false;
				break;
			}
		}
		

		
		int assembledContigCount = 0;
		int nonAssembledContigCount = 0;
		int juncPermCount = 0;
		
		if (isRegionOk) {
			List<String> bams = new ArrayList<String>(Arrays.asList(this.inputSams));
			
			// Get reference sequence matching current region (pad by 2 read lengths on each side)
			int chromosomeLength = c2r.getReferenceLength(region.getSeqname());
			int refSeqStart = Math.max((int) region.getStart() - this.readLength*2, 1);
			int refSeqLength = Math.min((int) region.getLength() + this.readLength*4, chromosomeLength-1);
			
			String refSeq = c2r.getSequence(region.getSeqname(), refSeqStart, refSeqLength);
			
			ContigAligner ssw = new ContigAligner(refSeq, region.getSeqname(), refSeqStart, this.readLength, minAnchorLen, maxAnchorMismatches);
			
			List<ContigAligner> junctionAligners = new ArrayList<ContigAligner>();
			
//			List<List<Feature>> junctionPermutations = JunctionUtils.combineJunctions(junctions, this.readLength);
			
			List<List<Feature>> junctionPermutations = new ArrayList<List<Feature>>();
			try {
				junctionPermutations = JunctionUtils.combineJunctions(region, junctions, new HashSet<Feature>(), MAX_REGION_LENGTH, this.readLength);
			} catch (TooManyJunctionPermutationsException e) {
				Logger.warn("TOO_MANY_POTENTIAL_JUNCTION_PERMUTATIONS: " + region.getDescriptor());
			}
			
			Logger.debug("NUM_JUNCTION_PERMUTATIONS:\t%d\t%s", junctionPermutations.size(), region);
			
			if (junctionPermutations.size() > JunctionUtils.MAX_JUNCTION_PERMUTATIONS) {
				Logger.warn("TOO_MANY_JUNCTION_PERMUTATIONS: " + region.getDescriptor() + "\t" + junctionPermutations.size());
			} else {
			
				juncPermCount = junctionPermutations.size();
				
				for (List<Feature> junctionPerm : junctionPermutations) {
					ContigAligner aligner = getContigAlignerForJunctionPermutation(junctionPerm, region, chromosomeLength);
					if (aligner != null) {
						Logger.debug("JUNC_REF_SEQ:\t%s\t%d", region.getDescriptor(), aligner.ref.length());
						junctionAligners.add(aligner);
					}
				}
							
				// Assemble contigs
				if (this.isSkipAssembly || region.getKmer() > this.readLength-15) {
					Logger.debug("Skipping assembly of region: " + region.getDescriptor() + " - " + region.getKmer());
				} else {
					
					//
					// Restrict # of contigs if junction count grows high.
					// TODO: Paramaterize
					int maxCombos = 1024;
					int maxNumContigs = junctionPermutations.size() == 0 ? maxAssembledContigs : Math.min(maxAssembledContigs, maxCombos/junctionPermutations.size());
					
					if (maxNumContigs != maxAssembledContigs) {
						Logger.info("MAX_ASSEM_CONTIG\t%s\t%d", region, maxNumContigs);
					}
					
					List<ContigAlignerResult> results = new ArrayList<ContigAlignerResult>();
					boolean shouldRetry = assemble(results, region, refSeq, bams, readsList, ssw, junctionAligners,
							assemblerSettings.getMinNodeFrequncy(), assemblerSettings.getMinBaseQuality(),
							assemblerSettings.getMinEdgeRatio(), junctions, chromosomeLength, maxNumContigs);
					
					if (shouldRetry) {
						Logger.debug("RETRY_ASSEMBLY: %s", region);
						// Indel near edge of contig indicates that we may have a low coverage indel encountered.
						// Try to reassemble using less stringent pruning to see if we can get greater coverage.
						results.clear();
						assemble(results, region, refSeq, bams, readsList, ssw, junctionAligners,
								assemblerSettings.getMinNodeFrequncy()/2, assemblerSettings.getMinBaseQuality()/2,
								assemblerSettings.getMinEdgeRatio()/2.0, junctions, chromosomeLength, maxNumContigs);
					}
					
					for (ContigAlignerResult sswResult : results) {
						mappedContigs.put(new SimpleMapper(sswResult.getSequence(), maxMismatchRate), sswResult);
					}
					
					assembledContigCount = mappedContigs.size();
				}
				
				if (useSoftClippedReads || useObservedIndels || (knownVariants != null && knownVariants.size() > 0)) {
					Logger.debug("Processing non-assembled contigs for region: [" + region + "]");
					// Go through artificial contig generation using indels observed in the original reads
					AltContigGenerator altContigGenerator = new AltContigGenerator(softClipParams[0], softClipParams[1], softClipParams[2], softClipParams[3],
							useObservedIndels, useSoftClippedReads, useConsensusSeq, minMappingQuality);
					Collection<String> altContigs = altContigGenerator.getAltContigs(readsList, c2r, readLength, junctionPermutations.size(), region, knownVariants);
					
					nonAssembledContigCount = altContigs.size();
					
					for (String contig : altContigs) {
						// TODO: Check to see if this contig is already in the map before aligning
						
						ContigAlignerResult sswResult = alignContig(region, contig, ssw, junctionAligners, junctions, chromosomeLength);
						if (sswResult != null && sswResult != ContigAlignerResult.INDEL_NEAR_END) {
							// Set as secondary for remap prioritization
							sswResult.setSecondary(true);
							// Store for read mapping
							mappedContigs.put(new SimpleMapper(sswResult.getSequence(), maxMismatchRate), sswResult);
						}
					}
				}
			}
		}
		
		long stop = System.currentTimeMillis();
		
		synchronized(this.getClass()) {
			Logger.info("PROCESS_REGION_MSECS:\t%s\t%d\t%d\t%d\t%d", region.getDescriptor(), (stop-start), 
					assembledContigCount, nonAssembledContigCount, juncPermCount);
		}
		
		return mappedContigs;
	}
	
	// Pair up junctions that could be spanned by a single read
	protected List<Pair<Feature, Feature>> pairJunctions(List<Feature> junctions, int maxDist) {
		List<Pair<Feature, Feature>> junctionPairs = new ArrayList<Pair<Feature, Feature>>();
		
		for (Feature junc1 : junctions) {
			for (Feature junc2 : junctions) {
				if (junc1.getEnd() < junc2.getStart() && junc1.getEnd() + maxDist >= junc2.getStart()) {
					junctionPairs.add(new Pair<Feature, Feature>(junc1, junc2));
				}
			}
		}
		
		return junctionPairs;
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
	
	private List<Feature> getRegionsNoBed(int readLength, SAMFileHeader header) throws IOException {
		
		List<Feature> regions = new ArrayList<Feature>();
		
		List<SAMSequenceRecord> refSeq = header.getSequenceDictionary().getSequences();
		
		for (SAMSequenceRecord seq : refSeq) {
			if (!chromosomeSkipRegex.matches(seq.getSequenceName())) {
				Feature region = new Feature(seq.getSequenceName(), 1, seq.getSequenceLength());
				regions.add(region);
			}
		}
		
		regions = RegionLoader.collapseRegions(regions, readLength);
		regions = splitRegions(regions);
		
		return regions;
	}
		
	private void loadRegions() throws IOException {
		
		if (regionsBed != null) {
			Logger.info("Loading target regions from : " + regionsBed);
			this.regions = getRegions(regionsBed, readLength, hasPresetKmers);
		} else {
			Logger.info("No target bed file specified.  Gathering regions using SAM header");
			this.regions = getRegionsNoBed(readLength, this.samHeaders[0]);
		}
		
		Logger.info("Num regions: " + regions.size());
		
		if (Logger.LEVEL == Logger.Level.TRACE) {
			for (Feature region : regions) {
				Logger.trace("%s\t%d\t%d\t%d", region.getSeqname(), region.getStart(), region.getEnd(), region.getKmer());
			}
		}
	}
	
	private void loadJunctions() throws IOException {
		if (this.gtfJunctionFile != null) {
			this.junctions = JunctionUtils.loadJunctionsFromGtf(gtfJunctionFile);
		}
		
		if (this.junctionFile != null) {
			RegionLoader loader = new RegionLoader();
			List<Feature> observedJunctions = loader.load(junctionFile, false);
			Logger.info("Loaded " + observedJunctions.size() + " observed junctions");
			junctions.addAll(observedJunctions);
		}
		
		Logger.info("Total junctions input: " + junctions.size());
	}

	public void setRegionsBed(String bedFile) {
		this.regionsBed = bedFile;
	}

	private void getSamHeaderAndReadLength() throws IOException {
		
		Logger.info("Identifying header and determining read length");
		
		this.samHeaders = new SAMFileHeader[this.inputSams.length];
		
		for (int i=0; i<this.inputSams.length; i++) {
		
			SamReader reader = SAMRecordUtils.getSamReader(inputSams[i]);
			try {
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
					
					cnt += 1;
				}
				
				// Allow some fudge in insert length
				minInsertLength = Math.max(minInsertLength - 2*readLength, 0);
				maxInsertLength = maxInsertLength + 2*readLength;
				
			} finally {
				reader.close();
			}
		}
		
		Logger.info("Min insert length: " + minInsertLength);
		Logger.info("Max insert length: " + maxInsertLength);
				
		Logger.info("Max read length is: " + readLength);
		if (assemblerSettings.getMinContigLength() < 1) {
			assemblerSettings.setMinContigLength(Math.max(readLength+1, MIN_CONTIG_LENGTH));
		}
		Logger.info("Min contig length: " + assemblerSettings.getMinContigLength());
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
		
		int maxKmerSize = readLength-5; 
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
		assem.setSkipUnmappedTrigger(this.isSkipUnmappedTrigger);

		return assem;
	}
	
	private void deleteOnExit(File file) {
		if (!isKeepTmp) {
			file.deleteOnExit();
		}
	}

	private String init() throws IOException {
		
		if (tmpDir == null) {
			tmpDir = System.getProperty("java.io.tmpdir");
		} else {
			System.setProperty("java.io.tmpdir", tmpDir);
		}
		
		this.chromosomeSkipRegex = new ChromosomeRegex(chromosomesToSkipRegex);
		
		ContigAligner.init(swScoring);
		
		Set<PosixFilePermission> perms = new HashSet<PosixFilePermission>();
        perms.add(PosixFilePermission.OWNER_READ);
        perms.add(PosixFilePermission.OWNER_WRITE);
        perms.add(PosixFilePermission.OWNER_EXECUTE);
        perms.add(PosixFilePermission.GROUP_READ);
        perms.add(PosixFilePermission.GROUP_EXECUTE);

		Path tempDir = Files.createTempDirectory("abra2_" + UUID.randomUUID(), PosixFilePermissions.asFileAttribute(perms));
		deleteOnExit(tempDir.toFile());
		
		Logger.info("Using temp directory: " + tempDir.toString());
		
		new NativeLibraryLoader().load(tempDir.toString(), NativeLibraryLoader.ABRA, false);
//		new NativeLibraryLoader().load(tempDir.toString(), NativeLibraryLoader.SSW, false);
//		new NativeLibraryLoader().load(tempDir.toString(), NativeLibraryLoader.SSW_JNI, false);
		new NativeLibraryLoader().load(tempDir.toString(), NativeLibraryLoader.DEFLATOR, true);
		
		threadManager = new ThreadManager(numThreads);
		
		if (inputVcf != null) {
			this.knownVariants = Variant.loadFromFile(inputVcf);
		}
		
		return tempDir.toString();
	}
	
	public List<Variant> getKnownVariants(String chromosome) {
		List<Variant> variants = null;
		if (knownVariants != null) {
			variants = knownVariants.get(chromosome);
		}
		
		if (variants == null) {
			variants = Collections.emptyList();
		}
		
		return variants;
	}
	
	public void setReference(String reference) {
		this.reference = reference;
	}

	public void setAssemblerSettings(AssemblerSettings settings) {
		this.assemblerSettings = settings;
	}
	
	public void setNumThreads(int numThreads) {
		this.numThreads = numThreads;
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
	
	private static String getVersion() {
		String version = "unknown";
		String metaFile = "/META-INF/maven/abra2/abra2/pom.properties";
		Properties prop = new Properties();
		try {
			URL url = NativeLibraryLoader.class.getResource(metaFile);
			InputStream input = url.openStream();
			prop.load(input);
			input.close();
			version = prop.getProperty("version");
		} catch (Exception e) {
			e.printStackTrace();
			Logger.error("Error reading version from pom.properties");
		}
		
		return version;
	}
	
	private static String getCommandLine(String[] args) {
		String jar = "";
		CodeSource cs = Abra.class.getProtectionDomain().getCodeSource();
		if (cs != null) {
			jar = cs.getLocation().toString();
			if (jar.startsWith("file:")) {
				jar = jar.replaceFirst("file:", "");
			}
		}
		
		StringBuffer cl = new StringBuffer();
		cl.append(jar);
		for (String arg : args) {
			cl.append(' ');
			cl.append(arg);
		}
		
		return cl.toString();
	}

	public static void run(String[] args) throws Exception {
		
		String version = getVersion();
		Logger.info("Abra version: " + version);
		
		String cl = getCommandLine(args);
		
		Logger.info("Abra params: [" + cl + "]");
		
		ReAlignerOptions options = new ReAlignerOptions();
		options.parseOptions(args);

		if (options.isValid()) {
			
			Logger.setLevel(options.getLoggerLevel());

			AssemblerSettings assemblerSettings = new AssemblerSettings();

			assemblerSettings.setKmerSize(options.getKmerSizes());
			assemblerSettings.setMinContigLength(options.getMinContigLength());
			assemblerSettings.setMinNodeFrequncy(options.getMinNodeFrequency());
			assemblerSettings.setMinBaseQuality(options.getMinBaseQuality());
			assemblerSettings.setMinReadCandidateFraction(options.getMinReadCandidateFraction());
			assemblerSettings.setMaxAverageDepth(options.getMaxAverageRegionDepth());
			assemblerSettings.setMinEdgeRatio(options.getMinEdgeRatio());
			assemblerSettings.setMaxNodes(options.getMaxNodes());

			ReAligner realigner = new ReAligner();
			realigner.setReference(options.getReference());
			realigner.setRegionsBed(options.getTargetRegionFile());
			realigner.setAssemblerSettings(assemblerSettings);
			realigner.setNumThreads(options.getNumThreads());
			realigner.isPairedEnd = options.isPairedEnd();
			realigner.minMappingQuality = options.getMinimumMappingQuality();
			realigner.maxMismatchRate = options.getMaxMismatchRate();
			realigner.maxReadsInRegion = options.getMaxReadsInRegion();
			realigner.hasPresetKmers = options.hasPresetKmers();
			realigner.isSkipAssembly = options.isSkipAssembly();
			realigner.isSkipUnmappedTrigger = options.isSkipUnmappedAssemblyTrigger();
			realigner.useObservedIndels = options.useObservedIndels();
			realigner.shouldSort = options.shouldSort();
			realigner.maxRealignDist = options.getMaxRealignDist();
			realigner.maxAssembledContigs = options.getMaxAssembledContigs();
			realigner.useConsensusSeq = options.useConsensusSequence();
			realigner.isKeepTmp = options.isKeepTmp();
			realigner.tmpDir = options.getTmpDir();
			realigner.useSoftClippedReads = options.useSoftClippedReads();
			realigner.junctionFile = options.getJunctionFile();
			realigner.gtfJunctionFile = options.getGtfJunctionFile();
			realigner.contigFile = options.getContigFile();
			realigner.swScoring = options.getSmithWatermanScoring();
			realigner.softClipParams = options.getSoftClipParams();
			realigner.maxCachedReads = options.getMaxCachedReads();
			realigner.finalCompressionLevel = options.getCompressionLevel();
			realigner.minAnchorLen = options.getContigAnchor()[0];
			realigner.maxAnchorMismatches = options.getContigAnchor()[1];
			realigner.chromosomesToSkipRegex = options.getChromosomesToSkipRegex();
			realigner.shouldUnsetDuplicates = options.shouldUnsetDuplicates();
			realigner.inputVcf = options.getInputVcf();
			realigner.shouldCreateIndex = options.shouldCreateIndex();
			realigner.shouldUseGkl = options.shouldUseGkl();
			realigner.ambiguousMapq = options.getAmbiguousMapq();
			realigner.maxReadNoise = options.getMaxReadNoise();
			realigner.maxReadsInRamForSort = options.getMaxReadsInRamForSort();
			MAX_REGION_LENGTH = options.getWindowSize();
			MIN_REGION_REMAINDER = options.getWindowOverlap();
			REGION_OVERLAP = options.getWindowOverlap();
			
			realigner.cl = cl.toString();
			realigner.version = version;
			
			long s = System.currentTimeMillis();
			
			realigner.reAlign(options.getInputFiles(), options.getOutputFiles());

			long e = System.currentTimeMillis();

			Logger.info("Elapsed seconds: " + (e - s) / 1000);
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


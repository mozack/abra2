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
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import abra.ReadEvaluator.Alignment;
import abra.SSWAligner.SSWAlignerResult;
import abra.SimpleMapper.Orientation;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;
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
	private boolean isSkipAssembly;
	private boolean isSkipNonAssembly;
	
	// If true, the input target file specifies kmer values
	private boolean hasPresetKmers = false;
	
	// RNA specific
	private String junctionFile;
	private List<Feature> junctions = new ArrayList<Feature>();
	
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
		loadJunctions();
		
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

		// Spawn thread for each chromosome
		// TODO: Validate identical sequence dictionary for each input file
		for (SAMSequenceRecord seqRecord : this.samHeaders[0].getSequenceDictionary().getSequences()) {
			String chromosome = seqRecord.getSequenceName();
			this.spawnChromosomeThread(chromosome);
		}
		
		log("Waiting for all threads to complete");
		threadManager.waitForAllThreadsToComplete();
		
		contigWriter.close();
		
		clock.stopAndPrint();		
		
		for (SAMFileWriter writer : this.writers) {
			writer.close();
		}
		
		System.err.println("Done.");
	}
	
	void debug(SAMRecord read, String msg) {
		if (read.getReadName().equals("UNC9-SN296:440:C5F7CACXX:5:2108:19995:43952")) {
			System.err.println("UNC9-SN296:440:C5F7CACXX:5:2108:19995:43952 - " + msg);
		}
	}
	
	
	
	void processChromosome(String chromosome) throws Exception {
		
		System.err.println("Processing chromosome: " + chromosome);
		
		MultiSamReader reader = new MultiSamReader(this.inputSams, this.minMappingQuality, this.isPairedEnd, chromosome);
		
		List<List<SAMRecordWrapper>> currReads = new ArrayList<List<SAMRecordWrapper>>();
		for (int i=0; i<this.inputSams.length; i++) {
			currReads.add(new ArrayList<SAMRecordWrapper>());
		}
		
		List<List<SAMRecordWrapper>> outOfRegionReads = new ArrayList<List<SAMRecordWrapper>>();
		for (int i=0; i<this.inputSams.length; i++) {
			outOfRegionReads.add(new ArrayList<SAMRecordWrapper>());
		}
		
		Map<Feature, Map<SimpleMapper, SSWAlignerResult>> regionContigs = new HashMap<Feature, Map<SimpleMapper, SSWAlignerResult>>();
		int readCount = 0;
		
		List<Feature> chromosomeRegions = new ArrayList<Feature>();
		for (Feature region : regions) {
			if (region.getSeqname().equals(chromosome)) {
				chromosomeRegions.add(region);
			}
		}
		
		List<Feature> chromosomeJunctions = new ArrayList<Feature>();
		for (Feature junction : junctions) {
			if (junction.getSeqname().equals(chromosome)) {
				chromosomeJunctions.add(junction);
			}
		}
		
//		Map<Feature, List<Feature>> regionJunctions = new HashMap<Feature, List<Feature>>();
//		for (Feature region : chromosomeRegions) {
//			regionJunctions.put(region, new ArrayList<Feature>());
//			for (Feature junction : chromosomeJunctions) {
//				if (region.containsEitherEnd(junction, MAX_REGION_LENGTH)) {
//					regionJunctions.get(region).add(junction);
//				}
//			}
//		}
		
		Map<Feature, List<Feature>> regionJunctions = JunctionUtils.getRegionJunctions(chromosomeRegions, chromosomeJunctions, readLength, MAX_REGION_LENGTH);
		
		Set<Integer> regionsToProcess = new TreeSet<Integer>();
	
		int currRegionIdx = -1;
		
		for (SAMRecordWrapper record : reader) {
			int regionIdx = Feature.findFirstOverlappingRegion(reader.getSAMFileHeader(), record, chromosomeRegions, currRegionIdx);
						
			// Identify next region that is a candidate for processing
			// Note: Splicing can cause reads to go in and out of a region
			if (regionIdx >= 0) {
				regionsToProcess.add(regionIdx);
				
				// Cache read for processing at end of region
				debug(record.getSamRecord(), "Adding read.");
				currReads.get(record.getSampleIdx()).add(record);
			}
			
			Iterator<Integer> regionIter = regionsToProcess.iterator();
			if (regionIter.hasNext()) {
				currRegionIdx = regionIter.next();
				// If start position for current read is beyond current region, trigger assembly
				Feature currRegion = chromosomeRegions.get(currRegionIdx);
				if (record.getAdjustedAlignmentStart() > currRegion.getEnd() + this.readLength*2) {
					System.err.println("Processing region: " + currRegion);
					Map<SimpleMapper, SSWAlignerResult> mappedContigs = processRegion(currRegion, currReads, regionJunctions.get(currRegion));
					System.err.println("Region: " + currRegion + " assembled: " + mappedContigs.keySet().size() + " contigs");
					regionContigs.put(currRegion, mappedContigs);
					// Remove curr region from list of regions to process
					regionIter.remove();
				}
			}
			
			
			// TODO: Consider dropping this...  Reads are out of scope when we've moved beyond them via standard processing?
			if (regionIdx < 0) {
				
				// Process out of region read and output if ready.
				List<SAMRecordWrapper> outOfRegionReadsForSample = outOfRegionReads.get(record.getSampleIdx());
				debug(record.getSamRecord(), "Out of Region!!!");
				outOfRegionReadsForSample.add(record);
				
				if (outOfRegionReads.get(record.getSampleIdx()).size() > 2500) {
					synchronized(this.writers[record.getSampleIdx()]) {
						for (SAMRecordWrapper outOfRegionRead : outOfRegionReadsForSample) {
							this.writers[record.getSampleIdx()].addAlignment(outOfRegionRead.getSamRecord());
						}
					}
					
					outOfRegionReadsForSample.clear();
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
							debug(record.getSamRecord(), "Ready to remap");
							sampleReadsToRemap.add(read);
							iter.remove();
						} else {
							debug(record.getSamRecord(), "Not ready to remap");
						}
					}					
				}

				// Remap out of scope reads
				remapReads(regionContigs, readsToRemap);
				
				// Remove out of scope region assemblies
				List<Feature> regionsToRemove = new ArrayList<Feature>();
				for (Feature region : regionContigs.keySet()) {
					if (getFirstStartPos(currReads)-region.getStart() > MAX_READ_RANGE) {
						regionsToRemove.add(region);
					}
				}
				
				for (Feature region : regionsToRemove) {
					System.err.println("Removing contigs for region: " + region);
					regionContigs.remove(region);
				}

				String logPrefix = record.getSamRecord().getReferenceName() + ":" + record.getSamRecord().getAlignmentStart() + " : ";
				
				if (regionContigs.size() > 10) {
					System.err.println(logPrefix + "regionContigs size: " + regionContigs.size());
				}
				
				int currReadsCount = 0;
				for (List<SAMRecordWrapper> reads : currReads) {
					currReadsCount += reads.size();
				}

				if (currReadsCount > 10000) {
					System.err.println(logPrefix + "Curr reads size: " + currReadsCount);
				}
				
				int outOfRegionCount = 0;
				for (List<SAMRecordWrapper> reads : outOfRegionReads) {
					outOfRegionCount += reads.size();
				}

				if (outOfRegionCount > 10000) {
					System.err.println(logPrefix + "Out of region reads size: " + outOfRegionCount);
				}
			}
			
			readCount += 1;
		}
		
		// Attempt to process last region if applicable
		Iterator<Integer> regionIter = regionsToProcess.iterator();
		if (regionIter.hasNext()) {
			currRegionIdx = regionIter.next();
			
			// We've moved beyond the current region
			// Assemble reads
			Feature region = chromosomeRegions.get(currRegionIdx);
			System.err.println("Processing region: " + region);
			Map<SimpleMapper, SSWAlignerResult> mappedContigs = processRegion(region, currReads, regionJunctions.get(region));
			System.err.println("Region: " + region + " assembled: " + mappedContigs.keySet().size() + " contigs");
			regionContigs.put(region, mappedContigs);
		}
		
		// Remap remaining reads
		remapReads(regionContigs, currReads);
		currReads.clear();
		regionContigs.clear();
		
		// Output remaining out of region reads
		for (int i=0; i<outOfRegionReads.size(); i++) {
			List<SAMRecordWrapper> outOfRegionReadsForSample = outOfRegionReads.get(i);
			synchronized(this.writers[i]) {
				for (SAMRecordWrapper outOfRegionRead : outOfRegionReadsForSample) {
					this.writers[i].addAlignment(outOfRegionRead.getSamRecord());
				}
			}
			
			outOfRegionReadsForSample.clear();
		}
		
		reader.close();
		
		System.err.println("Chromosome: " + chromosome + " done.");
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
		System.err.println("isSkipAssembly: " + isSkipAssembly);
		System.err.println("isSkipNonAssembly: " + isSkipNonAssembly);
		
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
		
	private void spawnChromosomeThread(String chromosome) throws InterruptedException {
		ReAlignerRunnable thread = new ReAlignerRunnable(threadManager, this, chromosome);
		System.err.println("Spawning thread for chromosome: " + chromosome);
		threadManager.spawnThread(thread);
	}
	
	private synchronized void appendContigs(String contigs) throws IOException {
		contigWriter.write(contigs);
	}
	
	private void remapRead(ReadEvaluator readEvaluator, SAMRecord read, int origEditDist) {
		
		
		Alignment alignment = readEvaluator.getImprovedAlignment(origEditDist, read.getReadString(), read);
		if (alignment != null) {
			
			int readPos = alignment.pos;
			
			// Set contig alignment info for all reads that map to contigs (even if read is unchanged)
			String ya = alignment.chromosome + ":" + alignment.contigPos + ":" + alignment.contigCigar;
			
			// If no change to alignment, just record the YA tag
			if (!read.getReadUnmappedFlag() && read.getAlignmentStart() == readPos && read.getCigarString().equals(alignment.cigar)) {
				read.setAttribute("YA", ya);
			}
			
			// If the read has actually moved to an improved alignment, update
			if (origEditDist > alignment.numMismatches && (read.getReadUnmappedFlag() || read.getAlignmentStart() != readPos || !read.getCigarString().equals(alignment.cigar))) {
				
				read.setAttribute("YA", ya);

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
				read.setMappingQuality(Math.min(read.getMappingQuality()+10, 60));
			}
		}
	}
	
	private void remapReads(Map<Feature, Map<SimpleMapper, SSWAlignerResult>> mappedContigs, List<List<SAMRecordWrapper>> readsList) throws Exception {
		
		int numContigs = 0;
		for (Feature region : mappedContigs.keySet()) {
			numContigs += mappedContigs.get(region).size();
		}

		if (readsList.get(0).size() > 0) {
			System.err.println("** REMAPPING [" + readsList.get(0).size() + "] reads to [" + numContigs + "] contigs");
		}
		
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
				
				System.err.println("Read edit dist: " + read.getReadName() + " : " + origEditDist);
				
				if (origEditDist > 0 || SAMRecordUtils.getNumSplices(read) > 0) {
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
	
	private SSWAlignerResult alignContig(String contig, SSWAligner ssw, List<SSWAligner> sswJunctions) {
		SSWAlignerResult bestResult = null;
		int bestScore = -1;
		
		SSWAlignerResult sswResult;
		for (SSWAligner sswJunc : sswJunctions) {
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


		if (bestResult != null) {
			System.err.println("BEST_SSW: " + bestResult.getGenomicPos() + " : " + bestResult.getCigar() + " : " + bestResult.getRefPos() + " : " + bestResult.getScore() + " : " + bestResult.getSequence());
		} else {
			System.err.println("NO_SSW: " + contig);
		}
		
		//TODO: Check for tie scores with different final alignment
		
		return bestResult;
		
//		mappedContigs.put(new SimpleMapper(bestResult.getSequence()), bestResult);
	}
	
	public Map<SimpleMapper, SSWAlignerResult> processRegion(Feature region, List<List<SAMRecordWrapper>> reads, List<Feature> junctions) throws Exception {
		if (isDebug) {
			log("Processing region: " + region.getDescriptor());
		}
		
		if (region.getLength() > 10000) {
			throw new IllegalArgumentException("Region too big: [" + region + "]");
		}
		
		Map<SimpleMapper, SSWAlignerResult> mappedContigs = new HashMap<SimpleMapper, SSWAlignerResult>();
		
		List<List<SAMRecordWrapper>> readsList = subsetReads(region, reads);
		
		try {
			String contigsFasta = tempDir + "/" + region.getDescriptor() + "_contigs.fasta";
			
			List<String> bams = new ArrayList<String>(Arrays.asList(this.inputSams));
			
			// Get reference sequence matching current region (pad by 2 read lengths on each side)
			int chromosomeLength = c2r.getReferenceLength(region.getSeqname());
			int refSeqStart = Math.max((int) region.getStart() - this.readLength*2, 1);
			int refSeqLength = Math.min((int) region.getLength() + this.readLength*4, chromosomeLength-1);
			
			System.out.println("Getting reference for: " + region.getSeqname() + ":" + refSeqStart + ", length: " + refSeqLength);
			
			String refSeq = c2r.getSequence(region.getSeqname(), refSeqStart, refSeqLength);
			
			SSWAligner ssw = new SSWAligner(refSeq, region.getSeqname(), refSeqStart);
			
			List<SSWAligner> sswJunctions = new ArrayList<SSWAligner>();
			
			List<List<Feature>> junctionPermutations = JunctionUtils.combineJunctions(junctions, this.readLength);
			
			System.err.println("NUM_JUNCTION_PERMUTATIONS:\t" + junctionPermutations.size() + "\t" + region);
			
			for (List<Feature> junctionPerm : junctionPermutations) {
				System.err.println("NUM_JUNCTIONS:\t" + junctionPerm.size() + "\t" + region);
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
					String rightSeq = c2r.getSequence(region.getSeqname(), rightStart, rightStop-rightStart);
					juncSeq.append(rightSeq);
					// Junction pos and length should already be added
					
					SSWAligner sswJunc = new SSWAligner(juncSeq.toString(), region.getSeqname(), refStart, junctionPos, junctionLengths);
					sswJunctions.add(sswJunc);
				}
			}
						
			// Assemble contigs
			if (this.isSkipAssembly || region.getKmer() > this.readLength-15) {
				System.err.println("Skipping assembly of region: " + region.getDescriptor() + " - " + region.getKmer());
			} else {
				NativeAssembler assem = (NativeAssembler) newAssembler(region);
				List<Feature> regions = new ArrayList<Feature>();
				regions.add(region); 
				String contigs = assem.assembleContigs(bams, contigsFasta, tempDir, regions, region.getDescriptor(), true, this, c2r, readsList);
				
				if (!contigs.equals("<ERROR>") && !contigs.equals("<REPEAT>") && !contigs.isEmpty()) {
					
					// TODO: Turn this off by default
					appendContigs(contigs);
					
					List<ScoredContig> scoredContigs = ScoredContig.convertAndFilter(contigs);
					
					System.err.println("# SCORED CONTIGS: " + scoredContigs.size());
					
					// Map contigs to reference
					for (ScoredContig contig : scoredContigs) {
						// Filter contigs that match the reference
						if (!refSeq.contains(contig.getContig())) {
							
							SSWAlignerResult sswResult = alignContig(contig.getContig(), ssw, sswJunctions);
							
							if (sswResult != null) {
								// TODO: In multi-region processing, check to ensure identical contigs have identical mappings
								mappedContigs.put(new SimpleMapper(sswResult.getSequence()), sswResult);
							}
						}
					}
				} 
			}
			
			if (!this.isSkipNonAssembly) {
				System.err.println("Processing non-assembled contigs for region: [" + region + "]");
				// Go through artificial contig generation using indels observed in the original reads
				AltContigGenerator altContigGenerator = new AltContigGenerator();
				Collection<String> altContigs = altContigGenerator.getAltContigs(readsList, c2r, readLength);
				
				for (String contig : altContigs) {
					// TODO: Check to see if this contig is already in the map before aligning
					SSWAlignerResult sswResult = ssw.align(contig);
					if (sswResult != null) {
						//TODO: Introduce penalty for non-assembled contigs?
						mappedContigs.put(new SimpleMapper(sswResult.getSequence()), sswResult);
					}
				}
			}
		}
		catch (Exception e) {
			e.printStackTrace();
			throw e;
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
		
	private void loadRegions() throws IOException {
		this.regions = getRegions(regionsBed, readLength, hasPresetKmers);
		
		System.err.println("Num regions: " + regions.size());
		if (isDebug) {
			for (Feature region : regions) {
				System.err.println(region.getSeqname() + "\t" + region.getStart() + "\t" + region.getEnd() + "\t" + region.getKmer());
			}
		}
	}
	
	private void loadJunctions() throws IOException {
		if (this.junctionFile != null) {
			RegionLoader loader = new RegionLoader();
			this.junctions = loader.load(junctionFile, false);
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
					
					cnt += 1;
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
		
		int maxKmerSize = this.readLength-5; 
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
			realigner.isSkipAssembly = options.isSkipAssembly();
			realigner.isSkipNonAssembly = options.isSkipNonAssembly();
			realigner.junctionFile = options.getJunctionFile();

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


/* Copyright 2013 University of North Carolina at Chapel Hill.  All rights reserved. */
package abra;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Random;
import java.util.Set;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.ValidationStringency;

/**
 * Handles regional assembly by invoking the native assembler.
 * 
 * @author Lisle E. Mose (lmose at unc dot edu)
 */
public class NativeAssembler {
	
	//TODO: Calc dynamically
	private static final int MAX_READ_LENGTHS_PER_REGION = 6;
	
	public static final int CYCLE_KMER_LENGTH_THRESHOLD = 43;
	
	private static final int MIN_CANDIDATE_BASE_QUALITY = 10;
	
	private boolean truncateOnRepeat;
	private int maxContigs = 5000;
	private int maxPathsFromRoot;
	private int readLength;
	private int[] kmers;
	private int minKmerFrequency;
	private int minBaseQuality;
	private double minReadCandidateFraction;
	private int maxAverageDepth;
	private boolean isCycleExceedingThresholdDetected = false;
	private double minEdgeRatio;
	private int maxNodes;
	private boolean isSkipUnmappedTrigger = false;
	private int maxReadLength; // Includes merged reads

	private native String assemble(String input, String output, String prefix,
			int truncateOnRepeat, int maxContigs, int maxPathsFromRoot, int readLength, 
			int kmerSize, int minKmerFreq, int minBaseQuality, double minEdgeRatio, int debug,
			int maxNodes);
	
	
	private boolean isHardClipped(SAMRecord read) {
		return read.getCigarString().contains("H");
	}
		
	//
	//  Require at lest <fraction> number of candidate reads per average region depth
	//  inclusive of overlapping reads.
	private int minCandidateCount(int numReads, Feature region) {
		double fraction = minReadCandidateFraction;
		
		int minCount = (int) ((double) numReads / readLengthsPerRegion(region) * fraction);
		
		// Always require at least 2 candidate reads.
		return Math.max(minCount, 2);
	}
	
	// Calc number of read lengths per region, padding by 2 to account for reads overlapping region ends.
	private double readLengthsPerRegion(Feature region) {
		return (double) region.getLength() / (double) readLength + 2;
	}
	
	private double readLengthsForAllRegions(List<Feature> regions) {
		double lengths = 0;
		for (Feature region : regions) {
			lengths += readLengthsPerRegion(region);
		}
		return lengths;
	}
	
	//
	//  Calc desired number of reads per file.
	private int desiredNumberOfReads(List<Feature> regions) {
		return (int) (readLengthsForAllRegions(regions) * (double) maxAverageDepth); 
	}
	
	private int countGaps(Cigar cigar) {
		int count = 0;
		for (CigarElement elem : cigar.getCigarElements()) {
			if (elem.getOperator() == CigarOperator.D || elem.getOperator() == CigarOperator.I || elem.getOperator() == CigarOperator.N) {
				count += 1;
			}
		}
		
		return count;
	}
	
	private int firstIndelLength(Cigar cigar) {
		int len = 0;
		
		for (CigarElement elem : cigar.getCigarElements()) {
			if (elem.getOperator() == CigarOperator.D || elem.getOperator() == CigarOperator.I) {
				len = elem.getLength();
				break;
			}
		}
		
		return len;
	}
	
	// Get length of longest soft clip element that also overlaps the current region.
	private int maxSoftClipLength(SAMRecord read, Feature region) {
		int len = 0;
		if (read.getCigarLength() > 1) {
			
			CigarElement elem = read.getCigar().getCigarElement(0);
			if (elem.getOperator() == CigarOperator.S && read.getAlignmentStart() >= region.getStart()-readLength) {
				len = elem.getLength();
			}
			
			elem = read.getCigar().getCigarElement(read.getCigarLength()-1);
			if (elem.getOperator() == CigarOperator.S && elem.getLength() > len && read.getAlignmentEnd() <= region.getEnd()+readLength) {
				len = elem.getLength();
			}			
		}
		
		return len;
	}
	
	private int getInsertBases(SAMRecord read) {
		int numInsertBases = 0;
		
		for (CigarElement element : read.getCigar().getCigarElements()) {
			if (element.getOperator() == CigarOperator.I) {
				numInsertBases += element.getLength();
			}
		}
		
		return numInsertBases;	
	}
	
	private boolean isAssemblyTriggerCandidate(SAMRecord read, CompareToReference2 c2r, Feature region) {
		
		// High quality unmapped read anchored by mate
		if (!isSkipUnmappedTrigger && read.getReadUnmappedFlag() && !read.getMateUnmappedFlag() &&
			read.getReadLength() >= readLength * .9 && 
			SAMRecordUtils.getNumHighQualBases(read, MIN_CANDIDATE_BASE_QUALITY) >= readLength * .9) {
			return true;
		}
		
		int numGaps = countGaps(read.getCigar());
		
		int insertBases = getInsertBases(read);
		
		// More than one indel and/or splice in read
//		if (numGaps > 1) {
//			return true;
//		}
		
//		if (numGaps > 0) {
//			int qualAdjustedEditDist = c2r.numHighQualityMismatches(read, MIN_CANDIDATE_BASE_QUALITY) + SAMRecordUtils.getNumIndelBases(read);
//			if (qualAdjustedEditDist > (readLength * .25)) {
//				return true;
//			}			
//		}
		
		if (insertBases > readLength * .15) {
			return true;
		}
		
		if (maxSoftClipLength(read, region) > readLength * .25) {
			if (SAMRecordUtils.getNumHighQualBases(read, MIN_CANDIDATE_BASE_QUALITY) >= readLength * .9) {
				return true;
			}
		}
		
		// Increment candidate count for substantial high quality soft clipping
		// TODO: Higher base quality threshold?
//		if (read.getCigarString().contains("S") && (c2r.numHighQualityMismatches(read, MIN_CANDIDATE_BASE_QUALITY) > (readLength/10))) {
//			return true;
//		}


		
		// if indel is at least 10% of read length
		// TODO: Longer threshold (especially for deletions)?
//		if (numGaps == 1 && firstIndelLength(read.getCigar()) >= (readLength/10)) {
//			return true;
//		}
		
		// Read contains indel / intron and SNV (inclusive of soft clip mismatch)
		// TODO: Higher base qual threshold?
		// TODO: Minimum edit dist (1 base indel + 1 mismatch insufficient)
//		if (numGaps > 0 && c2r.numHighQualityMismatches(read, MIN_CANDIDATE_BASE_QUALITY) > 0) {
//			return true;
//		}
		
		// Increment candidate count for indels
//		if (read.getCigarString().contains("I") || read.getCigarString().contains("D") || read.getCigarString().contains("N")) {
//			isCandidate = true;
//		}
		
		// Increment candidate count if read contains at least 3 high quality mismatches
//		if (SAMRecordUtils.getIntAttribute(read, "NM") >= 3) {
//			if (c2r.numHighQualityMismatches(read, MIN_CANDIDATE_BASE_QUALITY) > 3) {
//				isCandidate = true;
//			}
//		}

		return false;
	}
	
	public String assembleContigs(List<String> inputFiles, List<Feature> regions, String prefix,
			boolean checkForDupes, ReAligner realigner, CompareToReference2 c2r, List<List<SAMRecordWrapper>> readsList,
			int mnf, int mbq, double mer, StringBuffer readBuffer) {
		
				
		String contigs = "";
		
		long start = System.currentTimeMillis();
		
		int readCount = 0;

		boolean isAssemblyCandidate = false;
		
		try {
			
			int[] unfilteredReads = new int[readsList.size()];
			
			int sampleIdx = 0;
			
			maxReadLength = readLength;
			
			for (List<SAMRecordWrapper> reads : readsList) {
				int candidateReadCount = 0;
			
				for (SAMRecordWrapper read : reads) {
					
					if (read.shouldAssemble()) {
						if (!isAssemblyCandidate && isAssemblyTriggerCandidate(read.getSamRecord(), c2r, regions.get(0))) {
//							System.err.println("trigger: " + sampleIdx + " : " + read.getSamRecord().getSAMString());
							candidateReadCount++;
						}
						unfilteredReads[sampleIdx] += 1;
						
						maxReadLength = Math.max(maxReadLength, read.getReadLength());
					}
				}
				
				if (candidateReadCount >= minCandidateCount(unfilteredReads[sampleIdx], regions.get(0))) {
					isAssemblyCandidate = true;
				}
				
				readCount += unfilteredReads[sampleIdx];
				sampleIdx += 1;
			}
			
//			StringBuffer readBuffer = new StringBuffer();
			
			if (isAssemblyCandidate) {
				
				Logger.debug("ASSEMBLY_TRIGGERED\t%s", regions.get(0));
				
				if ((kmers.length == 0) || (kmers[0] < KmerSizeEvaluator.MIN_KMER)) {
					KmerSizeEvaluator kmerEval = new KmerSizeEvaluator();
					int kmer = kmerEval.identifyMinKmer(maxReadLength, c2r, regions);
					// Cap max kmer size below buffer size in C code
					this.kmers = realigner.toKmerArray(kmer, Math.min(maxReadLength,199));
				}
				
				int downsampleTarget = desiredNumberOfReads(regions);
				
				char sampleId = 1;
				
				sampleIdx = 0;
				for (List<SAMRecordWrapper> reads : readsList) {
					// Default to always keep
					double keepProbability = 1.1;
					
					if (reads.size() > downsampleTarget) {
						keepProbability = (double) downsampleTarget / (double) unfilteredReads[sampleIdx];
					}
					
					Set<String> mergedReadIds = new HashSet<String>();
					
					Random random = new Random(1);
					//Random random = new Random(System.currentTimeMillis());
					
					for (SAMRecordWrapper readWrapper : reads) {
						
						SAMRecord read = readWrapper.getSamRecord();
						
						if (readWrapper.hasMergedSeq() && mergedReadIds.contains(read.getReadName())) {
							// Only include merged sequence once
							continue;
						}
						
						mergedReadIds.add(read.getReadName());
												
						if (readWrapper.shouldAssemble() && random.nextDouble() < keepProbability) {
							
							readBuffer.append(sampleId);
							readBuffer.append(read.getReadNegativeStrandFlag() ? "1" : "0");
							
							String seq;
							String qual;
							
							if (readWrapper.getReadLength() == maxReadLength) {
								seq = readWrapper.getSeq();
								qual = readWrapper.getQual();
								
//								readBuffer.append(readWrapper.getSeq());
//								readBuffer.append(readWrapper.getQual());
							} else {
								StringBuffer basePadding = new StringBuffer();
								StringBuffer qualPadding = new StringBuffer();
								
								for (int i=0; i<maxReadLength-readWrapper.getReadLength(); i++) {
									basePadding.append('N');
									qualPadding.append('!');
								}

								seq = readWrapper.getSeq() + basePadding.toString();
								qual = readWrapper.getQual() + qualPadding.toString();
//								readBuffer.append(readWrapper.getSeq() + basePadding.toString());
//								readBuffer.append(readWrapper.getQual() + qualPadding.toString());							
							}
							
							if (seq.length() != maxReadLength) {
								String msg = String.format("Invalid seq length [%d] for region [%s] read [%s] seq [%s]", seq.length(), regions.get(0), read.getReadName(), seq);
								Logger.error(msg);
								throw new RuntimeException(msg);
							}
							
							if (qual.length() != maxReadLength) {
								String msg = String.format("Invalid qual length [%d] for region [%s] read [%s] qual [%s]", qual.length(), regions.get(0), read.getReadName(), qual);
								Logger.error(msg);
								throw new RuntimeException(msg);
							}
							
							readBuffer.append(seq);
							readBuffer.append(qual);
						}
					}
					
					sampleIdx += 1;
					sampleId += 1;
				}
			}
			
			if (isAssemblyCandidate) {
				for (int kmer : kmers) { 
				
					//TODO: Not really an output file anymore.  Cleanup.
					String outputFile = prefix + "_k" + kmer;
					
//					System.out.println(readBuffer.toString());
					
					contigs = assemble(
							readBuffer.toString(),
							outputFile, 
							prefix, 
							truncateOnRepeat ? 1 : 0,
							maxContigs,
							maxPathsFromRoot,
							maxReadLength,
							kmer,
							Math.max(mnf, 1),
							Math.max(mbq, 2),
							Math.max(mer, .0001),
							Logger.LEVEL == Logger.Level.DEBUG || Logger.LEVEL == Logger.Level.TRACE ? 1 : 0,
							maxNodes);
					
					if (!contigs.equals("<REPEAT>")) {
						break;
					} else {
						if (kmer >= readLength/2 || kmer >= CYCLE_KMER_LENGTH_THRESHOLD) {
							isCycleExceedingThresholdDetected = true;
						}
					}
				}
			} else {
//				System.out.println("Skipping assembly for: " + prefix);
			}

		} catch (Exception e) {
			e.printStackTrace();
			throw new RuntimeException(e);
		}
				
		long end = System.currentTimeMillis();
		
		int kmer = readLength + 1;
		if (kmers.length > 0) {
			kmer = kmers[0];
		}
		
		Logger.debug("Elapsed_msecs_in_NativeAssembler\tRegion:\t%s\tLength:\t%d\tReadCount:\t%d\tElapsed\t%d\tAssembled\t%s\t%d",
				regions.get(0).getDescriptor(), regions.get(0).getLength(), readCount, (end-start), isAssemblyCandidate, kmer);
		
		return contigs;
	}
	
	String nativeAssemble(String input, String output, String prefix, int truncateOnRepeat, int maxContigs, int maxPathsFromRoot, int readLength, int[] kmers,
			int minKmerFreq, int minBaseQuality, double minEdgeRatio, int debug, int maxNodes) {
		String result = "";
		for (int kmer : kmers) {
			result = assemble(input, output, prefix, truncateOnRepeat, maxContigs, maxPathsFromRoot, readLength, kmer, minKmerFreq, minBaseQuality, minEdgeRatio, debug,
					maxNodes);
			if (!result.equals("<REPEAT>")) {
				break;
			}
		}
		return result;
	}
	
	private boolean hasLowQualityBase(SAMRecord read) {
		//TODO: Don't hardcode phred33
		for (int i=0; i<read.getBaseQualityString().length(); i++) {
			if ((read.getBaseQualityString().charAt(i) - '!') < 20) {
				return true;
			}
		}
		
		return false;
	}

	public boolean isTruncateOnRepeat() {
		return truncateOnRepeat;
	}

	public void setTruncateOutputOnRepeat(boolean truncateOnRepeat) {
		this.truncateOnRepeat = truncateOnRepeat;
	}

	public int getMaxPathsFromRoot() {
		return maxPathsFromRoot;
	}

	public void setMaxPathsFromRoot(int maxPathsFromRoot) {
		this.maxPathsFromRoot = maxPathsFromRoot;
	}
	
	public void setReadLength(int readLength) {
		this.readLength = readLength;
	}
	
	public void setKmer(int[] kmers) {
		this.kmers = kmers;
	}
	
	public void setMinKmerFrequency(int frequency) {
		this.minKmerFrequency = frequency;
	}
	
	public void setMinEdgeRatio(double minEdgeRatio) {
		this.minEdgeRatio = minEdgeRatio;
	}
	
	public void setMinBaseQuality(int minBaseQuality) {
		this.minBaseQuality = minBaseQuality;
	}
	
	public void setMaxNodes(int maxNodes) {
		this.maxNodes = maxNodes;
	}
	
	public void setMaxAverageDepth(int maxAverageDepth) {
		this.maxAverageDepth = maxAverageDepth;
	}
	
	public void setMinReadCandidateFraction(double minReadCandidateFraction) {
		this.minReadCandidateFraction = minReadCandidateFraction;
	}
	
	public void setSkipUnmappedTrigger(boolean shouldSkip) {
		this.isSkipUnmappedTrigger = shouldSkip;
	}
		
	public boolean isCycleExceedingThresholdDetected() {
		return isCycleExceedingThresholdDetected;
	}
		
	static class Position implements Comparable<Position> {
		private String chromosome;
		private int position;
		
		Position(String chromosome, int position) {
			this.chromosome = chromosome;
			this.position = position;
		}
		
		String getChromosome() {
			return chromosome;
		}
		
		int getPosition() {
			return position;
		}
		
		public String toString() {
			return chromosome + ":" + position;
		}

		@Override
		public int compareTo(Position that) {
			int compare = this.chromosome.compareTo(that.chromosome);
			if (compare == 0) {
				compare = this.position - that.position;
			}
			return compare;
		}
	}
	
	public static void main(String[] args) throws Exception {
		
		/*
		NativeLibraryLoader l = new NativeLibraryLoader();
		l.load("/home/lmose/code/abra/target");
		
		NativeAssembler assem = new NativeAssembler();
		assem.setTruncateOutputOnRepeat(true);
		assem.setMaxContigs(5000);
		assem.setMaxPathsFromRoot(5000);
		assem.setKmer(new int[] { 43 });
		assem.setReadLength(100);
		assem.setMinKmerFrequency(2);
		assem.setMaxAverageDepth(400);
		assem.setShouldSearchForSv(true);
		
//		String bam1 = args[0];
		String bam1 = "/home/lmose/dev/abra/sv/test.bam";
		List<String> inputFiles = new ArrayList<String>();
		inputFiles.add(bam1);
		
		//String output = args[2];
		String output = "/home/lmose/dev/abra/sv/output.txt";
		//chr18:60,793,358-60,793,758
		Feature region = new Feature("chr18", 60793358, 60793758);
		String prefix = "pre";
		boolean checkForDupes = true;
		ReAligner realigner = new ReAligner();
		CompareToReference2 c2r = new CompareToReference2();
		//c2r.init(args[3]);
		c2r.init("/home/lmose/reference/chr18/chr18.fa");
		
		List<Feature> regions = new ArrayList<Feature>();
		regions.add(region);
		String contigs = assem.assembleContigs(inputFiles, output, "asm_temp", regions, prefix, checkForDupes, realigner, c2r);
		System.err.println(contigs);
		
		System.err.println("-------------------------");
		
		List<BreakpointCandidate> svCandidates = assem.getSvCandidateRegions();
		for (BreakpointCandidate svCandidate : svCandidates) { 
			System.err.println("SV: " + region.getDescriptor() + "-->" + svCandidate.getRegion().getDescriptor());
		}
		*/
		
//		assem.assembleContigs(args[0], args[1], "contig");
		
//		for (int i=0; i<10; i++) {
//			run(args[0], args[1] + "_" + i);
//		}
		
//		run(args[0], args[1]);
		
//		assem.assembleContigs("/home/lmose/code/abra/src/main/c/1810_reads.txt",
//				"/home/lmose/code/abra/src/main/c/1810.fa", "bar");
	}

	
}

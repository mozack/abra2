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

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMFileReader.ValidationStringency;

/**
 * Handles regional assembly by invoking the native assembler.
 * 
 * @author Lisle E. Mose (lmose at unc dot edu)
 */
public class NativeAssembler {
	
	private boolean truncateOnRepeat;
	private int maxContigs;
	private int maxPathsFromRoot;
	private int readLength;
	private int[] kmers;
	private int minKmerFrequency;
	private int minBaseQuality;
	private double minReadCandidateFraction;
	private Set<String> readIds;
	private int maxAverageDepth;
	List<Position> svCandidates = new ArrayList<Position>();
	List<Feature> svCandidateRegions = new ArrayList<Feature>();
	private boolean shouldSearchForSv = false;

	private native String assemble(String input, String output, String prefix, int truncateOnRepeat, int maxContigs, int maxPathsFromRoot, int readLength, int kmerSize, int minKmerFreq, int minBaseQuality);
	
	private String getIdentifier(SAMRecord read) {
		String id = read.getReadName();
		
		if (read.getReadPairedFlag() && read.getSecondOfPairFlag()) {
			id += "_2";
		}
		
		return id;
	}
	
	private boolean isHardClipped(SAMRecord read) {
		return read.getCigarString().contains("H");
	}
	
	private void filterPositionList(List<Integer> positions, int currentPos) {
		Iterator<Integer> iter = positions.iterator();
		while (iter.hasNext()) {
			if (iter.next() < currentPos-readLength) {
				iter.remove();
			}
		}
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
	
	public String assembleContigs(List<String> inputFiles, String output, String tempDir, List<Feature> regions, String prefix,
			boolean checkForDupes, ReAligner realigner, CompareToReference2 c2r) {
		
		String contigs = "";
		
		long start = System.currentTimeMillis();
		
		int readCount = 0;
		
		int minReadCount = Integer.MAX_VALUE;
		
		try {
			
			// if c2r is null, this is the unaligned region.
			boolean isAssemblyCandidate = c2r == null ? true : false;
			
			List<List<SAMRecord>> readsList = new ArrayList<List<SAMRecord>>();

			for (String input : inputFiles) {
				readIds = new HashSet<String>();
				List<SAMRecord> reads = new ArrayList<SAMRecord>();
				readsList.add(reads);
				
				
				for (Feature region : regions) {
					SAMFileReader reader = new SAMFileReader(new File(input));
					reader.setValidationStringency(ValidationStringency.SILENT);
					
					Iterator<SAMRecord> iter;
					if (region != null) {
						iter = reader.queryOverlapping(region.getSeqname(), (int) region.getStart(), (int) region.getEnd());
					} else {
						iter = reader.iterator();
					}
					
					int candidateReadCount = 0;
					
					while (iter.hasNext()) {
						
						SAMRecord read = iter.next();
						readCount++;
						
						if (read.getReadLength() > readLength) {
							reader.close();
							throw new IllegalArgumentException(
									"Read length exceeds expected value of: " + readLength + " for read [" +
									read.getSAMString() + "]");
						}
						
						// Don't allow same read to be counted twice.
						if ( (!realigner.isFiltered(read)) && 
							 (!read.getDuplicateReadFlag()) && 
							 (!read.getReadFailsVendorQualityCheckFlag()) &&
							 //(Sam2Fastq.isPrimary(read)) &&
							 (!isHardClipped(read)) &&
							 ((!checkForDupes) || (!readIds.contains(getIdentifier(read))))) {
	
							Integer numBestHits = (Integer) read.getIntegerAttribute("X0");
							boolean hasAmbiguousInitialAlignment = numBestHits != null && numBestHits > 1;
							
							if (!hasAmbiguousInitialAlignment) {
								if (!checkForDupes) {
									readIds.add(getIdentifier(read));
								}
								
								// Increment candidate count for indels
								if (!isAssemblyCandidate && read.getCigarString().contains("I") || read.getCigarString().contains("D")) {
									candidateReadCount++;
								}
	
								// Increment candidate count for substantial high quality soft clipping
								// TODO: Check for chimera directly?
								if (!isAssemblyCandidate && (read.getCigarString().contains("S"))) {
									if (c2r.numHighQualityMismatches(read, minBaseQuality) > (readLength/10)) {
										candidateReadCount++;
									}
								}
								
								// Increment candidate count if read contains at least 3 high quality mismatches
								if (!isAssemblyCandidate && (SAMRecordUtils.getIntAttribute(read, "NM") >= 3)) {
									if (c2r.numHighQualityMismatches(read, minBaseQuality) > 3) {
										candidateReadCount++;
									}
								}
								
								reads.add(read);
								
								if (shouldSearchForSv && isSvCandidate(read, region)) {
									svCandidates.add(new Position(read.getMateReferenceName(), read.getMateAlignmentStart()));
								}
							}
						}
					}
					
					reader.close();
					
					if (candidateReadCount > minCandidateCount(reads.size(), region)) {
						isAssemblyCandidate = true;
					}
				}
				
				if (reads.size() < minReadCount) {
					minReadCount = reads.size();
				}
			}
			
			readIds = null;
			
			StringBuffer readBuffer = new StringBuffer();
			
			if (isAssemblyCandidate) {
				
				int downsampleTarget = desiredNumberOfReads(regions);
				
				for (List<SAMRecord> reads : readsList) {
					// Default to always keep
					double keepProbability = 1.1;
					
					if (reads.size() > downsampleTarget) {
						keepProbability = (double) downsampleTarget / (double) reads.size();
					}
					
					Random random = new Random(1);
					
					for (SAMRecord read : reads) {
						if (random.nextDouble() < keepProbability) {
							readBuffer.append(read.getReadNegativeStrandFlag() ? "1" : "0");
							
							if (read.getReadLength() == readLength) {
								readBuffer.append(read.getReadString());
								readBuffer.append(read.getBaseQualityString());
							} else {
								StringBuffer basePadding = new StringBuffer();
								StringBuffer qualPadding = new StringBuffer();
								
								for (int i=0; i<readLength-read.getReadLength(); i++) {
									basePadding.append('N');
									qualPadding.append('!');
								}
								
								readBuffer.append(read.getReadString() + basePadding.toString());
								readBuffer.append(read.getBaseQualityString() + qualPadding.toString());							
							}
						}
					}
					
					// Make this set of reads eligible for GC
					reads.clear();
				}
			}
			
			readsList.clear();
			
			long end1 = System.currentTimeMillis();
			
			System.out.println("Elapsed msecs collection data to assemble" + (end1-start));
			
			if (isAssemblyCandidate) {
				for (int kmer : kmers) { 
				
					String outputFile = output + "_k" + kmer;
					
					contigs = assemble(
							readBuffer.toString(),
							outputFile, 
							prefix, 
							truncateOnRepeat ? 1 : 0,
							maxContigs,
							maxPathsFromRoot,
							readLength,
							kmer,
							minKmerFrequency,
							minBaseQuality);
					
					if (!contigs.equals("<REPEAT>")) {
						break;
					}
				}
			} else {
				System.out.println("Skipping assembly for: " + prefix);
			}

		} catch (Exception e) {
			e.printStackTrace();
			throw new RuntimeException(e);
		}
		
		if (this.shouldSearchForSv) {
			
			Collections.sort(this.svCandidates);
			Position last = null;
			String currentFeatureChr = null;
			int currentFeatureStart = -1;
			int currentFeatureStop = -1;
			int currentFeatureCount = 0;
			
			// TODO: Calc this dynamically
			int windowSize = 500;
			
			for (Position pos : this.svCandidates) {
				if ((last != null) && pos.getChromosome().equals(last.getChromosome()) && 
					 Math.abs(pos.getPosition()-last.getPosition()) < windowSize) {
					
					if (currentFeatureChr == null) {
						currentFeatureChr = pos.getChromosome();
						currentFeatureStart = last.getPosition();
						currentFeatureStop = pos.getPosition() + readLength;
						currentFeatureCount = 1;
					} else {
						currentFeatureStop = pos.getPosition() + readLength;
						currentFeatureCount++;
					}
				} else {
					if (currentFeatureChr != null) {
						if (currentFeatureCount > minReadCount * .02) {
							this.svCandidateRegions.add(new Feature(currentFeatureChr, currentFeatureStart-readLength, currentFeatureStop+readLength));
						}
						currentFeatureChr = null;
						currentFeatureStart = -1;
						currentFeatureStop = -1;
						currentFeatureCount = 0;
					} else {
						currentFeatureChr = pos.getChromosome();
						currentFeatureStart = pos.getPosition();
						currentFeatureStop = pos.getPosition() + readLength;
						currentFeatureCount = 1;
					}
				}
				last = pos;
			}
		}
		
		long end = System.currentTimeMillis();
		
		System.out.println("Elapsed_msecs_in_NativeAssembler\tRegion:\t" + regions.get(0).getDescriptor() + "\tLength:\t" + regions.get(0).getLength() + "\tReadCount:\t" + readCount + "\tElapsed\t" + (end-start));
		
		return contigs;
	}
	
	private boolean isSvCandidate(SAMRecord read, Feature region) {
		boolean isCandidate = false;
		if (!read.getProperPairFlag() && !read.getMateUnmappedFlag()) {
			if (!read.getReferenceName().equals(read.getMateReferenceName())) {
				isCandidate = true;
			} else if (Math.abs(read.getAlignmentStart() - read.getMateAlignmentStart()) > CombineChimera3.MAX_GAP_LENGTH &&
					!region.overlaps(read.getMateReferenceName(), read.getMateAlignmentStart()-(int)region.getLength(), read.getMateAlignmentStart()+read.getReadLength()+(int)region.getLength())) {
				isCandidate = true;
			}
		}
		return isCandidate;
	}
	
	public boolean shouldSearchForSv() {
		return shouldSearchForSv;
	}

	public void setShouldSearchForSv(boolean shouldSearchForSv) {
		this.shouldSearchForSv = shouldSearchForSv;
	}
	
	public List<Feature> getSvCandidateRegions() {
		return this.svCandidateRegions;
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

	public int getMaxContigs() {
		return maxContigs;
	}

	public void setMaxContigs(int maxContigs) {
		this.maxContigs = maxContigs;
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
	
	public void setMinBaseQuality(int minBaseQuality) {
		this.minBaseQuality = minBaseQuality;
	}
	
	public void setMaxAverageDepth(int maxAverageDepth) {
		this.maxAverageDepth = maxAverageDepth;
	}
	
	public void setMinReadCandidateFraction(double minReadCandidateFraction) {
		this.minReadCandidateFraction = minReadCandidateFraction;
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
		System.out.println(contigs);
		
		System.out.println("-------------------------");
		
		List<Feature> svCandidates = assem.getSvCandidateRegions();
		for (Feature svCandidate : svCandidates) {
			System.out.println("SV: " + region.getDescriptor() + "-->" + svCandidate.getDescriptor());
		}
		
//		assem.assembleContigs(args[0], args[1], "contig");
		
//		for (int i=0; i<10; i++) {
//			run(args[0], args[1] + "_" + i);
//		}
		
//		run(args[0], args[1]);
		
//		assem.assembleContigs("/home/lmose/code/abra/src/main/c/1810_reads.txt",
//				"/home/lmose/code/abra/src/main/c/1810.fa", "bar");
	}

	
}

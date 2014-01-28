/* Copyright 2013 University of North Carolina at Chapel Hill.  All rights reserved. */
package abra;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
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
public class NativeAssembler implements Assembler {
	
	private boolean truncateOnRepeat;
	private int maxContigs;
	private int maxPathsFromRoot;
	private int readLength;
	private int[] kmers;
	private int minKmerFrequency;
	private int minBaseQuality;
	private Set<String> readIds;

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
	
	public String assembleContigs(List<String> inputFiles, String output, String tempDir, Feature region, String prefix,
			boolean checkForDupes, ReAligner realigner, CompareToReference2 c2r) {
		
		String contigs = "";
		
		long start = System.currentTimeMillis();
		
		List<String> outputFiles = new ArrayList<String>();
		
		int count = 0;
		
		int readCount = 0;
		
		try {
			
			// if c2r is null, this is the unaligned region.
			boolean isAssemblyCandidate = c2r == null ? true : false;
			
//			int indelReadCount = 0;
			List<List<SAMRecord>> readsList = new ArrayList<List<SAMRecord>>();

			for (String input : inputFiles) {
				readIds = new HashSet<String>();
				List<SAMRecord> reads = new ArrayList<SAMRecord>();
				readsList.add(reads);
				
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
	//					boolean hasAmbiguousBases = read.getReadString().contains("N");
						Integer numBestHits = (Integer) read.getIntegerAttribute("X0");
						boolean hasAmbiguousInitialAlignment = numBestHits != null && numBestHits > 1;
						//TODO: Stampy ambiguous read (mapq < 4)
						
	//					if (!hasAmbiguousBases && !hasAmbiguousInitialAlignment && !hasLowQualityBase(read)) {
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
							// TODO: # mismatches should be function of read length
							if (!isAssemblyCandidate && (SAMRecordUtils.getIntAttribute(read, "NM") >= 3)) {
								if (c2r.numHighQualityMismatches(read, minBaseQuality) > 3) {
									candidateReadCount++;
								}
							}
							
							reads.add(read);
							
						}
					}
				}
				
				reader.close();
				
				//TODO: Calc read fraction based upon read length & region size
				if ((candidateReadCount > reads.size() * .01 / 4.0) && (candidateReadCount >= 2)) {
					isAssemblyCandidate = true;
				}
			}
			
			readIds = null;
			
			StringBuffer readBuffer = new StringBuffer();
			
			if (isAssemblyCandidate) {
				
				for (List<SAMRecord> reads : readsList) {
					// Default to always keep
					double keepProbability = 1.1;
					
					//TODO : Change constant to be a function of read length & region size.
					if (reads.size() > 600) {
						keepProbability = (double) 600 / (double) reads.size();
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
		
		long end = System.currentTimeMillis();
		
		System.out.println("Elapsed_msecs_in_NativeAssembler\tRegion:\t" + region.getDescriptor() + "\tLength:\t" + region.getLength() + "\tReadCount:\t" + readCount + "\tElapsed\t" + (end-start));
		
		return contigs;
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
	
	public static void main(String[] args) throws Exception {
		NativeAssembler assem = new NativeAssembler();
		assem.setTruncateOutputOnRepeat(true);
		assem.setMaxContigs(50000);
		assem.setMaxPathsFromRoot(100000);
		assem.setKmer(new int[] { 43 });
		assem.setReadLength(100);
		assem.setMinKmerFrequency(2);
		
		String bam1 = args[0];
		String bam2 = args[1];
		List<String> inputFiles = new ArrayList<String>();
		inputFiles.add(bam1);
		inputFiles.add(bam2);
		
		String output = args[2];
		
		Feature region = new Feature("chr1", 6162053, 6162453);
		String prefix = "pre";
		boolean checkForDupes = true;
		ReAligner realigner = new ReAligner();
		CompareToReference2 c2r = new CompareToReference2();
		c2r.init(args[3]);
		
		assem.assembleContigs(inputFiles, output, "asm_temp", region, prefix, checkForDupes, realigner, c2r);
		
//		assem.assembleContigs(args[0], args[1], "contig");
		
//		for (int i=0; i<10; i++) {
//			run(args[0], args[1] + "_" + i);
//		}
		
//		run(args[0], args[1]);
		
//		assem.assembleContigs("/home/lmose/code/abra/src/main/c/1810_reads.txt",
//				"/home/lmose/code/abra/src/main/c/1810.fa", "bar");
	}

	
}

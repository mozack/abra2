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

	private native int assemble(String input, String output, String prefix, int truncateOnRepeat, int maxContigs, int maxPathsFromRoot, int readLength, int kmerSize, int minKmerFreq, int minBaseQuality);
	
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
	
	public List<String> assembleContigs(List<String> inputFiles, String output, String tempDir, Feature region, String prefix, boolean checkForDupes, ReAligner realigner) {
		
		long start = System.currentTimeMillis();
		
		List<String> outputFiles = new ArrayList<String>();
		
		int count = 0;
		
		readIds = new HashSet<String>();
		
		try {
			
			String readFile;
			if (region != null) {
				readFile = tempDir + "/" + region.getDescriptor() + ".reads";
			} else {
				readFile = tempDir + "/" + "unaligned.reads";
			}
			
			BufferedWriter writer = new BufferedWriter(new FileWriter(readFile, false));

			for (String input : inputFiles) {
				SAMFileReader reader = new SAMFileReader(new File(input));
				reader.setValidationStringency(ValidationStringency.SILENT);
	
				Iterator<SAMRecord> iter;
				if (region != null) {
					iter = reader.queryOverlapping(region.getSeqname(), (int) region.getStart(), (int) region.getEnd());
				} else {
					iter = reader.iterator();
				}
				
				while (iter.hasNext()) {
					
					SAMRecord read = iter.next();
					
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
							
							writer.write(read.getReadNegativeStrandFlag() ? "1\n" : "0\n");
							
							if (read.getReadLength() == readLength) {
								writer.write(read.getReadString() + "\n");
								writer.write(read.getBaseQualityString() + "\n");
							} else {
								StringBuffer basePadding = new StringBuffer();
								StringBuffer qualPadding = new StringBuffer();
								
								for (int i=0; i<readLength-read.getReadLength(); i++) {
									basePadding.append('N');
									qualPadding.append('!');
								}
								
//								writer.write(read.getReadNegativeStrandFlag() ? "1\n" : "0\n");
								writer.write(read.getReadString() + basePadding.toString() + "\n");
								writer.write(read.getBaseQualityString() + qualPadding.toString() + "\n");							
							}
						}
					}
				}
				
				reader.close();
			}
			
			readIds = null;
			writer.close();
			
			for (int kmer : kmers) { 
			
				String outputFile = output + "_k" + kmer;
				
				count = assemble(
						readFile,
						outputFile, 
						prefix, 
						truncateOnRepeat ? 1 : 0,
						maxContigs,
						maxPathsFromRoot,
						readLength,
						kmer,
						minKmerFrequency,
						minBaseQuality);
				
				if (count > 0) {
					outputFiles.add(outputFile);
					break;
				} else {
					File fileToDelete = new File(outputFile);
					fileToDelete.delete();
				}
			}
			
			File inputReadFile = new File(readFile);
			inputReadFile.delete();

		} catch (IOException e) {
			e.printStackTrace();
			throw new RuntimeException(e);
		}
		
		long end = System.currentTimeMillis();
		
		System.out.println("Elapsed msecs in NativeAssembler: " + (end-start));
		
		return outputFiles;
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
	
//	public static void run(String input, String output) {
//		NativeAssembler assem = new NativeAssembler();
//		assem.setTruncateOutputOnRepeat(false);
//		assem.setMaxContigs(2000000);
//		assem.setMaxPathsFromRoot(5000);
//		
//		assem.assembleContigs(input, output, "contig", true);
//		
//	}
	
	/*
	public static void main(String[] args) {
		NativeAssembler assem = new NativeAssembler();
		assem.setTruncateOutputOnRepeat(true);
		assem.setMaxContigs(50000);
		assem.setMaxPathsFromRoot(100000);
		assem.setKmer(43);
		assem.setReadLength(76);
		assem.setMinKmerFrequency(2);
		
		String bam1 = args[0];
		String bam2 = args[1];
		List<String> inputFiles = new ArrayList<String>();
		inputFiles.add(bam1);
		inputFiles.add(bam2);
		String output = 
		
		assem.assembleContigs(inputFiles, output, tempDir, region, prefix, checkForDupes, realigner)
		
//		assem.assembleContigs(args[0], args[1], "contig");
		
//		for (int i=0; i<10; i++) {
//			run(args[0], args[1] + "_" + i);
//		}
		
//		run(args[0], args[1]);
		
//		assem.assembleContigs("/home/lmose/code/abra/src/main/c/1810_reads.txt",
//				"/home/lmose/code/abra/src/main/c/1810.fa", "bar");
	}
*/
}

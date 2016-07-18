/* Copyright 2013 University of North Carolina at Chapel Hill.  All rights reserved. */
package abra;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.ValidationStringency;

public class CompareToReference {
	
	private StringBuffer reference;
	private String refFileName;
	private String currSeqName = "";
	private String cachedRefLine = null;
	private BufferedReader refReader;

	public void compare(String sam, String refFileName, int maxDiff) throws IOException, FileNotFoundException {
		refReader = new BufferedReader(new FileReader(refFileName));
		
		SAMFileReader reader = new SAMFileReader(new File(sam));
		reader.setValidationStringency(ValidationStringency.SILENT);
		
		int count = 0;
		int totalMapped = 0;
		
		for (SAMRecord read : reader) {
			if (!read.getReadUnmappedFlag()) {
				String seq = read.getReferenceName();
				if (!seq.equals(currSeqName)) {
					loadSeqRef(seq);
				}
				
//				String refStr = reference.substring(read.getAlignmentStart()-1, read.getAlignmentStart()+99);
//				if (!refStr.equals(read.getReadString())) {
//				StringBuffer diffStr = new StringBuffer();
				if (numDifferences(read) > maxDiff) {
					System.err.println("------------");
					System.err.println("read: " + read.getSAMString());
//					System.out.println("ref: " + refStr);
//					System.out.println("dif: " + diffStr.toString());
					count += 1;
				}
				
				totalMapped += 1;
			}
		}
		
		System.out.println("count: " + count + " out of: " + totalMapped);
		
		reader.close();
	}
	
	public void init(String reference) throws FileNotFoundException {
		refReader = new BufferedReader(new FileReader(reference));
	}
	
	public void cleanup() throws IOException {
		refReader.close();
	}
	
	public int numMismatches(SAMRecord read) throws IOException {
		int mismatches = 0;
		
		if (!read.getReadUnmappedFlag()) {
			String seq = read.getReferenceName();
			if (!seq.equals(currSeqName)) {
				loadSeqRef(seq);
			}
			
			mismatches = numDifferences(read);
		}

		return mismatches;
	}
	
	private int numDifferences(SAMRecord read) {
		int diffs = 0;
		int readIdx = 0;
		int refIdx = read.getAlignmentStart()-1;
		for (CigarElement element : read.getCigar().getCigarElements()) {
			if (element.getOperator() == CigarOperator.M) {
				for (int i=0; i<element.getLength(); i++) {
					char readBase = Character.toUpperCase(read.getReadString().charAt(readIdx));
					char refBase  = Character.toUpperCase(reference.charAt(refIdx));
					if ((readBase != refBase) && (readBase != 'N') && (refBase != 'N')) {
						diffs++;
					}
					
					readIdx++;
					refIdx++;
				}
			} else if (element.getOperator() == CigarOperator.I) {
				readIdx += element.getLength();
			} else if (element.getOperator() == CigarOperator.D) {
				refIdx += element.getLength();
			} else if (element.getOperator() == CigarOperator.S) {
				readIdx += element.getLength();
			}
		}
//		int diffs = 0;
//		int i = 0;
//		while ((i < s1.length()) && (i < s2.length())) {
//			int ch1 = Character.toUpperCase(s1.charAt(i));
//			int ch2 = Character.toUpperCase(s2.charAt(i));
//			if (ch1 != ch2 && ch1 != 'N' && ch2 != 'N') {
//				diffs += 1;
//				desc.append('*');
//			} else {
//				desc.append(' ');
//			}
//			
//			i +=1;
//		}
		
		return diffs;
	}
	
	private void loadSeqRef(String seq) throws IOException {
		System.out.println("Loading: " + seq);
		
		this.currSeqName = seq;
		String line = getRefLine();
		
		while ((line != null) && (!line.equals(">" + seq))) {
			line = getRefLine();
		}
		
		if (line != null) {
			line = getRefLine();
		}
		
		reference = new StringBuffer();
		while ((line != null) && (!line.startsWith(">"))) {
			reference.append(line);
			line = getRefLine();
			if ((line != null) && (line.startsWith(">"))) {
				cachedRefLine = line;
			}
		}
		
		System.out.println("Loaded: " + seq);
	}
	
	private String getRefLine() throws IOException {
		String line = null;
		if (cachedRefLine != null) {
			line = cachedRefLine;
			cachedRefLine = null;
		} else {
			line = refReader.readLine();
		}
		
		return line;
	}
	
	public static void main(String[] args) throws Exception {
		
		String ref = args[0];
		String sam = args[1];
		int max = Integer.parseInt(args[2]);
		
//		String ref = "/home/lmose/reference/chr16/chr16.fa";
//		String sam = "/home/lmose/dev/ayc/sim/sim80/sorted_rr4.bam";
//		String sam = "/home/lmose/dev/ayc/sim/sim80/sorted_rr.bam";
		//String sam = "/home/lmose/dev/ayc/sim/sim80/chr16.bam";
//		String sam = "/home/lmose/dev/ayc/sim/sim80/rchr16.bam";
		CompareToReference c2r = new CompareToReference();
		c2r.compare(sam, ref, max);
	}
}

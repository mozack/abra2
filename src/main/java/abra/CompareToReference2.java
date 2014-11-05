/* Copyright 2013 University of North Carolina at Chapel Hill.  All rights reserved. */
package abra;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMRecord;

/**
 * Utility class used to compare sequence to genomic reference.
 * 
 * @author Lisle E. Mose (lmose at unc dot edu)
 */
public class CompareToReference2 {
	
	private Random random = new Random(1);
	private String refFileName;
	private BufferedReader refReader;
	private Map<String, byte[]> refMap;
	private boolean is2Bit = true;

	/**
	 * Memory efficient reference storage using 2 bits per base.
	 * N's are converted to random bases.
	 */
	public void init(String reference) throws FileNotFoundException, IOException {
		this.refFileName = reference;
		loadRefMap();
	}
	
	/**
	 * Reference storage using the exact content of the reference file.  i.e. N's are stored explicitly.
	 * Non 2-bit storage is not fully tested!
	 */
	public void init8bit(String reference) throws FileNotFoundException, IOException {
		is2Bit = false;
		this.refFileName = reference;
		loadRefMap();
	}
	
	public boolean containsChromosome(String chromosome) {
		return refMap.containsKey(chromosome);
	}
	
	public List<String> getChromosomes() {
		return new ArrayList<String>(refMap.keySet());
	}
	
	public void cleanup() throws IOException {
		refReader.close();
	}
	
	public int numMismatches(SAMRecord read) {
		int mismatches = 0;
		
		if (!read.getReadUnmappedFlag()) {
			
			mismatches = numDifferences(read, 0);
		}

		return mismatches;
	}
	
	/**
	 * Returns # of read mismatches with bases exceeding minBaseQual + indel lengths.
	 * Soft clipped bases are included in comparison to reference. 
	 */
	public int numHighQualityMismatches(SAMRecord read, int minBaseQual) {
		int mismatches = 0;
		
		if (!read.getReadUnmappedFlag()) {
			
			mismatches = numDifferences(read, minBaseQual);
		}

		return mismatches;		
	}
	
	public List<Integer> mismatchPositions(SAMRecord read) {
		return mismatchPositions(read, -1);
	}
	
	private long getRefLength(String refName) {
		return refMap.get(refName.trim()).length * 4;
	}
	
	public String getAlternateReference(SAMRecord read, Cigar cigar) {
		String alt = null;
		
		if (read.getAlignmentEnd() < getRefLength(read.getReferenceName())) {
			
			StringBuffer altBuf = new StringBuffer(read.getReadLength());
		
			int readIdx = 0;
			int refIdx = read.getAlignmentStart()-1;
			for (CigarElement element : cigar.getCigarElements()) {
				if (element.getOperator() == CigarOperator.M) {
					
					for (int i=0; i<element.getLength(); i++) {

						if (refIdx >= getRefLength(read.getReferenceName())) {
							// You're off the edge of the map matey.  Monsters be here!
							// This read has aligned across chromosomes.  Do not proceed.
							return null;
						}
						
						char refBase = getRefBase(refIdx, read.getReferenceName());
						
						altBuf.append(refBase);
											
						readIdx++;
						refIdx++;
					}
				} else if (element.getOperator() == CigarOperator.I) {
					altBuf.append(read.getReadString().substring(readIdx, readIdx + element.getLength()));
					readIdx += element.getLength();
				} else if (element.getOperator() == CigarOperator.D) {
					refIdx += element.getLength();
				} else if (element.getOperator() == CigarOperator.S) {
					readIdx += element.getLength();
				}
			}
			
			alt = altBuf.toString();
		}
		
		return alt;
	}
	
	public List<Integer> mismatchPositions(SAMRecord read, int maxMismatches) {
		if (read.getReadUnmappedFlag()) {
			return Collections.emptyList();
		}
		
		List<Integer> mismatches = new ArrayList<Integer>();
		
		int readIdx = 0;
		int refIdx = read.getAlignmentStart()-1;
		for (CigarElement element : read.getCigar().getCigarElements()) {
			if (element.getOperator() == CigarOperator.M) {
				for (int i=0; i<element.getLength(); i++) {
					char readBase = getReadBase(read, readIdx);
					char refBase = getRefBase(refIdx, read.getReferenceName());
					if ((readBase != refBase) && (readBase != 'N') && (refBase != 'N')) {
						mismatches.add(readIdx);
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
			
			if ((maxMismatches > 0) && (mismatches.size() > maxMismatches)) {
				break;
			}
		}

		return mismatches;
	}
	
	private char getReadBase(SAMRecord read, int index) {
		return (char) read.getReadBases()[index];
	}
	
	private int getBaseQuality(SAMRecord read, int index) {
		return (char) read.getBaseQualities()[index];
	}
	
	private int numDifferences(SAMRecord read, int minBaseQual) {
		
		int diffs = 0;
		if (refMap.get(read.getReferenceName().trim()) != null) {
			int readIdx = 0;
			int refIdx = read.getAlignmentStart()-1;
			int elementIdx = 0;
			for (CigarElement element : read.getCigar().getCigarElements()) {
				if (element.getOperator() == CigarOperator.M) {
					for (int i=0; i<element.getLength(); i++) {
						char readBase = getReadBase(read, readIdx);
						char refBase = getRefBase(refIdx, read.getReferenceName());
						if ((readBase != refBase) && (readBase != 'N') && (refBase != 'N')) {
							if (minBaseQual == 0 || getBaseQuality(read, readIdx) >= minBaseQual) {
								diffs++;
							}
						}
						
						readIdx++;
						refIdx++;
					}
				} else if (element.getOperator() == CigarOperator.I) {
					readIdx += element.getLength();
				} else if (element.getOperator() == CigarOperator.D) {
					refIdx += element.getLength();
				} else if (element.getOperator() == CigarOperator.S) {
					
					if (elementIdx == 0) {
						refIdx -= element.getLength();
					}
					
					//TODO: Should this always be included?
					for (int i=0; i<element.getLength(); i++) {
						if ((refIdx >= 0) && (refIdx < getRefLength(read.getReferenceName())-1)) {
//							char readBase = Character.toUpperCase(read.getReadString().charAt(readIdx));
							char readBase = getReadBase(read, readIdx);
							char refBase = getRefBase(refIdx, read.getReferenceName());
							if ((readBase != refBase) && (readBase != 'N') && (refBase != 'N')) {
								if (minBaseQual == 0 || getBaseQuality(read, readIdx) >= minBaseQual) {
									diffs++;
								}
							}
						} else {
							if (minBaseQual > 0) {
								diffs++;
							}
						}
						
						readIdx++;
						refIdx++;
					}
				}
				
				elementIdx++;
			}
		}
		
		return diffs;
	}
	
	// Convert input StringBuffer to 2 bit representation.
	private byte[] getBytes(StringBuffer buf) {
		if (is2Bit) {
			int numBytes = buf.length()/4;
			if (buf.length() % 4 > 0) {
				numBytes += 1;
			}
			byte[] bytes = new byte[numBytes];
			
			int subIdx = 0;
			int byteIdx = 0;
			for (int i=0; i<buf.length(); i++) {
				
				byte base = getBase(buf.charAt(i));
				byte shifted = (byte) (base << (6-subIdx*2));
				
				bytes[byteIdx] = (byte) (bytes[byteIdx] | shifted);
				
				subIdx++;
				if (subIdx == 4) {
					subIdx = 0;
					byteIdx++;
				}
			}
			return bytes;
		} else {
			byte[] bytes = new byte[buf.length()];
			for (int i=0; i<buf.length(); i++) {
				bytes[i] = (byte) Character.toUpperCase(buf.charAt(i));
			}
			return bytes;
		}
	}
	
	private byte getBase(char ch) {
		switch (Character.toUpperCase(ch)) {
		case 'A':
			return 0;
		case 'T':
			return 1;
		case 'C':
			return 2;
		case 'G':
			return 3;
		default:
			return randomBase();
		}
	}
	
	private byte randomBase() {
		byte val;
		double rand = random.nextDouble();
		
		if (rand < .25) {
			val = 0;
		} else if (rand < .5) {
			val = 1;
		} else if (rand < .75) {
			val = 2;
		} else {
			val = 3;
		}
		
		return val;
	}
	
	private void loadRefMap() throws IOException {
		System.err.println("Loading reference map:  " + this.refFileName);
		long s = System.currentTimeMillis();
		this.refMap = new HashMap<String, byte[]>();
		
		BufferedReader reader = new BufferedReader(new FileReader(refFileName));
		
		String line = reader.readLine();
		StringBuffer sequence = new StringBuffer();
		String currSeqName = null;
		while (line != null) {
			if (line.startsWith(">")) {
				if (currSeqName != null) {
					System.err.println("\tChromosome: " + currSeqName + " length: " + sequence.length());
					refMap.put(currSeqName, getBytes(sequence));
				}
				
				sequence = new StringBuffer();
				
				currSeqName = line.substring(1, line.length()).trim();
				int spaceIdx = currSeqName.indexOf(' ');
				if (spaceIdx > 0) {
					currSeqName = currSeqName.substring(0, spaceIdx);
				}
				int tabIdx = currSeqName.indexOf('\t');
				if (tabIdx > 0) {
					currSeqName = currSeqName.substring(0, tabIdx);
				}
				
			} else {
				line.toUpperCase();
				sequence.append(line);
			}
			
			line = reader.readLine();
		}
		
		System.err.println("\tChromosome: " + currSeqName + " length: " + sequence.length());
		refMap.put(currSeqName, getBytes(sequence));
		
		sequence = null;
		
		reader.close();
		
		long e = System.currentTimeMillis();
		System.err.println("Done loading ref map.  Elapsed secs: " + (e-s)/1000);
	}
	
	private char getRefBase(int idx, String ref) {
		return getBaseAsChar(idx, refMap.get(ref.trim()));
	}
	
	private char getBaseAsChar(int idx, byte[] ref) {
		int byteIdx = idx / 4;
		int bitShift = (3-(idx % 4)) * 2;
		byte b = ref[byteIdx];
		byte shifted = (byte) (b >>> bitShift);
		byte val = (byte) (shifted & 3);
		
		switch (val) {
			case 0:
				return 'A';
			case 1:
				return 'T';
			case 2:
				return 'C';
			case 3:
				return 'G';
			default:
				throw new IllegalArgumentException("Invalid base value");
		}
	}
	
	public String getSequence(String chromosome, int position, int length) {
		byte[] ref = refMap.get(chromosome);
		
		if (ref == null) {
			System.err.println("No ref for chromosome: " + chromosome);
		}
		
		position -= 1;
		
		if (is2Bit) {
			StringBuffer buf = new StringBuffer(length);
			int start = Math.max(position, 0);
			int stop = Math.min(position+length, ref.length*4+1);
			
			for (int i=start; i<stop; i++) {
				buf.append(getBaseAsChar(i, ref));
			}
			return buf.toString();
		} else {
			byte[] sub = Arrays.copyOfRange(ref, Math.max(position,0), Math.min(position+length, ref.length));
			
			return new String(sub);
		}
	}
	
	/**
	 * Returns length of reference for input chromosome (give or take a few bases)
	 */
	public int getReferenceLength(String chromosome) {
		byte[] ref = refMap.get(chromosome);
		if (is2Bit) {
			return ref.length * 4;
		} else {
			return ref.length;
		}
	}
	
	/*
	public static void main(String[] args) {
		String foo = "ATCGNatcgn";
		System.out.println(foo);
		byte[] bytes = foo.getBytes();
		System.out.println("num bytes: " + bytes.length);
		for (int i=0; i<bytes.length; i++) {
			System.out.println("b: " + bytes[i] + "-" + ((char) bytes[i]));
		}
		
	}
	*/
	
	/*
	public static void main(String[] args) throws Exception {
		CompareToReference2 c2r = new CompareToReference2();
		
		StringBuffer sb = new StringBuffer();
		sb.append("ATCGGCT");
		
		byte[] bytes = c2r.getBytes(sb);
		
		System.out.println("len: " + bytes.length);
		for (byte b : bytes) {
			System.out.println(b);
		}
	}
	*/
	
	public static void main(String[] args) throws Exception {
		CompareToReference2 c2r = new CompareToReference2(); 
		c2r.init("/home/lmose/reference/chr3/chr3.fa");
//		c2r.init("/home/lmose/reference/test/test.fa");
		System.out.println(c2r.getSequence("chr3", 178948013, 100));

//		c2r.init("/home/lmose/reference/t2/t2.fa");
//		System.out.println(c2r.getSequence("x1", 1, 12));
	}
	
	/*
	public static void main(String[] args) throws Exception {
		
//		String ref = args[0];
//		String sam = args[1];
		
//		String ref = "/home/lmose/reference/chr7/chr7.fa";
//		String sam = "/home/lmose/dev/ayc/24/26/test";
		
		String ref = "/home/lmose/reference/chr1/1.fa";
		String sam = "/home/lmose/dev/abra_wxs/4/ttest1.bam";

		
		CompareToReference2 c2r = new CompareToReference2();
		c2r.init(ref);
		
//		Thread.sleep(100000);
		
		SAMFileReader rdr = new SAMFileReader(new File(sam));
		
		for (SAMRecord read : rdr) {
			int mismatches = c2r.numMismatches(read);
			
			System.out.println("mismatches: " + mismatches);
		}
		
		rdr.close();
	}
	*/
}

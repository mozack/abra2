package abra.utils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.ValidationStringency;

public class CompareToReference3 {
	
	private String refFileName;
	private String currSeqName = "";
	private String cachedRefLine = null;
	private BufferedReader refReader;
	private Map<String, byte[]> refMap;
	private int minBaseQuality = Integer.MIN_VALUE;
	
/*
	public void compare(String sam, String refFileName, int maxDiff) throws IOException, FileNotFoundException {
		loadRefMap();
		
		SAMFileReader reader = new SAMFileReader(new File(sam));
		reader.setValidationStringency(ValidationStringency.SILENT);
		
		int count = 0;
		int totalMapped = 0;
		
		for (SAMRecord read : reader) {
			if (!read.getReadUnmappedFlag()) {
				
				if (numDifferences(read) > maxDiff) {
					System.out.println("------------");
					System.out.println("read: " + read.getSAMString());
					count += 1;
				}
				
				totalMapped += 1;
			}
		}
		
		System.out.println("count: " + count + " out of: " + totalMapped);
		
		reader.close();
	}
	*/
	
	public void init(String reference) throws FileNotFoundException, IOException {
		this.refFileName = reference;
		loadRefMap();
	}
	
	public void cleanup() throws IOException {
		refReader.close();
	}
		
	public List<Integer> mismatchPositions(SAMRecord read) {
		return mismatchPositions(read, -1);
	}
	
	public String getAlternateReference(SAMRecord read, Cigar cigar) {
		String alt = null;
		
		byte[] reference = refMap.get(read.getReferenceName().trim());
		
		if (read.getAlignmentEnd() < reference.length) {
			
			StringBuffer altBuf = new StringBuffer(read.getReadLength());
		
			int readIdx = 0;
			int refIdx = read.getAlignmentStart()-1;
			for (CigarElement element : cigar.getCigarElements()) {
				if (element.getOperator() == CigarOperator.M) {
					
					for (int i=0; i<element.getLength(); i++) {

						if (refIdx >= reference.length) {
							// You're off the edge of the map matey.  Monsters be here!
							// This read has aligned across chromosomes.  Do not proceed.
							return null;
						}
						
						char refBase  = Character.toUpperCase((char) reference[refIdx]);
						
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
		
		byte[] reference = refMap.get(read.getReferenceName().trim());
		
		int readIdx = 0;
		int refIdx = read.getAlignmentStart()-1;
		for (CigarElement element : read.getCigar().getCigarElements()) {
			if (element.getOperator() == CigarOperator.M) {
				for (int i=0; i<element.getLength(); i++) {
					//char readBase = Character.toUpperCase(read.getReadString().charAt(readIdx));
					char readBase = getReadBase(read, readIdx);
					//char readBase = (char) read.getReadBases()[readIdx];
					char refBase  = Character.toUpperCase((char) reference[refIdx]);
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
		//return Character.toUpperCase(read.getReadString().charAt(index));
		return (char) read.getReadBases()[index];
	}
	
	private int getReadBaseQuality(SAMRecord read, int index) {
		return read.getBaseQualities()[index];
	}
	
	public int noiseAroundIndel(SAMRecord read, CigarOperator indelType, int indelPos, int indelLength) {
		
		int diffs = 0;
		byte[] reference = refMap.get(read.getReferenceName().trim());
		if (reference != null) {
			int readIdx = 0;
			int refIdx = read.getAlignmentStart()-1;
			int elementIdx = 0;
			for (CigarElement element : read.getCigar().getCigarElements()) {
				if (element.getOperator() == CigarOperator.M) {
					for (int i=0; i<element.getLength(); i++) {
						int readBaseQuality = getReadBaseQuality(read, readIdx);
						
						if (readBaseQuality >= minBaseQuality) {
						
							char readBase = getReadBase(read, readIdx);
							char refBase  = Character.toUpperCase((char) reference[refIdx]);
							if ((readBase != refBase) && (readBase != 'N') && (refBase != 'N')) {
								diffs++;
							}
						}
						
						readIdx++;
						refIdx++;
					}
				} else if (element.getOperator() == CigarOperator.I) {
					if ((indelType != CigarOperator.I) || (indelPos != refIdx) || (element.getLength() != indelLength)) {
						diffs++;
					}
					readIdx += element.getLength();
				} else if (element.getOperator() == CigarOperator.D) {
					if ((indelType != CigarOperator.D) || (indelPos != refIdx) || (element.getLength() != indelLength)) {
						diffs++;
					}
					refIdx += element.getLength();
				} else if (element.getOperator() == CigarOperator.S) {
//					readIdx += element.getLength();
					
					if (elementIdx == 0) {
						refIdx -= element.getLength();
					}
					
					//TODO: Should this always be included?
					for (int i=0; i<element.getLength(); i++) {
						int readBaseQuality = getReadBaseQuality(read, readIdx);
						
						if (readBaseQuality >= minBaseQuality) {
						
							if ((refIdx >= 0) && (refIdx < reference.length-1)) {
	//							char readBase = Character.toUpperCase(read.getReadString().charAt(readIdx));
								
								char readBase = getReadBase(read, readIdx);
								char refBase  = Character.toUpperCase((char) reference[refIdx]);
								if ((readBase != refBase) && (readBase != 'N') && (refBase != 'N')) {
									diffs++;
								}
							} else {
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
	
	private byte[] getBytes(StringBuffer buf) {
		byte[] bytes = new byte[buf.length()];
		for (int i=0; i<buf.length(); i++) {
			bytes[i] = (byte) buf.charAt(i);
		}
		return bytes;
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
		
		refMap.put(currSeqName, getBytes(sequence));
		
		sequence = null;
		
		reader.close();
		
		long e = System.currentTimeMillis();
		System.err.println("Done loading ref map.  Elapsed secs: " + (e-s)/1000);
	}
	
	public void setMinBaseQuality(int minBaseQuality) {
		this.minBaseQuality = minBaseQuality;
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
		
//		String ref = args[0];
//		String sam = args[1];
		
//		String ref = "/home/lmose/reference/chr7/chr7.fa";
//		String sam = "/home/lmose/dev/ayc/24/26/test";
		
		String ref = "/home/lmose/reference/chr1/1.fa";
		String sam = "/home/lmose/dev/abra_wxs/4/ttest1.bam";

		
		CompareToReference3 c2r = new CompareToReference3();
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

package abra;

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

import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMFileReader.ValidationStringency;

public class CompareToReference2 {
	
	private String refFileName;
	private String currSeqName = "";
	private String cachedRefLine = null;
	private BufferedReader refReader;
	private Map<String, StringBuffer> refMap;
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
	
	public int numMismatches(SAMRecord read) {
		int mismatches = 0;
		
		if (!read.getReadUnmappedFlag()) {
			
			mismatches = numDifferences(read);
		}

		return mismatches;
	}
	
	public List<Integer> mismatchPositions(SAMRecord read) {
		if (read.getReadUnmappedFlag()) {
			return Collections.emptyList();
		}
		
		List<Integer> mismatches = new ArrayList<Integer>();
		
		StringBuffer reference = refMap.get(read.getReferenceName().trim());
		
		int readIdx = 0;
		int refIdx = read.getAlignmentStart()-1;
		for (CigarElement element : read.getCigar().getCigarElements()) {
			if (element.getOperator() == CigarOperator.M) {
				for (int i=0; i<element.getLength(); i++) {
					char readBase = Character.toUpperCase(read.getReadString().charAt(readIdx));
					char refBase  = Character.toUpperCase(reference.charAt(refIdx));
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
		}

		return mismatches;
	}
	
	private int numDifferences(SAMRecord read) {
		
		int diffs = 0;
		StringBuffer reference = refMap.get(read.getReferenceName().trim());
		if (reference != null) {
			int readIdx = 0;
			int refIdx = read.getAlignmentStart()-1;
			int elementIdx = 0;
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
//					readIdx += element.getLength();
					
					if (elementIdx == 0) {
						refIdx -= element.getLength();
					}
					
					//TODO: Should this always be included?
					for (int i=0; i<element.getLength(); i++) {
						if ((refIdx >= 0) && (refIdx < reference.length()-1)) {
							char readBase = Character.toUpperCase(read.getReadString().charAt(readIdx));
							char refBase  = Character.toUpperCase(reference.charAt(refIdx));
							if ((readBase != refBase) && (readBase != 'N') && (refBase != 'N')) {
								diffs++;
							}
						} else {
							diffs++;
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
	
	private void loadRefMap() throws IOException {
		this.refMap = new HashMap<String, StringBuffer>();
		
		BufferedReader reader = new BufferedReader(new FileReader(refFileName));
		
		String line = reader.readLine();
		StringBuffer sequence = new StringBuffer();
		String currSeqName = null;
		while (line != null) {
			if (line.startsWith(">")) {
				if (currSeqName != null) {
					refMap.put(currSeqName, sequence);
				}
				
				currSeqName = line.substring(1, line.length()).trim();
				sequence = new StringBuffer();
			} else {
				sequence.append(line);
			}
			
			line = reader.readLine();
		}
		
		refMap.put(currSeqName, sequence);
		
		reader.close();
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
	public static void main(String[] args) throws Exception {
		
		String ref = args[0];
		String sam = args[1];
		int max = Integer.parseInt(args[2]);
		
//		String ref = "/home/lmose/reference/chr16/chr16.fa";
//		String sam = "/home/lmose/dev/ayc/sim/sim80/sorted_rr4.bam";
//		String sam = "/home/lmose/dev/ayc/sim/sim80/sorted_rr.bam";
		//String sam = "/home/lmose/dev/ayc/sim/sim80/chr16.bam";
//		String sam = "/home/lmose/dev/ayc/sim/sim80/rchr16.bam";
		CompareToReference2 c2r = new CompareToReference2();
		c2r.compare(sam, ref, max);
	}
	*/
}

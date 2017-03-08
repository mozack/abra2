/* Copyright 2013 University of North Carolina at Chapel Hill.  All rights reserved. */
package abra;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import abra.Logger.Level;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.TextCigarCodec;

/**
 * Utility class for shifting Indels into leftmost position.
 * 
 * @author Lisle E. Mose (lmose at unc dot edu)
 */
public class IndelShifter {
	
	public Cigar shiftAllIndelsLeft(int refStart, int refEnd, String chromosome, Cigar cigar, String seq, CompareToReference2 c2r) {
		
		try {
		
			List<CigarElement> elems = new ArrayList<CigarElement>(cigar.getCigarElements());
			
			int elemSize = elems.size();
			
			for (int i=1; i<elemSize-1; i++) {
				
				int refOffset = 0;
				int readOffset = 0;
				
				for (int j=0; j<i-1; j++) {
					refOffset += getRefOffset(elems.get(j));
					readOffset += getReadOffset(elems.get(j));
				}
				
				CigarElement prev = elems.get(i-1);
				CigarElement elem = elems.get(i);
				CigarElement next = elems.get(i+1);
				
				if ((elem.getOperator() == CigarOperator.DELETION || elem.getOperator() == CigarOperator.INSERTION) &&
					prev.getOperator() == CigarOperator.MATCH_OR_MISMATCH && next.getOperator() == CigarOperator.MATCH_OR_MISMATCH) {
					
					// subset cigar here and attempt to shift
					Cigar subCigar = new Cigar(Arrays.asList(prev, elem, next));
					String subSeq = seq.substring(readOffset, readOffset+subCigar.getReadLength());
					Cigar newCigar = shiftIndelsLeft(refStart + refOffset, refStart + refOffset + subCigar.getReferenceLength(),
							chromosome, subCigar, subSeq, c2r);
					
					//TODO: Merge abutting indels if applicable
					elems.set(i-1, newCigar.getCigarElement(0));
					elems.set(i, newCigar.getCigarElement(1));
					if (newCigar.getCigarElements().size() == 3) {
						elems.set(i+1, newCigar.getCigarElement(2));
					} else {
						elems.remove(i+1);
						elemSize -= 1;
					}
				}
			}
			
			mergeAbuttingElements(elems);
			return new Cigar(elems);
		} catch (RuntimeException e) {
			e.printStackTrace();
			Logger.error("Error on: " + refStart + ", " + refEnd + "," + chromosome + ", " + cigar + ", " + seq);
			throw e;
		}
	}
	
	// Merge neighboring Cigar elements of same type
	// mutates input
	private void mergeAbuttingElements(List<CigarElement> elems) {
		int origSize = elems.size();
		List<CigarElement> orig = Logger.LEVEL == Level.TRACE ? new ArrayList<CigarElement>(elems) : null;
		
		for (int i=1; i<elems.size(); i++) {
			if (elems.get(i).getOperator() == elems.get(i-1).getOperator()) {
				int length = elems.get(i).getLength() + elems.get(i-1).getLength();
				elems.set(i-1, new CigarElement(length, elems.get(i-1).getOperator()));
				elems.remove(i);
				i = i-1;
			}
		}
		
		if (Logger.LEVEL == Level.TRACE && origSize != elems.size()) {
			Cigar old = new Cigar(orig);
			Cigar curr = new Cigar(elems);
			Logger.trace("Cigar abutting elements merged: [%s] -> [%s]", old, curr);
		}
	}
	
	private int getRefOffset(CigarElement elem) {
		int offset = -1;
		switch(elem.getOperator()) {
			case M:
			case D:
				offset = elem.getLength();
				break;
			case I:
				offset = 0;
				break;
			default:
				throw new IllegalArgumentException("Unexpected Cigar operator: " + elem.getOperator());
		}
		
		return offset;
	}

	private int getReadOffset(CigarElement elem) {
		int offset = -1;
		switch(elem.getOperator()) {
			case M:
			case I:
				offset = elem.getLength();
				break;
			case D:
				offset = 0;
				break;
			default:
				throw new IllegalArgumentException("Unexpected Cigar operator: " + elem.getOperator());
		}
		
		return offset;
	}
	
	public Cigar shiftIndelsLeft(int refStart, int refEnd, String chromosome, Cigar cigar, String seq, CompareToReference2 c2r) {
		try {
			
			if (containsIndel(cigar)) {
				int indelPos = firstIndelOffset(cigar);
				String origReadAltRef = c2r.getAlternateReference(refStart, refEnd, chromosome, seq, cigar);
				
//				System.out.println("o: " + origReadAltRef);
				
				if (origReadAltRef != null) {
					for (int i=indelPos; i>0; i--) {
						Cigar newCigar = shiftCigarLeft(cigar, i);
						
						String shiftedReadAltRef = c2r.getAlternateReference(refStart, refEnd, chromosome, seq, newCigar);
						
						if ((shiftedReadAltRef != null) && (origReadAltRef.equals(shiftedReadAltRef))) {
							return newCigar;
						}
					}
				}
			}
		} catch (RuntimeException e) {
			Logger.error("Error processing: " + seq + ", " + cigar.toString());
			throw e;
		}
		
		return cigar;
	}
	
	private SAMRecord cloneRead(SAMRecord read) {
		try {
			return (SAMRecord) read.clone();
		} catch (CloneNotSupportedException e) {
			// Infamous "this should never happen" comment here.
			e.printStackTrace();
			throw new RuntimeException(e);
		}
	}
	
	protected Cigar shiftCigarLeft(Cigar cigar, int positionsToShift) {
		Cigar newCigar = new Cigar();
		
		for (int i=0; i<cigar.getCigarElements().size(); i++) {
			CigarElement elem = cigar.getCigarElement(i);
			
			if (isFirstNonSoftClippedElem(i, cigar)) {
				int newLen = elem.getLength() - positionsToShift;
				
				if (newLen > 0) {
					CigarElement newElem = new CigarElement(newLen, elem.getOperator());
					newCigar.add(newElem);
				}
			} else if (isLastNonSoftClippedElem(i, cigar)) {
				if (elem.getOperator() == CigarOperator.M) {
					CigarElement newElem = new CigarElement(elem.getLength() + positionsToShift, CigarOperator.M);
					newCigar.add(newElem);
				} else {
					CigarElement newElem = new CigarElement(positionsToShift, CigarOperator.M);
					newCigar.add(elem);
					newCigar.add(newElem);
				}
			} else {
				newCigar.add(elem);
			}
		}
		
		return newCigar;
	}
	
	private boolean isFirstNonSoftClippedElem(int idx, Cigar cigar) {
		// First element, not soft clipped
		if ((idx == 0) && (cigar.getCigarElement(0).getOperator() != CigarOperator.S)) {
			return true;
		}
		
		// Second element, with first element soft clipped
		if ((idx == 1) && (cigar.getCigarElement(0).getOperator() == CigarOperator.S)) {
			return true;
		}
		
		return false;
	}
	
	private boolean isLastNonSoftClippedElem(int idx, Cigar cigar) {
		int numElems = cigar.getCigarElements().size();
		
		// Last element, not soft clipped.
		if ((idx == numElems-1) && (cigar.getCigarElement(idx).getOperator() != CigarOperator.S)) {
			return true;
		}
		
		// Second to last element, with last element soft clipped.
		if ((idx == numElems-2) && (cigar.getCigarElement(numElems-1).getOperator() == CigarOperator.S)) {
			return true;
		}

		return false;
	}
	
	private boolean containsIndel(Cigar cigar) {
		for (CigarElement elem : cigar.getCigarElements()) {
			if ((elem.getOperator() == CigarOperator.D) || (elem.getOperator() == CigarOperator.I)) {
				return true;
			}
		}
		
		return false;
	}
	
	private int firstIndelOffset(Cigar cigar) {
		int pos = 0;
		
		for (CigarElement elem : cigar.getCigarElements()) {
			if ((elem.getOperator() == CigarOperator.D) || (elem.getOperator() == CigarOperator.I)) {
				return pos;
			} else if (elem.getOperator() != CigarOperator.S) {
				pos += elem.getLength();
			}
		}
		
		throw new IllegalArgumentException("No indel for record: [" + cigar + "]");
	}
	
	/*
	public static void main(String[] args) throws IOException {
		String in = args[0];
		String out = args[1];
		String ref = args[2];

//		String in = "/home/lmose/dev/abra/leftalign.sam";
//		String out = "/home/lmose/dev/abra/la.out.sam";
//		String ref = "/home/lmose/reference/chr1/chr1.fa";
		
		SamReader reader = SAMRecordUtils.getSamReader(in);
		
		SAMFileWriter writer = new SAMFileWriterFactory().makeSAMOrBAMWriter(
				reader.getFileHeader(), false, new File(out));
		
		CompareToReference2 c2r = new CompareToReference2();
		c2r.init(ref);
		
		IndelShifter indelShifter = new IndelShifter();

		for (SAMRecord read : reader) {
//			SAMRecord shiftedRead = indelShifter.shiftIndelsLeft(read, c2r);
//			writer.addAlignment(shiftedRead);
		}
		
		writer.close();
		reader.close();
	}
	*/
	
	public static void main(String[] args) throws IOException {
		CompareToReference2 c2r = new CompareToReference2();
		c2r.init("/home/lmose/dev/reference/hg19/chr7.fa");
		
		String seq = "CGCTGGCGGCCTTCCTGCCCGACCTGTCCCTGGCCAAGAGGAAGACCTGGCGGAACCCTTGGCGGCGTTCCTACCACAAGCGCAAGAAGGCCGAGTGGGAGGTGGACCCTGACAACTGCGAGGAGGTGAAGCAGACACCGCCCTACGACAGCAGCCACCGCATCCTGGACGTCATGGACATGACGATCTTCGACTTCCTCATGGGAAACATGGACCGTCACCACTACGAGACTTGTGAGAAGTTTGGGAATGAAACGTTCATCATCCACTTAGACAATGGAAGAGGGATCCGGA";
		IndelShifter is = new IndelShifter();
		Cigar orig = TextCigarCodec.decode("7M88I1M10D109M241D82M60D7M");
		
		Cigar cigar = is.shiftAllIndelsLeft(296604, 297121, "chr7", orig, 
				seq, 
				c2r);
		
		System.out.println("old: " + orig + ", " + orig.getReadLength());
		System.out.println("new: " + cigar + ", " + cigar.getReadLength());
		//296604, 297121,chr7, 7M88I1M10D109M241D82M60D7M, CGCTGGCGGCCTTCCTGCCCGACCTGTCCCTGGCCAAGAGGAAGACCTGGCGGAACCCTTGGCGGCGTTCCTACCACAAGCGCAAGAAGGCCGAGTGGGAGGTGGACCCTGACAACTGCGAGGAGGTGAAGCAGACACCGCCCTACGACAGCAGCCACCGCATCCTGGACGTCATGGACATGACGATCTTCGACTTCCTCATGGGAAACATGGACCGTCACCACTACGAGACTTGTGAGAAGTTTGGGAATGAAACGTTCATCATCCACTTAGACAATGGAAGAGGGATCCGGA
	}
}

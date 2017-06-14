/* Copyright 2013 University of North Carolina at Chapel Hill.  All rights reserved. */
package abra;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.DefaultSAMRecordFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.samtools.ValidationStringency;

/**
 * Utility methods for dealing with SAMRecord
 * 
 * @author Lisle E. Mose (lmose at unc dot edu)
 */
public class SAMRecordUtils {

	public static int getNumHighQualBases(SAMRecord read, int minBq) {
		int count = 0;
		for (byte bq : read.getBaseQualities()) {
			if (bq >= minBq) {
				count += 1;
			}
		}
		return count;
	}
	
	/**
	 * Replace hard clips with soft clips.
	 */
	public static void replaceHardClips(SAMRecord read) {
		Cigar cigar = read.getCigar();
		
		if (cigar.getCigarElements().size() > 0) {
			CigarElement firstElement = cigar.getCigarElement(0);
			CigarElement lastElement  = cigar.getCigarElement(cigar.numCigarElements()-1);
			
			if ((firstElement.getOperator() == CigarOperator.H) ||
				(lastElement.getOperator() == CigarOperator.H)) {
				
				Cigar newCigar = new Cigar();
				
				boolean isFirst = true;
				
				for (CigarElement element : cigar.getCigarElements()) {
					if (element.getOperator() != CigarOperator.H) {
						newCigar.add(element);
					} else {
						CigarElement softClip = new CigarElement(element.getLength(), CigarOperator.S);
						newCigar.add(softClip);
						
						if (isFirst) {
							read.setReadString(padBases(element.getLength()) + read.getReadString());
						} else {
							read.setReadString(read.getReadString() + padBases(element.getLength()));							
						}
					}
					
					isFirst = false;
				}
				
				read.setCigar(newCigar);
			}
		}
	}
	
	private static String padBases(int length) {
		StringBuffer buf = new StringBuffer(length);
		for (int i=0; i<length; i++) {
			buf.append('N');
		}
		return buf.toString();
	}
	
	/**
	 * Remove leading or trailing soft clips from the input read.
	 * Does not modify a read entirely comprised of soft clips.
	 */
	public static void removeSoftClips(SAMRecord read) {
		
		Cigar cigar = read.getCigar();
		
		CigarElement firstElement = cigar.getCigarElement(0);
		CigarElement lastElement  = cigar.getCigarElement(cigar.numCigarElements()-1);
		
		if ((firstElement.getOperator() == CigarOperator.S) ||
			(lastElement.getOperator() == CigarOperator.S)) {
		
			Cigar newCigar = new Cigar();
			
			String bases = read.getReadString();
			//String qualities = read.getBaseQualityString();
					
			if (firstElement.getOperator() == CigarOperator.S) {
				bases = bases.substring(firstElement.getLength(), bases.length());
				//qualities = qualities.substring(firstElement.getLength(), qualities.length()-1);
			} else {
				newCigar.add(firstElement);
			}
			
			for (int i=1; i<cigar.numCigarElements()-1; i++) {
				newCigar.add(cigar.getCigarElement(i));
			}
			
			if (lastElement.getOperator() == CigarOperator.S) {
				bases = bases.substring(0, bases.length() - lastElement.getLength());
				//qualities = qualities.substring(0, qualities.length() - lastElement.getLength() - 1);
			} else {
				newCigar.add(lastElement);
			}
			
			read.setCigar(newCigar);
			read.setReadString(bases);
			//read.setBaseQualityString(qualities);
		}
	}

	public static Cigar subset(Cigar cigar, int startIdx, int endIdx) {
		
		List<CigarElement> subElements = new ArrayList<CigarElement>();
		
		// Find first element and index into first element
		int len = 0;
		int elemIdx = -1;
		for (CigarElement elem : cigar.getCigarElements()) {

			// Treat deletions as zero length.
			int elemLength = elem.getOperator() == CigarOperator.D ? 0 : elem.getLength();
			
			if (elemIdx < 0) {
				
				// Find first element (Should never be a deletion)
				int elemStart = len;
				int elemStop  = len + elemLength;
				
				if ((startIdx >= elemStart) && (startIdx < elemStop)) {
					elemIdx = startIdx - elemStart;
					int elemLen = Math.min(elem.getLength()-elemIdx, endIdx-startIdx+1);
					CigarElement newElem = new CigarElement(elemLen, elem.getOperator());
					subElements.add(newElem);
				}
				
				len += elemLength;
				
			} else if ((len + elemLength) <= endIdx) {
				// Add this entire element
				subElements.add(elem);
				len += elemLength;
			} else if (len <= endIdx) {
				// Add part of last sub element (should never be a deletion)
				CigarElement newElem = new CigarElement(endIdx-len+1, elem.getOperator());
				subElements.add(newElem);
				break;
			} else {
				break;
			}
		}
		
		return new Cigar(subElements);
	}
	
	/**
	 * Calculates edit distance for the input read.
	 * If the input c2r is not null, compare to the actual reference.
	 * If c2r is null, check the input NM tag.
	 */
	public static int getEditDistance(SAMRecord read, CompareToReference2 c2r) {
		
		Integer distance = null;
		
		if (read.getReadUnmappedFlag()) {
			distance = read.getReadLength();
		} else if (c2r != null) {
			//distance = c2r.numMismatches(read) + getNumIndelBases(read);
			distance = c2r.numMismatches(read) + getNumIndelBases(read);
		} else {
			distance = read.getIntegerAttribute("NM");
			
			if (distance == null) {
				// STAR format
				distance = read.getIntegerAttribute("nM");
			}
			
			if (distance == null) {
				distance = read.getReadLength();
			}
		}
				
		return distance;
	}
	
	public static int getNumSplices(SAMRecord read) {
		int splices = 0;
		
		for (CigarElement element : read.getCigar().getCigarElements()) {
			if (element.getOperator() == CigarOperator.N) {
				splices += 1;
			}
		}
		
		return splices;
	}
		
	/**
	 *  Returns total length of deletions and insertions for the input read. 
	 */
	public static int getNumIndelBases(SAMRecord read) {
		int numIndelBases = 0;
		
		for (CigarElement element : read.getCigar().getCigarElements()) {
			if ((element.getOperator() == CigarOperator.D) || (element.getOperator() == CigarOperator.I)) {
				numIndelBases += element.getLength();
			}
		}
		
		return numIndelBases;
	}
	
	/**
	 *  Returns original edit distance as set in YX tag.
	 */
	public static int getOrigEditDistance(SAMRecord read) {
		
		Integer distance = null;
		
		if (read.getReadUnmappedFlag()) {
			distance = read.getReadLength();
		} else {
			distance = read.getIntegerAttribute("YX");
			
			if (distance == null) {
				distance = read.getReadLength();
			}
		}
				
		return distance;
	}
	
	/**
	 * Convenience method for retrieving int attribute
	 */
	public static int getIntAttribute(SAMRecord read, String attribute) {
		Integer num = read.getIntegerAttribute(attribute);
		
		if (num == null) {
			return 0;
		} else {
			return num;
		}
	}
	
	/**
	 * Return the number of insertions and deletions in a SAMRecord 
	 */
	public static int getNumIndels(SAMRecord read) {
		int numIndels = 0;
		
		for (CigarElement element : read.getCigar().getCigarElements()) {
			if ((element.getOperator() == CigarOperator.D) || (element.getOperator() == CigarOperator.I)) {
				numIndels += 1;
			}
		}
		
		return numIndels;
	}
	
	private static boolean isClip(CigarElement elem) {
		return elem.getOperator() == CigarOperator.S || elem.getOperator() == CigarOperator.H;
	}
	
	public static String getMappedReadPortion(SAMRecord read) {
		int start = 0; 
		
		List<CigarElement> elems = read.getCigar().getCigarElements();
		int i = 0;
		while (i < elems.size() && isClip(elems.get(i))) {
			if (elems.get(i).getOperator() == CigarOperator.S) {
				start += elems.get(i).getLength();
			}
			
			i += 1;
		}
		
		int stop = read.getReadLength();
		i = elems.size()-1;
		while (i >= 0 && isClip(elems.get(i))) {
			if (elems.get(i).getOperator() == CigarOperator.S) {
				stop -= elems.get(i).getLength();
			}
			
			i -= 1;
		}
		
		return read.getReadString().substring(start, stop);
	}
	
	public static String getLeadingClips(SAMRecord read) {
		List<CigarElement> elems = read.getCigar().getCigarElements();
		
		List<CigarElement> leading = new ArrayList<CigarElement>();
		
		for (CigarElement elem : elems) {
			if (isClip(elem)) {
				leading.add(elem);
			} else {
				break;
			}
		}
		
		String ret = "";
		if (leading.size() > 0) {
			Cigar cigar = new Cigar(leading);
			ret = TextCigarCodec.encode(cigar);
		}
		
		return ret;
	}
	
	public static String getTrailingClips(SAMRecord read) {
		List<CigarElement> elems = read.getCigar().getCigarElements();
		
		List<CigarElement> trailing = new ArrayList<CigarElement>();
		boolean isNonClippedReached = false;
		
		for (CigarElement elem : elems) {
			if (isClip(elem)) {
				if (isNonClippedReached) {
					trailing.add(elem);
				}
			} else {
				isNonClippedReached = true;
			}
		}

		String ret = "";
		if (trailing.size() > 0) {
			Cigar cigar = new Cigar(trailing);
			ret =  TextCigarCodec.encode(cigar);
		}
		
		return ret;
	}
	
	/**
	 *  Returns true if the updatedRead is essentially the same as the origRead
	 *  minus soft clipping. 
	 */
	public static boolean isSoftClipEquivalent(SAMRecord origRead, SAMRecord updatedRead) {
		
		boolean isEquivalent = false;
		
		if ((origRead.getCigarString().contains("S")) &&
			(origRead.getReferenceName().equals(updatedRead.getReferenceName())) &&
			(origRead.getReadNegativeStrandFlag() == updatedRead.getReadNegativeStrandFlag()) &&
			(origRead.getCigarLength() > 1)) {
			
			// Compare start positions
			int nonClippedOrigStart = origRead.getAlignmentStart();
			CigarElement firstElem = origRead.getCigar().getCigarElement(0); 
			if (firstElem.getOperator() == CigarOperator.S) {
				nonClippedOrigStart -= firstElem.getLength(); 
			}
			
			if (nonClippedOrigStart == updatedRead.getAlignmentStart()) {
				// Compare cigars
				List<CigarElement> elems = new ArrayList<CigarElement>(origRead.getCigar().getCigarElements());
				
				CigarElement first = elems.get(0);
				
				// If first element is soft clipped, lengthen the second element
				if (first.getOperator() == CigarOperator.S) {
					CigarElement second = elems.get(1);
					CigarElement newElem = new CigarElement(first.getLength() + second.getLength(), second.getOperator());
					elems.set(1,  newElem);
				}
				
				CigarElement last = elems.get(elems.size()-1);
				if (last.getOperator() == CigarOperator.S) {
					CigarElement nextToLast = elems.get(elems.size()-2);
					CigarElement newElem = new CigarElement(last.getLength() + nextToLast.getLength(), nextToLast.getOperator());
					elems.set(elems.size()-2, newElem);
				}
				
				List<CigarElement> newElems = new ArrayList<CigarElement>();

				for (CigarElement elem : elems) {
					if (elem.getOperator() != CigarOperator.S) {
						newElems.add(elem);
					}
				}
				
				Cigar convertedCigar = new Cigar(newElems);
				
				if (convertedCigar.equals(updatedRead.getCigar())) {
					isEquivalent = true;
				}
			}
		}
		
		return isEquivalent;
	}
	
	/**
	 * Returns a clone of the input read
	 */
	public static SAMRecord cloneRead(SAMRecord read) {
		try {
			return (SAMRecord) read.clone();
		} catch (CloneNotSupportedException e) {
			// Infamous "this should never happen" comment here.
			e.printStackTrace();
			throw new RuntimeException(e);
		}
	}

	/**
	 * Returns true if the input read should be filtered
	 */
	public static boolean isFiltered(boolean isPairedEnd, SAMRecord read) {
		// Filter out single end reads when in paired end mode.
		return ((isPairedEnd) && (!read.getReadPairedFlag()));
	}
	
	/**
	 * Returns true if the input read is primary.
	 * i.e. Bit flag not secondary 0x200 or supplemental 0x800
	 */
	public static boolean isPrimary(SAMRecord read) {
		return ((read.getFlags() & 0x800)  == 0) && (!read.getNotPrimaryAlignmentFlag());
	}
	
	public static int getMappedLength(SAMRecord read) {
		int length = read.getReadLength();
		for (CigarElement elem : read.getCigar().getCigarElements()) {
			if (elem.getOperator() == CigarOperator.S) {
				length -= elem.getLength();
			}
		}
		
		return length;
	}

	public static int sumBaseQuals(SAMRecord read) {
		int sum = 0;
		for (byte b : read.getBaseQualities()) {
			sum += b - '!';
		}
		return sum;
	}
	
	public static SamReader getSamReader(String filename) {

		return SamReaderFactory.make()
                        .validationStringency(ValidationStringency.SILENT)
                        .samRecordFactory(DefaultSAMRecordFactory.getInstance())
                        .open(new File(filename));
	}
	
	// From HTSJDK SAMUtils
	// Adapted to return a block for all non-clipped elems
	public static List<ReadBlock> getReadBlocks(final Cigar cigar, final int alignmentStart) {
        if (cigar == null) return Collections.emptyList();

        final List<ReadBlock> alignmentBlocks = new ArrayList<ReadBlock>();
        int readBase = 1;
        int refBase = alignmentStart;

        for (final CigarElement e : cigar.getCigarElements()) {
            switch (e.getOperator()) {
                case H:
                	alignmentBlocks.add(new ReadBlock(readBase, refBase, e.getLength()));
                    break; // ignore hard clips
                case P:
                	alignmentBlocks.add(new ReadBlock(readBase, refBase, e.getLength()));
                    break; // ignore pads
                case S:
                	alignmentBlocks.add(new ReadBlock(readBase, refBase, e.getLength()));
                    readBase += e.getLength();
                    break; // soft clip read bases
                case N:
                	alignmentBlocks.add(new ReadBlock(readBase, refBase, e.getLength()));
                    refBase += e.getLength();
                    break;  // reference skip
                case D:
                	alignmentBlocks.add(new ReadBlock(readBase, refBase, e.getLength()));
                    refBase += e.getLength();
                    break;
                case I:
                	alignmentBlocks.add(new ReadBlock(readBase, refBase, e.getLength()));
                    readBase += e.getLength();
                    break;
                case M:
                case EQ:
                case X:
                    final int length = e.getLength();
                    alignmentBlocks.add(new ReadBlock(readBase, refBase, length));
                    readBase += length;
                    refBase += length;
                    break;
                default:
                    throw new IllegalStateException("Case statement didn't deal with " + cigar.toString() + " op: " + e.getOperator());
            }
        }
        return Collections.unmodifiableList(alignmentBlocks);
    }
	
	public static class ReadBlock {
		private int readPos;
		private int refPos;
		private int length;
		
		public ReadBlock(int readPos, int refPos, int length) {
			this.readPos = readPos;
			this.refPos = refPos;
			this.length = length;
		}

		public int getReadPos() {
			return readPos;
		}

		public int getRefPos() {
			return refPos;
		}

		public int getLength() {
			return length;
		}
	}
}

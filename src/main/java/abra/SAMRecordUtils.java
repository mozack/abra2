/* Copyright 2013 University of North Carolina at Chapel Hill.  All rights reserved. */
package abra;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.lang.mutable.MutableFloat;

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
	
//	public static int getEditDistance(SAMRecord read) {
//		return getEditDistance(read, null);
//	}
	
	/**
	 * Calculates edit distance for the input read.
	 * If the input c2r is not null, compare to the actual reference.
	 * If c2r is null, check the input NM tag.
	 */
	public static int getEditDistance(SAMRecord read, CompareToReference2 c2r, boolean includeSoftClipping) {
		
		Integer distance = null;
		
		if (read.getReadUnmappedFlag()) {
			distance = read.getReadLength();
		} else if (c2r != null) {
			//distance = c2r.numMismatches(read) + getNumIndelBases(read);
			distance = c2r.numMismatches(read, includeSoftClipping) + getNumIndelBases(read);
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
	 * Return the number of insertions, deletions and introns in a SAMRecord 
	 */
	public static int getNumGaps(SAMRecord read) {
		int numGaps = 0;
		
		for (CigarElement element : read.getCigar().getCigarElements()) {
			if ((element.getOperator() == CigarOperator.D) || (element.getOperator() == CigarOperator.I) || 
				(element.getOperator() == CigarOperator.N)) {
				numGaps += 1;
			}
		}
		
		return numGaps;
	}
	
	public static int getNumIndels(SAMRecord read) {
		int numGaps = 0;
		
		for (CigarElement element : read.getCigar().getCigarElements()) {
			if ((element.getOperator() == CigarOperator.D) || (element.getOperator() == CigarOperator.I)) {
				numGaps += 1;
			}
		}
		
		return numGaps;
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
	
	private static Pair<String, String> consensusSeq(String s1, String s2, String qual1, String qual2, int maxMismatches,
			MutableFloat mismatchFrac) {
		StringBuffer consensus = new StringBuffer();
		StringBuffer consensusQual = new StringBuffer();
		int numMismatches = 0;
		
		for (int i=0; i<s1.length(); i++) {
			if (s1.charAt(i) != s2.charAt(i)) {

				numMismatches += 1;
				if (numMismatches > maxMismatches) {
					return null;
				} else {
					if (qual1.charAt(i) >= qual2.charAt(i) + 10) {
						consensus.append(s1.charAt(i));
						consensusQual.append(qual1.charAt(i));
					} else if (qual2.charAt(i) >= qual1.charAt(i) + 10) {
						consensus.append(s2.charAt(i));
						consensusQual.append(qual2.charAt(i));						
					} else {
						consensus.append('N');
						consensusQual.append('!');
					}
				}
				
			} else {
				consensus.append(s1.charAt(i));
				consensusQual.append((char) Math.max(qual1.charAt(i), qual2.charAt(i))); 
			}
		}
		
		mismatchFrac.setValue((float) numMismatches / (float) s1.length());
		
		return new Pair<String,String>(consensus.toString(), consensusQual.toString());
	}
	
	// TODO: Parameterize
	private static int MIN_OVERLAP = 10;
	
	public static Pair<String,String> mergeSequences(String seq1, String seq2, String qual1, String qual2) {
		
		if (seq1.length() < MIN_OVERLAP+10 || seq2.length() < MIN_OVERLAP+10) {
			return null;
		}
		
		// Identify first 1 kmers in seq2
		String[] head2Kmers = new String[1];
		for (int i=0; i<head2Kmers.length; i++) {
			head2Kmers[i] = seq2.substring(i, i+MIN_OVERLAP);
		}
		
		// Search for seq2 head kmers in seq1
		// Identify set of possible start positions in seq1
		Set<Integer> indices = new HashSet<Integer>();
		
		for (int i=0; i<head2Kmers.length; i++) {
			String head2Kmer = head2Kmers[i];
			int idx = seq1.indexOf(head2Kmer, i);
			while (idx > -1) {
				indices.add(idx - i);
				idx = seq1.indexOf(head2Kmer, idx+1);
			}
		}
		
		Pair<String,String> bestResult = null;
		float bestMismatchFrac = 2;
		
		for (int idx : indices) {
			
			String tail1 = seq1.substring(idx);
			int overlap = tail1.length();
			if (seq2.length() > overlap) {
				String head2 = seq2.substring(0, overlap);

				String qualTail1 = qual1.substring(idx);
				String qualHead2 = qual2.substring(0, overlap);
				
				MutableFloat mismatchFrac = new MutableFloat(0);
				
				// Require 90% of bases to match
				Pair<String,String> consensus = consensusSeq(tail1, head2, qualTail1, qualHead2, (int) (overlap * .10), mismatchFrac);

				if (consensus != null) {
					if (mismatchFrac.floatValue() < bestMismatchFrac) {
						bestMismatchFrac = mismatchFrac.floatValue();
						
						String mergedSeq  = seq1.substring(0, idx) + consensus.getFirst() + seq2.substring(overlap);
						String mergedQual = qual1.substring(0, idx) + consensus.getSecond() + qual2.substring(overlap);
						bestResult = new Pair<String,String>(mergedSeq,mergedQual);
					}
				}
			}
		}
		
		return bestResult;
	}
	
	public static int mergeReadPair(SAMRecordWrapper readWrapper, Map<String, SAMRecordWrapper> firstReads, Map<String, SAMRecordWrapper> secondReads) {
		
		int alignmentStart = -1;
		SAMRecord read = readWrapper.getSamRecord();
		
		if (read.getReadPairedFlag() && !read.getReadUnmappedFlag() && !read.getMateUnmappedFlag() &&
				read.getReadNegativeStrandFlag() != read.getMateNegativeStrandFlag()) {
			
			SAMRecordWrapper pair = null;
			
			if (read.getFirstOfPairFlag()) {
				pair = secondReads.get(read.getReadName() + "_" + read.getMateAlignmentStart());
			} else {
				pair = firstReads.get(read.getReadName() + "_" + read.getMateAlignmentStart());
			}
			
			if (pair != null) {
				SAMRecordWrapper first = null;
				SAMRecordWrapper second = null;
				if (read.getReadNegativeStrandFlag()) {
					first = pair;
					second = readWrapper;
				} else {
					first = readWrapper;
					second = pair;
				}

				
				if (first.getAdjustedAlignmentStart() < second.getAdjustedAlignmentStart() &&
						first.getAdjustedAlignmentEnd() > second.getAdjustedAlignmentStart() &&
						first.getSamRecord().getReadLength() > MIN_OVERLAP && 
						second.getSamRecord().getReadLength() > MIN_OVERLAP) {

					Pair<String, String> merged = mergeSequences(first.getSamRecord().getReadString(), second.getSamRecord().getReadString(),
							first.getSamRecord().getBaseQualityString(), second.getSamRecord().getBaseQualityString());
					
					if (merged != null) {
						readWrapper.setMerged(merged.getFirst(), merged.getSecond(), first.getAdjustedAlignmentStart(), second.getAdjustedAlignmentEnd());
						pair.setMerged(merged.getFirst(), merged.getSecond(), first.getAdjustedAlignmentStart(), second.getAdjustedAlignmentEnd());
						
						alignmentStart = first.getAdjustedAlignmentStart();
					}
				}
			}
		}
		
		return alignmentStart;
	}
	/*
	public static void mergeReadPairOld(SAMRecordWrapper readWrapper, Map<String, SAMRecordWrapper> firstReads, Map<String, SAMRecordWrapper> secondReads) {
		
		SAMRecord read = readWrapper.getSamRecord();
		
		if (read.getReadPairedFlag() && !read.getReadUnmappedFlag() && !read.getMateUnmappedFlag() &&
				read.getReadNegativeStrandFlag() != read.getMateNegativeStrandFlag()) {
			
			SAMRecordWrapper pair = null;
			
			if (read.getFirstOfPairFlag()) {
				pair = secondReads.get(read.getReadName() + "_" + read.getMateAlignmentStart());
			} else {
				pair = firstReads.get(read.getReadName() + "_" + read.getMateAlignmentStart());
			}
			
			if (pair != null) {
				SAMRecordWrapper first = null;
				SAMRecordWrapper second = null;
				if (read.getReadNegativeStrandFlag()) {
					first = pair;
					second = readWrapper;
				} else {
					first = readWrapper;
					second = pair;
				}
				
				int MIN_OVERLAP = 10;

				if (first.getAdjustedAlignmentStart() < second.getAdjustedAlignmentStart() &&
						first.getAdjustedAlignmentEnd() > second.getAdjustedAlignmentStart() &&
						first.getSamRecord().getReadLength() > MIN_OVERLAP && 
						second.getSamRecord().getReadLength() > MIN_OVERLAP) {
					
					SAMRecord firstRead = first.getSamRecord();
					SAMRecord secondRead = second.getSamRecord();
					String tail = firstRead.getReadString().substring(firstRead.getReadLength()-MIN_OVERLAP);
					int idx = secondRead.getReadString().indexOf(tail);
					if (idx >= 0) {
						int overlap = idx + MIN_OVERLAP; // Number of overlapping bases
						if (overlap < firstRead.getReadLength()) {
							tail = firstRead.getReadString().substring(firstRead.getReadLength()-overlap);
							String head = secondRead.getReadString().substring(0, overlap);
							
							// Requiring exact match here
							if (tail.equals(head)) {
								String mergeSeq = firstRead.getReadString() + secondRead.getReadString().substring(overlap);
								
								// TODO: Adjust base quality for merged sequence
								String mergedQual = firstRead.getBaseQualityString() + secondRead.getBaseQualityString().substring(overlap);
								readWrapper.setMergedSeqAndQual(mergeSeq, mergedQual);
								pair.setMergedSeqAndQual(mergeSeq, mergedQual);
								
								Logger.trace("OLD Merging: %s : %s", firstRead.getReadName(), mergeSeq);
							}
						}
					}
				}
			}
		}
	}
	*/
	
	public static boolean hasPossibleAdapterReadThrough(SAMRecord read, Map<String, SAMRecordWrapper> firstReads, Map<String, SAMRecordWrapper> secondReads) {
		
		boolean hasReadThrough = false;
		
		// Check for fragment read through in paired end
		if (read.getReadPairedFlag() && !read.getReadUnmappedFlag() && !read.getMateUnmappedFlag() &&
				read.getAlignmentStart() == read.getMateAlignmentStart() &&
				read.getReadNegativeStrandFlag() != read.getMateNegativeStrandFlag()) {
			
			SAMRecordWrapper pair = null;
			
			if (read.getFirstOfPairFlag()) {
				pair = secondReads.get(read.getReadName() + "_" + read.getMateAlignmentStart());
			} else {
				pair = firstReads.get(read.getReadName() + "_" + read.getMateAlignmentStart());
			}
			
			if (pair != null && read.getCigar().getCigarElements().size() > 0 && pair.getSamRecord().getCigar().getCigarElements().size() > 0) {
				
				// Looking for something like:
				//     -------->
				//  <--------
				SAMRecord first = null;
				SAMRecord second = null;
				if (read.getReadNegativeStrandFlag()) {
					first = read;
					second = pair.getSamRecord();
				} else {
					first = pair.getSamRecord();
					second = read;
				}
				
				CigarElement firstElement = first.getCigar().getFirstCigarElement();
				CigarElement lastElement = second.getCigar().getLastCigarElement();
				
				if (firstElement.getOperator() == CigarOperator.S && lastElement.getOperator() == CigarOperator.S) {
					// We likely have read through into adapter here.
					hasReadThrough = true;
				}
			}
		}

		return hasReadThrough;
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
	
	public static int sumBaseQuals(String quals) {
		int sum = 0;
		
		for (int i=0; i<quals.length(); i++) {
			char ch = (char) (quals.charAt(i) - '!');
			sum += ch;
		}
		
		return sum;
	}

	public static int sumBaseQuals(SAMRecord read) {
		int sum = 0;
		for (byte b : read.getBaseQualities()) {
			sum += b;
		}
		return sum;
	}
	
	public static int getUnclippedLength(SAMRecord read) {
		int softClipLen = 0;
		for (CigarElement elem : read.getCigar().getCigarElements()) {
			if (elem.getOperator() == CigarOperator.S) {
				softClipLen += elem.getLength();
			}
		}
		
		return read.getReadLength() - softClipLen;
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

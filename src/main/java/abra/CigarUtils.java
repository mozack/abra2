package abra;

import java.util.ArrayList;
import java.util.List;

import abra.ReadEvaluator.Alignment;


public class CigarUtils {
	
	/**
	 * Extract subset of cigar string based upon input position (index) into cigar and length. 
	 */
	public static int subsetCigarString(int pos, int length, String cigar, StringBuffer newCigar) {
		List<CigarBlock> cigarBlocks = getCigarBlocks(cigar);
		List<CigarBlock> newCigarBlocks = new ArrayList<CigarBlock>();
		int relativeRefPos = subsetCigarBlocks(cigarBlocks, pos, length, newCigarBlocks);
		
		for (CigarBlock block : newCigarBlocks) {
			newCigar.append(block.length);
			newCigar.append(block.type);
		}
		
		return relativeRefPos;
	}
	
	private static String cigarStringFromCigarBlocks(List<CigarBlock> blocks) {
		StringBuffer newCigar = new StringBuffer();
		
		for (CigarBlock block : blocks) {
			newCigar.append(block.length);
			newCigar.append(block.type);
		}
		
		return newCigar.toString();
	}
	
	public static String extendCigarWithMatches(String cigar, int leftPad, int rightPad) {
		List<CigarBlock> blocks = getCigarBlocks(cigar);
		
		if (blocks.get(0).type == 'M') {
			blocks.get(0).length += leftPad;
		} else {
			blocks.add(0, new CigarBlock(leftPad, 'M'));
		}
		
		int lastBlockIdx = blocks.size()-1;
		
		if (blocks.get(lastBlockIdx).type == 'M') {
			blocks.get(lastBlockIdx).length += rightPad;
		} else {
			blocks.add(new CigarBlock(rightPad, 'M'));
		}
		
		return cigarStringFromCigarBlocks(blocks);
	}
	
	public static String injectSplice(String cigar, int junctionPos, int junctionLength) {
		
		// Identify pos relative to reference and insert N element
		List<CigarBlock> blocks = getCigarBlocks(cigar);
		List<CigarBlock> newBlocks = new ArrayList<CigarBlock>();
		int refPos = 0;

		for (CigarBlock block : blocks) {
			if (block.type == 'M' || block.type == 'D') {
				if (refPos < junctionPos && refPos + block.length >= junctionPos) {
					// Split up current block into 2 blocks with splice block in between
					int blockLen1 = junctionPos - refPos;
					int blockLen2 = block.length - blockLen1;
					newBlocks.add(new CigarBlock(blockLen1, block.type));
					newBlocks.add(new CigarBlock(junctionLength, 'N'));
					if (blockLen2 > 0) {
						newBlocks.add(new CigarBlock(blockLen2, block.type));
					}
					
					refPos += block.length;
				} else {
					newBlocks.add(block);
					refPos += block.length;
				}
			} else {
				// Do not advance ref pos for insertions or introns
				newBlocks.add(block);
			}
		}
		
		return cigarStringFromCigarBlocks(newBlocks);
	}
	
	// Assumes input junctions are sorted by coordinate
	public static String injectSplices(String cigar, List<Integer> junctionPos, List<Integer> junctionLength) {
		
		for (int i=0; i<junctionPos.size(); i++) {
			cigar = injectSplice(cigar, junctionPos.get(i), junctionLength.get(i));
		}
		
		return cigar;
	}
	
	private static int selectPrimaryAlignment(Alignment alignment1, Alignment alignment2, int def) {
		int selection = def;

		if (!alignment1.isSecondary && alignment2.isSecondary) {
			selection = 1;
		} else if (alignment1.isSecondary && !alignment2.isSecondary) {
			selection = 2;
		}
		
		return selection;
	}
	
	// Return 0 if cigars are not equivalent (treating deletions the same as junctions)
	// Return 1 if cigar1 has more junctions
	// Return 2 if cigar2 has more junctions
	public static int testEquivalenceAndSelectIntronPreferred(Alignment alignment1, Alignment alignment2) {
		
		String cigar1 = alignment1.cigar;
		String cigar2 = alignment2.cigar;
		
		// Cigars are equal, pick non-secondary or just the first.
		if (cigar1.equals(cigar2)) {
			return (selectPrimaryAlignment(alignment1, alignment2, 1));
		}
		
		// Cigars are different, pick non-secondary one or neither
		if (cigar1.length() != cigar2.length()) {
			return (selectPrimaryAlignment(alignment1, alignment2, 0));
		}
		
		int cigar1Introns = 0;
		int cigar2Introns = 0;
		
		for (int i=0; i<cigar1.length(); i++) {
			char ch1 = cigar1.charAt(i);
			char ch2 = cigar2.charAt(i);
			if (ch1 != ch2) {
				if ((ch1 != 'N' && ch1 != 'D') ||
					(ch2 != 'N' && ch2 != 'D')) {
					// Non-equivalent cigars
					return (selectPrimaryAlignment(alignment1, alignment2, 0));
				} else {
					if (ch1 == 'N') {
						cigar1Introns += 1;
					}
					if (ch2 == 'N') {
						cigar2Introns += 1;
					}
				}
			}
		}
		
		if (cigar1Introns != cigar2Introns) {
			return cigar1Introns >= cigar2Introns ? 1 : 2;
		}
		
		// Equivalent cigars.  Pick non-secondary or just the first.
		return (selectPrimaryAlignment(alignment1, alignment2, 0));		
	}
	
	private static List<CigarBlock> getCigarBlocks(String cigar) {
		
		List<CigarBlock> cigarBlocks = new ArrayList<CigarBlock>();
		try {
			StringBuffer len = new StringBuffer();
			for (int i=0; i<cigar.length(); i++) {
				char ch = cigar.charAt(i);
				if (Character.isDigit(ch)) {
					len.append(ch);
				} else {
					cigarBlocks.add(new CigarBlock(Integer.valueOf(len.toString()), ch));
					len.setLength(0);
				}
			}
		} catch (NumberFormatException e) {
			Logger.error("NumberFormatException: " + cigar);
			throw e;
		}
		
		return cigarBlocks;
	}
	
	private static int subsetCigarBlocks(List<CigarBlock> contigBlocks, int pos, int readLength, List<CigarBlock> readBlocks) {
		int currLen = 0;
		int contigPos = 0;
		int relativeRefPos = 0;
		boolean isReadPosReached = false;
//		List<CigarBlock> readBlocks = new ArrayList<CigarBlock>();
		
		for (CigarBlock block : contigBlocks) {
			
			int blockLength = block.length;
			
			// Identify the start point for subsetting
			if (!isReadPosReached) {
				if (!block.isGap()) {  // Never start in a deletion
					if (contigPos + block.length >= pos) {
						blockLength = contigPos + block.length - pos;
						isReadPosReached = true;
						
						if (block.type != 'I') {
							// Include partial block length for matches
							relativeRefPos += block.length - blockLength;
						}
					} else {
						contigPos += block.length;
						if (block.type != 'I') {
							// Include entire block for matches
							relativeRefPos += block.length;
						}
					}
				} else {
					// Include entire block for deletes
					relativeRefPos += block.length;
				}				
			} 
			
			if (isReadPosReached && blockLength > 0) {
				if (block.isGap()) {
					// Never start in a deletion
					if (!readBlocks.isEmpty()) {
						readBlocks.add(block);
					} else {
						// skip over leading deletion in reference position
						relativeRefPos += block.length;
					}
				}
				else if (blockLength < readLength-currLen) {
					currLen += blockLength;
					readBlocks.add(new CigarBlock(blockLength, block.type));
				} else {
					int len = readLength - currLen;
					currLen += len;
					readBlocks.add(new CigarBlock(len, block.type));
					break;
				}
			}
		}
		
		return relativeRefPos;
	}
	
	static class CigarBlock {
		int length;
		char type;
		
		CigarBlock(int length, char type) {
			this.length = length;
			this.type = type;
		}
		
		boolean isGap() {
			return type == 'D' || type == 'N';
		}
	}

}

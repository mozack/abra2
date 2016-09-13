package abra;

import java.util.ArrayList;
import java.util.List;


public class CigarUtils {
	
	/**
	 * Extract subset of cigar string based upon input position (index) into cigar and length. 
	 */
	public static String subsetCigarString(int pos, int length, String cigar) {
		List<CigarBlock> cigarBlocks = getCigarBlocks(cigar);
		List<CigarBlock> newCigarBlocks = subsetCigarBlocks(cigarBlocks, pos, length);
		
		StringBuffer newCigar = new StringBuffer();
		for (CigarBlock block : newCigarBlocks) {
			newCigar.append(block.length);
			newCigar.append(block.type);
		}
		
		return newCigar.toString();
	}
	
	private static List<CigarBlock> getCigarBlocks(String cigar) {
		
		List<CigarBlock> cigarBlocks = new ArrayList<CigarBlock>();
		
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
		
		return cigarBlocks;
	}
	
	private static List<CigarBlock> subsetCigarBlocks(List<CigarBlock> contigBlocks, int pos, int readLength) {
		int currLen = 0;
		int contigPos = 0;
		boolean isReadPosReached = false;
		List<CigarBlock> readBlocks = new ArrayList<CigarBlock>();
		
		for (CigarBlock block : contigBlocks) {
			
			int blockLength = block.length;
			
			// Identify the start point for subsetting
			if (!isReadPosReached) {
				if (block.type != 'D') {  // Never start in a deletion
					if (contigPos + block.length >= pos) {
						blockLength = contigPos + block.length - pos;
						isReadPosReached = true;
					} else {
						contigPos += block.length;
					}
				}
			} 
			
			if (isReadPosReached && blockLength > 0) {
				if (block.type == 'D') {
					// Never start in a deletion
					if (!readBlocks.isEmpty()) {
						readBlocks.add(block);
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
		
		return readBlocks;
	}
	
	static class CigarBlock {
		int length;
		char type;
		
		CigarBlock(int length, char type) {
			this.length = length;
			this.type = type;
		}
	}

}

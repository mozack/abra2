package abra.cadabra;

import abra.ReverseComplementor;

public class HomopolymerRun {
	
	private int length;
	private char base;
	private int pos;

	/**
	 *  Search for short HP run on either side of position. 
	 */
	public static HomopolymerRun find(String seq) {
	
//		String seq = c2r.getSequence(chromosome, position-9, 20);
		
		int maxLen = 0;
		char maxBase = '0';
		int maxPos = -1;
		
		int length = 1;
		char prev = '0';
		int i = 0;
		
		while (i < seq.length()) {
			if (seq.charAt(i) == prev) {
				length += 1;
				
				if (length > maxLen) {
					maxLen = length;
					maxBase = seq.charAt(i);
					maxPos = i - (length-1);
				}
			} else {
				length = 1;
			}
			
			prev = seq.charAt(i);
			i += 1;
		}
		
		HomopolymerRun hrun = null;
		
		if (maxLen >= 4) {
			hrun = new HomopolymerRun(maxLen, maxBase, maxPos);
		}
		
		return hrun;
	}
	
	public HomopolymerRun(int length, char base, int pos) {
		this.length = length;
		this.base = base;
		this.pos = pos;
	}

	public int getLength() {
		return length;
	}

	public char getBase() {
		return base;
	}
	
	public int getPos() {
		return pos;
	}
}

package abra;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class SemiGlobalAligner {
	
	public enum Direction { UP, LEFT, DIAG, NONE };
	
	private int match = 8;
	private int mismatch = -32;
	private int gapOpen = -48;
	private int gapExtend = -1;
	
	// Bit flags for backtrack matrix
	private static final byte LEFT_OPEN = 1;
	private static final byte UP_OPEN = 2;
	private static final byte DIR_DIAG = 4;
	private static final byte DIR_UP = 8;
	private static final byte DIR_LEFT = 16; 
		
	public SemiGlobalAligner() {
	}
	
	public SemiGlobalAligner(int match, int mismatch, int gapOpen, int gapExtend) {
		this.match = match;
		this.mismatch = mismatch;
		this.gapOpen = gapOpen;
		this.gapExtend = gapExtend;
	}

	public Result align(String seq1, String seq2) {
		int[][] matrix = new int[seq1.length()+1][seq2.length()+2];
		byte[][] bt = new byte[seq1.length()+1][seq2.length()+2];
		
		populate(matrix, bt, seq1, seq2);
		return backtrack(matrix, bt, seq1, seq2);
	}
	
	private void populate(int[][] matrix, byte[][] bt, String seq1, String seq2) {
		int col0Init = 0;
		for (int r=0; r<=seq1.length(); r++) {
			matrix[r][0] = col0Init;
			col0Init += gapOpen;
		}
		
		for (int r=1; r<=seq1.length(); r++) {
			for (int c=1; c<=seq2.length(); c++) {
				
				int prevCell = matrix[r-1][c-1];
				int prevBt;
				int diagScore = seq1.charAt(r-1) == seq2.charAt(c-1) ? prevCell + match : prevCell + mismatch;

				prevCell = matrix[r][c-1];
				prevBt = bt[r][c-1];
				int leftScore = (prevBt & LEFT_OPEN) == LEFT_OPEN ? prevCell + gapExtend : prevCell + gapOpen;

				prevCell = matrix[r-1][c];
				prevBt = bt[r-1][c];
				int upScore   = (prevBt & UP_OPEN) == UP_OPEN ? prevCell + gapExtend : prevCell + gapOpen;
				
				int max = max(diagScore, leftScore, upScore);				
				
				byte currBt = 0;

				if (diagScore == max) {
					currBt = DIR_DIAG;
				}
				
				if (upScore == max) {				
					currBt = (byte) (DIR_UP | UP_OPEN);
				}
				
				if (leftScore == max) {				
					currBt = (byte) (DIR_LEFT | LEFT_OPEN);
				}
				
				matrix[r][c] = max;
				bt[r][c] = currBt;
			}
		}
	}
	
	private void dumpMatrix(Cell matrix[][]) {
		for(int r=0; r<matrix.length; r++) {
			StringBuffer line = new StringBuffer();
			for (int c=0; c<matrix[0].length; c++) {
				line.append(matrix[r][c].score);
				line.append(' ');
				switch(matrix[r][c].prev) {
					case DIAG:
						line.append('\\');
						break;
					case UP:
						line.append('^');
						break;
					case LEFT:
						line.append('<');
					default:
						break;
				}
				line.append('\t');
			}
			System.out.println(line);
		}
	}
	
	private Element updateCurrElem(char op, Element currElem, List<Element> elems) {
		if (currElem.operator == op) {
			currElem.length += 1;
		} else {
			currElem = new Element(op, 1);
			elems.add(currElem);
		}
		
		return currElem;
	}
	
	private Result backtrack(int matrix[][], byte[][] bt, String seq1, String seq2) {
//		dumpMatrix(matrix);
		// Find best score in last row (end of seq1)
		int bestIdx = -1;
		int bestScore = -30000;
		int secondBestScore = -30000;
		int row = seq1.length();
		
		for (int c=1; c<=seq2.length(); c++) {
			if (matrix[row][c] > bestScore) {
				bestIdx = c;
				bestScore = matrix[row][c];
			} else if (matrix[row][c] > secondBestScore) {
				secondBestScore = matrix[row][c]; 
			}
		}
		
		int r = seq1.length();
		int c = bestIdx;
		int refEndIdx = c;
		
		List<Element> elems = new ArrayList<Element>();
		Element currElem = new Element('0', 0);
		
		while (r > 0 && c > 0) {
			byte currBt = bt[r][c];

			if ((currBt & DIR_DIAG) == DIR_DIAG) {
				r -= 1;
				c -= 1;
				currElem = updateCurrElem('M', currElem, elems);
			} else if ((currBt & DIR_LEFT) == DIR_LEFT) {
				c -= 1;
				currElem = updateCurrElem('D', currElem, elems);
			} else if ((currBt & DIR_UP) == DIR_UP) {
				r -= 1;
				currElem = updateCurrElem('I', currElem, elems);
			} else {
				break;
			}
		}
		
		int refIdx = c;
		
		Collections.reverse(elems);
		StringBuffer cigar = new StringBuffer();
		for (Element elem : elems) {
			cigar.append(elem.length);
			cigar.append(elem.operator);
		}
		
		return new Result(bestScore, secondBestScore, refIdx, refEndIdx, cigar.toString());
	}
	
	int max(int s1, int s2, int s3) {
		return Math.max(Math.max(s1, s2), s3);
	}
	
	static class Result {
		Result(int score, int secondBest, int position, int endPosition, String cigar) {
			this.score = score;
			this.secondBest = secondBest;
			this.position = position;
			this.endPosition = endPosition;
			this.cigar = cigar;
		}
		
		int score;
		int secondBest;
		int position;
		int endPosition;
		String cigar;
		
		public String toString() {
			return String.format("score: %d, secondBest: %d, pos: %d, endPos: %d, cigar: %s", score, secondBest, position, endPosition, cigar);
		}
	}
	
	static class Cell {
		int score;
		Direction prev;
		boolean leftOpen;
		boolean upOpen;
		
		Cell(int score, Direction prev, boolean leftOpen, boolean upOpen) {
			this.score = score;
			this.prev = prev;
			this.leftOpen = leftOpen;
			this.upOpen = upOpen;
		}
	}
	
	static class Element {
		char operator;
		int length;
		
		Element(char operator, int length) {
			this.operator = operator;
			this.length = length;
		}
	}
	
	public static void main(String[] args) throws Exception {

//		double i = 1;
//		i += -1e8;
		
//		System.out.println(i);
	
		long start = System.currentTimeMillis();
		Result res = null;
		for (int i=0; i<1000; i++) {
			String ref = "CCAGATCAGCCTAGGCAACATGGTGAAACCCCGTCTCTACCAAAAATAAAAAACTTAGCTGAGCGTGGTGGTGCACGCCTGTAGCCCCAGCTGCTGAGGAGCCTGAGCCCAGGGGGTGGAGGCTGCAGTGAGCCATGATCACACTACTGTACTCCAGCCTAGGTGACAGAGTGAGACCCTGTCTCAAAAAAATAAAAGAAAATAAAAATAAACAAAGAGAGAAGTGGAAGAAGAGGTGGAGTTTTGTATTTATGACTTGAATTTTGTATTCATGACTGGGTTGACACCCCAATCCACTCCATTTTTAGCCTTGAAACATGGCAAACAGTAACCATTAAAAGGATGGAAAAGAGAAGAAGGCATGGGTGGGAAACTGTGCCTCCCATTTTTGTGCATCTTTGTTGCTGTCCTTCCACTATACTGTACCTTTCAGCATTTTGACGGCAACCTGGATTGAGACTCCTGTTTTGCTAATTCCATAAGCTGTTGCGTTCATCACTTTTCCAAAAGCACCTGATCCTAGTACCTTCCCTGCAAAGACAAATGGTGAGTACGTGCATTTTAAAGATTTTCCAATGGAAAAGAAATGCTGCAGAAACATTTGGCACATTCCATTCTTACCAAACTCTAAATTTTCTCTTGGAAACTCCCATTTGAGATCATATTCATATTCTCTGAAATCAACGTAGAAGTACTCATTATCTGAGGAGCCGGTCACCTGTACCATCTGTAGCTGGCTTTCATACCTAAATTGCTTCAGAGATGAAATGATGAGTCAGTTAGGAATAGGCAGTTCTGCAGATAGAGGAAAGAATAATGAATTTTTACCTTTGCTTTTACCTTTTTGTACTTGTGACAAATTAGCAGGGTTAAAACGACAATGAAGAGGAGACAAACACCAATTGTTGCATAGAATGAGATGTTGTCTTGGATGAAAGGGAAGGGGCCTGCAACAAAAGAGTGTCACTCAGCGATGAAACAGAATTCCTGTGTGACATTATAAATAGTGGACAACTCATTATAATCTCTCACATCCTGTTTCAGTAATAATCATTTTCAGTCCTAACAACCACTCTACATATACTCTACTCCCCACAGACAATCAGGCAATGTCCCTGTAAAGGATACATTTCCTCCCTAGAAAATTGCGGATTATTCTCAATCCATTCTTTAAAACCATTTACTAGGGTAAATTTACAAGAATTACATCTGGTCCAGGCACGATGGCTCACGCCTGTAGTCCCAGCACTTTGGGAGGCCAAGATGGGAGGATCACTTGAGTCCAAGAATTAGACACCAGCCCAGGCAACACAGTGAAATCCCGTCTCTAAAAAAATTCAAAAATTAGCTGGGCGTGGTGGCAGGTGCCTGTAATCCCAGCTGCTCGGGAGGCTGAGGCAGGAG";
			String seq = "CCATTCTTACCAAACTCTAAATTTTCTCTTGGAAACTCCCATTTGAGATCATATTCATATTCTCTGAAATCAACGTAGAAGTACTCATTATCTGAGGAGCCGCACATTCCATTCTTGCCAAACTCTAGATTTTCTCTTGGAAACTCCCATTTGAGATCACATTCATATTCTCTGAAATCAACGTAGAAGTACTCATTATCTGAGGAGCCGGTCACCCGTACCATCTGTAGC";
			
			SemiGlobalAligner sga = new SemiGlobalAligner(8, -32, -48, -1);
			res = sga.align(seq, ref);
		}
		long end = System.currentTimeMillis();
		System.out.println("Elapsed: " + (end-start));
		System.out.println(String.format("score: %d, second: %d, position: %d, cigar: %s", res.score, res.secondBest, res.position, res.cigar));
		
//		Thread.sleep(300*1000);
	}
}

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
	
	public SemiGlobalAligner() {
	}
	
	public SemiGlobalAligner(int match, int mismatch, int gapOpen, int gapExtend) {
		this.match = match;
		this.mismatch = mismatch;
		this.gapOpen = gapOpen;
		this.gapExtend = gapExtend;
	}

	public Result align(String seq1, String seq2) {
		Cell[][] matrix = new Cell[seq1.length()+1][seq2.length()+1];
		
		populate(matrix, seq1, seq2);
		return backtrack(matrix, seq1, seq2);
	}
	
	private void populate(Cell[][] matrix, String seq1, String seq2) {
		int col0Init = 0;
		for (int r=0; r<=seq1.length(); r++) {
			matrix[r][0] = new Cell(col0Init, Direction.NONE, false, false);
			col0Init += gapOpen;
		}
		
		for (int c=0; c<=seq2.length(); c++) {
			matrix[0][c] = new Cell(0, Direction.NONE, false, false);
		}
		
		for (int r=1; r<=seq1.length(); r++) {
			for (int c=1; c<=seq2.length(); c++) {
				
				Cell prevCell = matrix[r-1][c-1];
				int diagScore = seq1.charAt(r-1) == seq2.charAt(c-1) ? prevCell.score + match : prevCell.score + mismatch;

				prevCell = matrix[r][c-1];
				int leftScore = prevCell.leftOpen ? prevCell.score + gapExtend : prevCell.score + gapOpen;

				prevCell = matrix[r-1][c];
				int upScore   = prevCell.upOpen ? prevCell.score + gapExtend : prevCell.score + gapOpen;
				
				int max = max(diagScore, leftScore, upScore);				
				Cell cell = null;
				
				Direction dir = null;
				boolean leftOpen = false;
				boolean upOpen = false;

				if (diagScore == max) {
					dir = Direction.DIAG; 
				}
				
				if (upScore == max) {
					dir = Direction.UP;
					upOpen = true;
				}
				
				if (leftScore == max) {
					dir = Direction.LEFT;
					leftOpen = true;
				}
				
				cell = new Cell(max, dir, leftOpen, upOpen);
				
				matrix[r][c] = cell;
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
	
	private Result backtrack(Cell matrix[][], String seq1, String seq2) {
//		dumpMatrix(matrix);
		// Find best score in last row (end of seq1)
		int bestIdx = -1;
		int bestScore = -30000;
		int secondBestScore = -30000;
		int row = seq1.length();
		
		for (int c=1; c<=seq2.length(); c++) {
			if (matrix[row][c].score > bestScore) {
				bestIdx = c;
				bestScore = matrix[row][c].score;
			} else if (matrix[row][c].score > secondBestScore) {
				secondBestScore = matrix[row][c].score; 
			}
		}
		
		int r = seq1.length();
		int c = bestIdx;
		int refEndIdx = c;
		
		List<Element> elems = new ArrayList<Element>();
		Element currElem = new Element('0', 0);
		
		while (r > 0 && c > 0) {
			Cell cell = matrix[r][c]; 

			if (cell.prev == Direction.DIAG) {
				r -= 1;
				c -= 1;
				currElem = updateCurrElem('M', currElem, elems);
			} else if (cell.prev == Direction.LEFT) {
				c -= 1;
				currElem = updateCurrElem('D', currElem, elems);
			} else if (cell.prev == Direction.UP) {
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
	
	public static void main(String[] args) {

		Result res = null;
//		for (int i=0; i<1000; i++) {
			String ref = "CCAGATCAGCCTAGGCAACATGGTGAAACCCCGTCTCTACCAAAAATAAAAAACTTAGCTGAGCGTGGTGGTGCACGCCTGTAGCCCCAGCTGCTGAGGAGCCTGAGCCCAGGGGGTGGAGGCTGCAGTGAGCCATGATCACACTACTGTACTCCAGCCTAGGTGACAGAGTGAGACCCTGTCTCAAAAAAATAAAAGAAAATAAAAATAAACAAAGAGAGAAGTGGAAGAAGAGGTGGAGTTTTGTATTTATGACTTGAATTTTGTATTCATGACTGGGTTGACACCCCAATCCACTCCATTTTTAGCCTTGAAACATGGCAAACAGTAACCATTAAAAGGATGGAAAAGAGAAGAAGGCATGGGTGGGAAACTGTGCCTCCCATTTTTGTGCATCTTTGTTGCTGTCCTTCCACTATACTGTACCTTTCAGCATTTTGACGGCAACCTGGATTGAGACTCCTGTTTTGCTAATTCCATAAGCTGTTGCGTTCATCACTTTTCCAAAAGCACCTGATCCTAGTACCTTCCCTGCAAAGACAAATGGTGAGTACGTGCATTTTAAAGATTTTCCAATGGAAAAGAAATGCTGCAGAAACATTTGGCACATTCCATTCTTACCAAACTCTAAATTTTCTCTTGGAAACTCCCATTTGAGATCATATTCATATTCTCTGAAATCAACGTAGAAGTACTCATTATCTGAGGAGCCGGTCACCTGTACCATCTGTAGCTGGCTTTCATACCTAAATTGCTTCAGAGATGAAATGATGAGTCAGTTAGGAATAGGCAGTTCTGCAGATAGAGGAAAGAATAATGAATTTTTACCTTTGCTTTTACCTTTTTGTACTTGTGACAAATTAGCAGGGTTAAAACGACAATGAAGAGGAGACAAACACCAATTGTTGCATAGAATGAGATGTTGTCTTGGATGAAAGGGAAGGGGCCTGCAACAAAAGAGTGTCACTCAGCGATGAAACAGAATTCCTGTGTGACATTATAAATAGTGGACAACTCATTATAATCTCTCACATCCTGTTTCAGTAATAATCATTTTCAGTCCTAACAACCACTCTACATATACTCTACTCCCCACAGACAATCAGGCAATGTCCCTGTAAAGGATACATTTCCTCCCTAGAAAATTGCGGATTATTCTCAATCCATTCTTTAAAACCATTTACTAGGGTAAATTTACAAGAATTACATCTGGTCCAGGCACGATGGCTCACGCCTGTAGTCCCAGCACTTTGGGAGGCCAAGATGGGAGGATCACTTGAGTCCAAGAATTAGACACCAGCCCAGGCAACACAGTGAAATCCCGTCTCTAAAAAAATTCAAAAATTAGCTGGGCGTGGTGGCAGGTGCCTGTAATCCCAGCTGCTCGGGAGGCTGAGGCAGGAG";
			String seq = "CCATTCTTACCAAACTCTAAATTTTCTCTTGGAAACTCCCATTTGAGATCATATTCATATTCTCTGAAATCAACGTAGAAGTACTCATTATCTGAGGAGCCGCACATTCCATTCTTGCCAAACTCTAGATTTTCTCTTGGAAACTCCCATTTGAGATCACATTCATATTCTCTGAAATCAACGTAGAAGTACTCATTATCTGAGGAGCCGGTCACCCGTACCATCTGTAGC";
			
			SemiGlobalAligner sga = new SemiGlobalAligner(8, -32, -48, -1);
			res = sga.align(seq, ref);
//		}
		System.out.println(String.format("score: %d, second: %d, position: %d, cigar: %s", res.score, res.secondBest, res.position, res.cigar));

	}
}

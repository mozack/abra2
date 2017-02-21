package abra;

import java.util.ArrayList;
import java.util.List;

public class SemiGlobalAligner {
	
	enum Direction { UP, LEFT, DIAG, NONE };
	
	private short match = 8;
	private short mismatch = -32;
	private short gapOpen = -48;
	private short gapExtend = -1;

//	private short match = 1;
//	private short mismatch = -4;
//	private short gapOpen = -6;
//	private short gapExtend = 0;
	
	public SemiGlobalAligner() {
	}
	
	public SemiGlobalAligner(int match, int mismatch, int gapOpen, int gapExtend) {
		this.match = (short) match;
		this.mismatch = (short) mismatch;
		this.gapOpen = (short) gapOpen;
		this.gapExtend = (short) gapExtend;
	}

	public Result align(String seq1, String seq2) {
		Cell[][] matrix = new Cell[seq1.length()+1][seq2.length()+1];
		
		short col0Init = 0;
		for (int r=0; r<=seq1.length(); r++) {
			matrix[r][0] = new Cell((short) col0Init, Direction.NONE);
			col0Init += gapOpen;
		}
		
		for (int c=0; c<=seq2.length(); c++) {
			matrix[0][c] = new Cell((short) 0, Direction.NONE);
		}
		
		for (int r=1; r<=seq1.length(); r++) {
			for (int c=1; c<=seq2.length(); c++) {
				
				Cell prevCell = matrix[r-1][c-1];
				short diagScore = (short) (seq1.charAt(r-1) == seq2.charAt(c-1) ? prevCell.score + match : prevCell.score + mismatch);

				prevCell = matrix[r][c-1];
				short leftScore = (short) (prevCell.prev == Direction.LEFT ? prevCell.score + gapExtend : prevCell.score + gapOpen);

				prevCell = matrix[r-1][c];
				short upScore   = (short) (prevCell.prev == Direction.UP ? prevCell.score + gapExtend : prevCell.score + gapOpen);
				
				short max = max(diagScore, leftScore, upScore);				
				Cell cell = null;
				
				if (diagScore == max) {
					cell = new Cell(diagScore, Direction.DIAG); 
				} else if (leftScore == max) {
					cell = new Cell(leftScore, Direction.LEFT);
				} else if (upScore == max) {
					cell = new Cell(upScore, Direction.UP);
				} else {
					throw new IllegalStateException(String.format("Somehow no max score found!: %s,%s", seq1, seq2));
				}
				
				matrix[r][c] = cell;
			}
		}
		
		return backtrack(matrix, seq1, seq2);
	}
	
	private void dumpMatrix(Cell matrix[][]) {
		for(int r=0; r<matrix.length; r++) {
			StringBuffer line = new StringBuffer();
			for (int c=0; c<matrix[0].length; c++) {
				line.append(matrix[r][c].score);
				line.append('\t');
			}
			System.out.println(line);
		}
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
		
		List<Cell> path = new ArrayList<Cell>();
		
		while (r > 0 && c > 0) {
			Cell cell = matrix[r][c]; 
			path.add(cell);
			if (cell.prev == Direction.DIAG) {
				r -= 1;
				c -= 1;
			} else if (cell.prev == Direction.LEFT) {
				c -= 1;
			} else if (cell.prev == Direction.UP) {
				r -= 1;
			} else {
				break;
			}
		}
		
		int refIdx = c;
		
		StringBuffer cigar = new StringBuffer();
		Direction dir = Direction.NONE;
		int elemLength = 0;

		for (int i=path.size()-1; i>=0; i--) {
			Cell cell = path.get(i);
			if (cell.prev == dir) {
				elemLength += 1;
			} else {
				if (elemLength > 0) {
					cigar.append(elemLength);
					if (dir == Direction.DIAG) {
						cigar.append('M');
					} else if (dir == Direction.LEFT) {
						cigar.append('D');
					} else if (dir == Direction.UP) {
						cigar.append('I');
					} else {
						cigar.append('Z');
					}
				}
				
				dir = cell.prev;
				elemLength = 1;
			}
		}
		
		// Get last element
		if (elemLength > 0) {
			cigar.append(elemLength);
			if (dir == Direction.DIAG) {
				cigar.append('M');
			} else if (dir == Direction.LEFT) {
				cigar.append('D');
			} else if (dir == Direction.UP) {
				cigar.append('I');
			} else {
				cigar.append('Z');
			}
		}
		
		return new Result((short) bestScore, (short) secondBestScore, (short) refIdx, (short) refEndIdx, cigar.toString());
	}
	
	short max(short s1, short s2, short s3) {
		return (short) Math.max(Math.max(s1, s2), s3);
	}
	
	static class Result {
		Result(short score, short secondBest, short position, short endPosition, String cigar) {
			this.score = score;
			this.secondBest = secondBest;
			this.position = position;
			this.endPosition = endPosition;
			this.cigar = cigar;
		}
		
		short score;
		short secondBest;
		short position;
		short endPosition;
		String cigar;
		
		public String toString() {
			return String.format("score: %d, secondBest: %d, pos: %d, endPos: %d, cigar: %s", score, secondBest, position, endPosition, cigar);
		}
	}
	
	static class Cell {
		Cell(short score, Direction prev) {
			this.score = score;
			this.prev = prev;
		}
		short score;
		Direction prev;
	}
	
	public static void main(String[] args) {
		String ref = "CCAGATCAGCCTAGGCAACATGGTGAAACCCCGTCTCTACCAAAAATAAAAAACTTAGCTGAGCGTGGTGGTGCACGCCTGTAGCCCCAGCTGCTGAGGAGCCTGAGCCCAGGGGGTGGAGGCTGCAGTGAGCCATGATCACACTACTGTACTCCAGCCTAGGTGACAGAGTGAGACCCTGTCTCAAAAAAATAAAAGAAAATAAAAATAAACAAAGAGAGAAGTGGAAGAAGAGGTGGAGTTTTGTATTTATGACTTGAATTTTGTATTCATGACTGGGTTGACACCCCAATCCACTCCATTTTTAGCCTTGAAACATGGCAAACAGTAACCATTAAAAGGATGGAAAAGAGAAGAAGGCATGGGTGGGAAACTGTGCCTCCCATTTTTGTGCATCTTTGTTGCTGTCCTTCCACTATACTGTACCTTTCAGCATTTTGACGGCAACCTGGATTGAGACTCCTGTTTTGCTAATTCCATAAGCTGTTGCGTTCATCACTTTTCCAAAAGCACCTGATCCTAGTACCTTCCCTGCAAAGACAAATGGTGAGTACGTGCATTTTAAAGATTTTCCAATGGAAAAGAAATGCTGCAGAAACATTTGGCACATTCCATTCTTACCAAACTCTAAATTTTCTCTTGGAAACTCCCATTTGAGATCATATTCATATTCTCTGAAATCAACGTAGAAGTACTCATTATCTGAGGAGCCGGTCACCTGTACCATCTGTAGCTGGCTTTCATACCTAAATTGCTTCAGAGATGAAATGATGAGTCAGTTAGGAATAGGCAGTTCTGCAGATAGAGGAAAGAATAATGAATTTTTACCTTTGCTTTTACCTTTTTGTACTTGTGACAAATTAGCAGGGTTAAAACGACAATGAAGAGGAGACAAACACCAATTGTTGCATAGAATGAGATGTTGTCTTGGATGAAAGGGAAGGGGCCTGCAACAAAAGAGTGTCACTCAGCGATGAAACAGAATTCCTGTGTGACATTATAAATAGTGGACAACTCATTATAATCTCTCACATCCTGTTTCAGTAATAATCATTTTCAGTCCTAACAACCACTCTACATATACTCTACTCCCCACAGACAATCAGGCAATGTCCCTGTAAAGGATACATTTCCTCCCTAGAAAATTGCGGATTATTCTCAATCCATTCTTTAAAACCATTTACTAGGGTAAATTTACAAGAATTACATCTGGTCCAGGCACGATGGCTCACGCCTGTAGTCCCAGCACTTTGGGAGGCCAAGATGGGAGGATCACTTGAGTCCAAGAATTAGACACCAGCCCAGGCAACACAGTGAAATCCCGTCTCTAAAAAAATTCAAAAATTAGCTGGGCGTGGTGGCAGGTGCCTGTAATCCCAGCTGCTCGGGAGGCTGAGGCAGGAG";
		String seq = "CCATTCTTACCAAACTCTAAATTTTCTCTTGGAAACTCCCATTTGAGATCATATTCATATTCTCTGAAATCAACGTAGAAGTACTCATTATCTGAGGAGCCGCACATTCCATTCTTGCCAAACTCTAGATTTTCTCTTGGAAACTCCCATTTGAGATCACATTCATATTCTCTGAAATCAACGTAGAAGTACTCATTATCTGAGGAGCCGGTCACCCGTACCATCTGTAGC";
		
		SemiGlobalAligner sga = new SemiGlobalAligner(8, -32, -48, -1);
		Result res = sga.align(seq, ref);
		System.out.println(String.format("score: %d, second: %d, position: %d, cigar: %s", res.score, res.secondBest, res.position, res.cigar)); 
	}
}

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
//	private static final byte LEFT_OPEN = 1;
//	private static final byte UP_OPEN = 2;
	private static final byte DIR_UP = 1;
	private static final byte DIR_DIAG = 2;
	private static final byte DIR_LEFT = 3;
	
	private static final int I = 0;
	private static final int M = 1;
	private static final int D = 2;
	
	private static final int MAX_REF_LEN = 5000;
	private static final int MAX_CONTIG_LEN = 2000;

	private int[][][] matrix = new int[MAX_CONTIG_LEN][MAX_REF_LEN][3];
	private byte[][][] bt = new byte[MAX_CONTIG_LEN][MAX_REF_LEN][3];
	
	public SemiGlobalAligner() {

	}
	
	public SemiGlobalAligner(int match, int mismatch, int gapOpen, int gapExtend) {
		this.match = match;
		this.mismatch = mismatch;
		this.gapOpen = gapOpen;
		this.gapExtend = gapExtend;
	}

	public Result align(String seq1, String seq2) {
//		int[][][] matrix = new int[seq1.length()+1][seq2.length()+1][3];
//		byte[][][] bt = new byte[seq1.length()+1][seq2.length()+1][3];
		
		populate(matrix, bt, seq1, seq2);
		return backtrack(matrix, bt, seq1, seq2);
	}
	
	private void populate(int[][][] matrix, byte[][][] bt, String seq1, String seq2) {
		for (int r=1; r<=seq1.length(); r++) {
			matrix[r][0][I] = gapOpen + (r*gapExtend);
			matrix[r][0][M] = gapOpen + (r*gapExtend);
			matrix[r][0][D] = gapOpen + (r*gapExtend);
		}
		
		for (int c=0; c<=seq2.length(); c++) {
			matrix[0][c][I] = gapOpen + (c*gapExtend);
			matrix[0][c][M] = 0;
			matrix[0][c][D] = gapOpen + (c*gapExtend);
		}
		
		for (int r=1; r<=seq1.length(); r++) {
			for (int c=1; c<=seq2.length(); c++) {
				
				//
				// Insertion (lower) matrix
				int insertExt = matrix[r-1][c][I] + gapExtend;
				int insertOpen = matrix[r-1][c][M] + gapOpen;
				
				if (insertExt >= insertOpen) {
					matrix[r][c][I] = insertExt;
					bt[r][c][I] = DIR_UP;
				} else {
					matrix[r][c][I] = insertOpen;
					bt[r][c][I] = DIR_DIAG;
				}
				
				// Deletion (upper) matrix
				int deleteExt = matrix[r][c-1][D] + gapExtend;
				int deleteOpen = matrix[r][c-1][M] + gapOpen;

				if (deleteExt >= deleteOpen) {
					matrix[r][c][D] = deleteExt;
					bt[r][c][D] = DIR_LEFT;
				} else {
					matrix[r][c][D] = deleteOpen;
					bt[r][c][D] = DIR_DIAG;
				}
				
				//
				// Match/mismatch (middle) matrix
				int insertClose = matrix[r][c][I]; 
				int baseMatch = seq1.charAt(r-1) == seq2.charAt(c-1) ? matrix[r-1][c-1][M] + match : matrix[r-1][c-1][M] + mismatch;
				int deleteClose = matrix[r][c][D];
				
				if (baseMatch >=insertClose && baseMatch >= deleteClose) {
					matrix[r][c][M] = baseMatch;
					bt[r][c][M] = DIR_DIAG;
				} else if (insertClose >= deleteClose) {
					matrix[r][c][M] = insertClose;
					bt[r][c][M] = DIR_UP;					
				} else {
					matrix[r][c][M] = deleteClose;
					bt[r][c][M] = DIR_LEFT;										
				}
			}
		}
		
//		dumpMatrix(seq1, seq2, matrix, bt);
	}
	
	String intString(int[] i) {
		return String.format("%d,%d,%d", i[0],i[1],i[2]);
	}
	
	String byteString(byte[] i) {
		return String.format("%d,%d,%d", i[0],i[1],i[2]);
	}
	
	private void dumpMatrix(String seq1, String seq2, int matrix[][][], byte bt[][][]) {
		
		StringBuffer line = new StringBuffer();
		line.append("\t\t");
		for (int i=0; i<seq2.length(); i++) {
			line.append(seq2.charAt(i));
			line.append('\t');
		}
		System.out.println(line);
		
		for(int r=0; r<matrix.length; r++) {
			line = new StringBuffer();
			if (r > 0) {
				line.append(seq1.charAt(r-1));
			}
			line.append('\t');
			for (int c=0; c<matrix[0].length; c++) {
				line.append(intString(matrix[r][c]));
				line.append(':');
				line.append(byteString(bt[r][c]));
				
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
	
	private Result backtrack(int matrix[][][], byte[][][] bt, String seq1, String seq2) {
//		dumpMatrix(matrix);
		// Find best score in last row (end of seq1)
		int bestIdx = -1;
		int bestScore = -300000000;
		int secondBestScore = -300000000;
		int row = seq1.length();
		
		for (int c=1; c<=seq2.length(); c++) {
			if (matrix[row][c][M] > bestScore) {
				bestIdx = c;
				bestScore = matrix[row][c][M];
			} else if (matrix[row][c][M] > secondBestScore) {
				secondBestScore = matrix[row][c][M]; 
			}
		}
		
		int r = seq1.length();
		int c = bestIdx;
		int refEndIdx = c;
		
		List<Element> elems = new ArrayList<Element>();
		Element currElem = new Element('0', 0);
		
		int level = M;
		
		while (r > 0 && c > 0) {
			byte currBt = bt[r][c][level];
			
//			System.out.println(String.format("bt -- r:%d c:%d bt: %d level: %d, score: %d", r,c,currBt,level,matrix[r][c][level]));

			if (currBt == DIR_DIAG) {
				if (level == M) {
					r -= 1;
					c -= 1;					
				} else if (level == I) {
					r -= 1;
				} else if (level == D) {
					c -= 1;
				}
				
				if (level == M) {
					// If moving back to M level from I or D, skip update.
					currElem = updateCurrElem('M', currElem, elems);
				}
				
				level = M;

			} else if (currBt == DIR_LEFT) {
				if (level == D) {
					c -= 1;	
				} else if (level == M) {
					// noop
				}
				
				level = D;
				
				currElem = updateCurrElem('D', currElem, elems);
			} else if (currBt == DIR_UP) {
				if (level == I) {
					r -= 1;
				} else if (level == M) {
					// noop
				}
				
				level = I;
				
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

//		String ref = "ATCGAATTCCGGGCTA";
//		String seq = "GAACCCCTTCCG";
		
//		String ref = "ACTCG";
//		String seq = "ACCG";
		
//		String ref = "ATGCCC";
//		String seq = "ATGGCC";
		
//		String ref = "ATCGAATTCCGGGCTA";
//		String seq = "AATTCTA";

//		SemiGlobalAligner sg = new SemiGlobalAligner();
//		
//		Result res = sg.align(seq, ref);
//		System.out.println(res);
		
//		double i = 1;
//		i += -1e8;
		
//		System.out.println(i);

		
		Result res = null;
		SemiGlobalAligner sga = new SemiGlobalAligner(8, -32, -48, -1);
		
		long start = System.currentTimeMillis();
		
		for (int i=0; i<1000; i++) {
			String ref = "CCAGATCAGCCTAGGCAACATGGTGAAACCCCGTCTCTACCAAAAATAAAAAACTTAGCTGAGCGTGGTGGTGCACGCCTGTAGCCCCAGCTGCTGAGGAGCCTGAGCCCAGGGGGTGGAGGCTGCAGTGAGCCATGATCACACTACTGTACTCCAGCCTAGGTGACAGAGTGAGACCCTGTCTCAAAAAAATAAAAGAAAATAAAAATAAACAAAGAGAGAAGTGGAAGAAGAGGTGGAGTTTTGTATTTATGACTTGAATTTTGTATTCATGACTGGGTTGACACCCCAATCCACTCCATTTTTAGCCTTGAAACATGGCAAACAGTAACCATTAAAAGGATGGAAAAGAGAAGAAGGCATGGGTGGGAAACTGTGCCTCCCATTTTTGTGCATCTTTGTTGCTGTCCTTCCACTATACTGTACCTTTCAGCATTTTGACGGCAACCTGGATTGAGACTCCTGTTTTGCTAATTCCATAAGCTGTTGCGTTCATCACTTTTCCAAAAGCACCTGATCCTAGTACCTTCCCTGCAAAGACAAATGGTGAGTACGTGCATTTTAAAGATTTTCCAATGGAAAAGAAATGCTGCAGAAACATTTGGCACATTCCATTCTTACCAAACTCTAAATTTTCTCTTGGAAACTCCCATTTGAGATCATATTCATATTCTCTGAAATCAACGTAGAAGTACTCATTATCTGAGGAGCCGGTCACCTGTACCATCTGTAGCTGGCTTTCATACCTAAATTGCTTCAGAGATGAAATGATGAGTCAGTTAGGAATAGGCAGTTCTGCAGATAGAGGAAAGAATAATGAATTTTTACCTTTGCTTTTACCTTTTTGTACTTGTGACAAATTAGCAGGGTTAAAACGACAATGAAGAGGAGACAAACACCAATTGTTGCATAGAATGAGATGTTGTCTTGGATGAAAGGGAAGGGGCCTGCAACAAAAGAGTGTCACTCAGCGATGAAACAGAATTCCTGTGTGACATTATAAATAGTGGACAACTCATTATAATCTCTCACATCCTGTTTCAGTAATAATCATTTTCAGTCCTAACAACCACTCTACATATACTCTACTCCCCACAGACAATCAGGCAATGTCCCTGTAAAGGATACATTTCCTCCCTAGAAAATTGCGGATTATTCTCAATCCATTCTTTAAAACCATTTACTAGGGTAAATTTACAAGAATTACATCTGGTCCAGGCACGATGGCTCACGCCTGTAGTCCCAGCACTTTGGGAGGCCAAGATGGGAGGATCACTTGAGTCCAAGAATTAGACACCAGCCCAGGCAACACAGTGAAATCCCGTCTCTAAAAAAATTCAAAAATTAGCTGGGCGTGGTGGCAGGTGCCTGTAATCCCAGCTGCTCGGGAGGCTGAGGCAGGAG";
			String seq = "CCATTCTTACCAAACTCTAAATTTTCTCTTGGAAACTCCCATTTGAGATCATATTCATATTCTCTGAAATCAACGTAGAAGTACTCATTATCTGAGGAGCCGCACATTCCATTCTTGCCAAACTCTAGATTTTCTCTTGGAAACTCCCATTTGAGATCACATTCATATTCTCTGAAATCAACGTAGAAGTACTCATTATCTGAGGAGCCGGTCACCCGTACCATCTGTAGC";
			
			res = sga.align(seq, ref);
		}
		long end = System.currentTimeMillis();
		System.out.println("Elapsed: " + (end-start));
		System.out.println(String.format("score: %d, second: %d, position: %d, cigar: %s", res.score, res.secondBest, res.position, res.cigar));
//		Thread.sleep(300*1000);
	}
}

package abra;

import abra.SemiGlobalAligner.Result;

public class NativeSemiGlobalAligner {

	private native String align(String seq1, String seq2, int match, int mismatch, int gapOpen, int gapExtend);
	
	private int match = 8;
	private int mismatch = -32;
	private int gapOpen = -48;
	private int gapExtend = -1;
	
	public NativeSemiGlobalAligner(int match, int mismatch, int gapOpen, int gapExtend) {
		this.match = match;
		this.mismatch = mismatch;
		this.gapOpen = gapOpen;
		this.gapExtend = gapExtend;
	}
	
	public Result align(String seq1, String seq2) {
		
		// Result returned in format score:secondBest:pos:endPos:cigar
		// 789:741:611:734:52M108I71M
		String res = align(seq1, seq2, match, mismatch, gapOpen, gapExtend);
		
		String[] results = res.split(":");
		Result result = new Result(Integer.valueOf(results[0]), Integer.valueOf(results[1]),
				Integer.valueOf(results[2]), Integer.valueOf(results[3]), results[4]);
		
		return result;
	}
}

/* Copyright 2013 University of North Carolina at Chapel Hill.  All rights reserved. */
package abra;

import static abra.SequenceUtil.MatchResult.EXACT_MATCH;
import static abra.SequenceUtil.MatchResult.UNMATCHED;

public class SequenceUtil {

	public static MatchResult isMatch(String seq1, String seq2, int allowedMismatches) {
		
		if (allowedMismatches == 0) {
			return seq1.equals(seq2) ? EXACT_MATCH : UNMATCHED;
		}
		
		int mismatches = 0;
		
		if (seq1.length() != seq2.length()) {
			return UNMATCHED;
		}
		
		for (int i=0; i<seq1.length(); i++) {
			if (seq1.charAt(i) != seq2.charAt(i)) {
				mismatches += 1;
			}
			
			if (mismatches > allowedMismatches) {
				return UNMATCHED;
			}
		}
		
		return new MatchResult(true, mismatches);
	}
	
	public static class MatchResult {
		private boolean isMatch;
		private int numMismatches = -1;
		
		public static final MatchResult EXACT_MATCH = new MatchResult(true, 0);
		public static final MatchResult UNMATCHED = new MatchResult(false, -1);
		
		private MatchResult(boolean isMatch, int numMismatches) {
			this.isMatch = isMatch;
			this.numMismatches = numMismatches;
		}
		
		public boolean isMatch() {
			return isMatch;
		}
		
		public int getNumMismatches() {
			return numMismatches;
		}
	}
}

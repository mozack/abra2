package abra;

import java.util.List;

public class ConsensusSequence {

	// Mutates input!!
	// Assumes all input contigs same length
	public String buildConsensus(List<ScoredContig> contigs) {
		
		StringBuffer consensus = new StringBuffer();
		// Subset to the contigs with highest qual scores
		contigs = ScoredContig.filter(contigs, 31);
		
		int maxLen = 0;
		for (ScoredContig contig : contigs) {
			maxLen = Math.max(maxLen, contig.getContig().length());
		}
		
		for (int i=0; i<maxLen; i++) {
			int a = 0;
			int c = 0;
			int t = 0;
			int g = 0;
			
			for (ScoredContig contig : contigs) {
				if (i < contig.getContig().length()) {
					switch (contig.getContig().charAt(i)) {
						case 'A':
							a++;
							break;
						case 'C':
							c++;
							break;
						case 'T':
							t++;
							break;
						case 'G':
							g++;
							break;
						default:
							break;
					}
				}
			}
			
			int max = Math.max(g, Math.max(t, Math.max(a, c)));
			if (max == a) {
				consensus.append('A');
			} else if (max == c) {
				consensus.append('C');
			} else if (max == t) {
				consensus.append('T');
			} else if (max == g) {
				consensus.append('G');
			} else {
				consensus.append('N');
			}
		}
		
		return consensus.toString();
	}
}

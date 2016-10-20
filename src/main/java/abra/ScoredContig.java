package abra;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class ScoredContig implements Comparable<ScoredContig> {
	
	// TODO: Parameterize and optimize.
	public static int MAX_CONTIGS = 200;

	private double score;
	private String contig;
	
	public ScoredContig(double score, String contig) {
		this.score = score;
		this.contig = contig;
	}

	public double getScore() {
		return score;
	}

	public String getContig() {
		return contig;
	}

	@Override
	public int compareTo(ScoredContig o) {
		if (this.score < o.score) {
			return -1;
		} else if (this.score > o.score) {
			return 1;
		} else {
			return 0;
		}
	}
	
	public static List<ScoredContig> convertAndFilter(String contigStrings) {
		List<ScoredContig> contigs = new ArrayList<ScoredContig>();
		
		double score = Double.NEGATIVE_INFINITY;
		String[] contigSeq = contigStrings.split("\n");
		for (String str : contigSeq) {
			if (str.startsWith(">")) {
				String[] fields = str.split("_");
				score = Double.parseDouble(fields[4]);
			}
			
			contigs.add(new ScoredContig(score, str));
		}
		
		if (contigs.size() > MAX_CONTIGS) {
			System.err.println("Shrinking eligible contigs from " + contigs.size() + " to " + MAX_CONTIGS);
			Collections.sort(contigs);
			// Subset to only the first MAX_CONTIGS
			contigs = contigs.subList(0, MAX_CONTIGS);
		}
		
		return contigs;
	}
}

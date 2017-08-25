package abra;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class ScoredContig implements Comparable<ScoredContig> {
	
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
	
	public String toString() {
		return String.valueOf(score);
	}

	@Override
	public int compareTo(ScoredContig o) {
		if (this.score < o.score) {
			return 1;
		} else if (this.score > o.score) {
			return -1;
		} else {
			return 0;
		}
	}
	
//	public static List<ScoredContig> convertAndFilter(String contigStrings) {
//		return convertAndFilter(contigStrings, MAX_CONTIGS);
//	}
	
	public static List<ScoredContig> convertAndFilter(String contigStrings, int maxContigs, StringBuffer readBuffer) {
		List<ScoredContig> contigs = new ArrayList<ScoredContig>();
		
		double score = Double.NEGATIVE_INFINITY;
		String[] contigSeq = contigStrings.split("\n");
		for (String str : contigSeq) {
			if (str.startsWith(">")) {
				String[] fields = str.split("_");
				try {
					score = Double.parseDouble(fields[4]);
				} catch (ArrayIndexOutOfBoundsException e) {
					Logger.error("Error parsing assembled contigs.  Line: [" + str + "]\n\nContigs: [\n" + contigStrings + "\n]");
					Logger.error("Read buffer: [\n" + readBuffer.toString() + "\n]");
					throw e;
				}
			} else {
				contigs.add(new ScoredContig(score, str));
			}
		}
		
		return filter(contigs, maxContigs);
	}
	
	public static List<ScoredContig> filter(List<ScoredContig> contigs, int maxContigs) {
		if (contigs.size() > maxContigs) {
			Logger.debug("Shrinking eligible contigs from %d to %d", contigs.size(), maxContigs);
			Collections.shuffle(contigs);
			Collections.sort(contigs);
			// Subset to only the first MAX_CONTIGS
			contigs = contigs.subList(0, maxContigs);
		}
		
		return contigs;
	}
}

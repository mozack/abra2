package abra;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Simple ungapped alignment based upon initial exact seed match
 * @author lmose
 */
public class SimpleMapper {
	
	//TODO: Parameterize
	private static final int KMER_SIZE = 25;
	private static final int MAX_MISMATCHES = 5;
	
	static final int UNMAPPED = -1;
	static final int HOMOLOGOUS_MAPPING = -2;
	
	// Represents a single contig
	private String ref;
	private String refRc;
	
	// kmer -> list of positions
	private Map<String, List<Integer>> kmerPositions = new HashMap<String, List<Integer>>();

	public SimpleMapper(String ref) {
		this.ref = ref;
		for (int i=0; i<ref.length()-KMER_SIZE; i++) {
			String kmer = ref.substring(i, i+KMER_SIZE);
			if (!kmerPositions.containsKey(kmer)) {
				kmerPositions.put(kmer, new ArrayList<Integer>());
			}
			kmerPositions.get(kmer).add(i);
		}
	}
	
	//TODO: reverse complement?
	public SimpleMapperResult map(String read) {
		
		// ref position -> # mismatches
		Map<Integer, Integer> posMismatches = new HashMap<Integer, Integer>();
		
		// Evaluate all reference positions where read and ref share a kmer
		for (int i=0; i<read.length()-KMER_SIZE; i++) {
			String kmer = read.substring(i, i+KMER_SIZE);
			
			if (kmerPositions.containsKey(kmer)) {
				for (int pos : kmerPositions.get(kmer)) {
					int refStartPos = pos - i;
					
					// Compare strings only if this position has not already been evaluated
					if (refStartPos >= 0 && !posMismatches.containsKey(refStartPos) && refStartPos < ref.length() - read.length()) {
						int mismatches = countMismatches(refStartPos, read);
						posMismatches.put(refStartPos, mismatches);
					}
				}
			}
		}
		
		// Find position with fewest mismatches
		int bestMismatches = read.length() + 1;
		int bestPos = UNMAPPED;
		
		for (int pos : posMismatches.keySet()) {
			if (posMismatches.get(pos) < bestMismatches) {
				bestMismatches = posMismatches.get(pos);
				bestPos = pos;
			} else if (posMismatches.get(pos) == bestMismatches) {
				bestPos = HOMOLOGOUS_MAPPING;
			}
		}
		
		if (bestMismatches > MAX_MISMATCHES) {
			bestPos = UNMAPPED;
		}
		
		return new SimpleMapperResult(bestPos, bestMismatches);
	}
	
	private int countMismatches(int refPosition, String read) {
		int mismatches = 0;
		for (int i=0; i<read.length(); i++) {
			if (read.charAt(i) != ref.charAt(refPosition+i)) {
				mismatches += 1;
				if (mismatches > MAX_MISMATCHES) {
					break;
				}
			}
		}
		
		return mismatches;
	}
	
	static class SimpleMapperResult {
		private int pos;
		private int mismatches;
		
		SimpleMapperResult(int pos, int mismatches) {
			this.pos = pos;
			this.mismatches = mismatches;
		}
		
		public int getPos() {
			return pos;
		}
		public int getMismatches() {
			return mismatches;
		}
	}
	
	public static void main(String[] args) {
		String contig = "TTCAACTAGAGAGAGGTAAAAATTTTTCTAGAACATGAATTGCCCACTCCCCTCATTCCTTCTCAGAAACTAACTGAATTCCAGTGGGTGTGCCTGGCAAACCCAAAAGCAGTTTCTGTTCAGGATGCTGGTCTTACCTGTGAAGGCGTTCATGAACGTGGAGAGGGACCGGTTCAACATTTTGAAGAAAGGGTCTCTGCACGGATATTTCTGAGACCCACAAAGGACGGTATGCTCAAGAATGTGAGGAACACCAGTACTGTCCATGGGAGTGGTACGGAACTGCACGCTAGGGAAGAGAGAGGAATGGCACGCTAGGGAAGGCGAATGACCAGAACGCAAAAGGTTCAGCTTAGTGCTGCGGACACAGTTCCCAGATGCATCATCACCTCAGGCTACTAGAAATCATCATTCTGACACCACAATCCTCCAGCACAGGGTTTTCCAACTATA";
		String read = "TACTGTCCATGGGAGTGGTACGGAACTGCACGCTAGGGAAGAGAGAGGAATGGCACGCTAGGGAAGGCGAATGACCAGAACGCAAAAGGTTCAGCTTAGTG";
		
		SimpleMapper sm = new SimpleMapper(contig);
		
		SimpleMapperResult result = sm.map(read);
		System.out.println("pos: " + result.getPos() + ". mismatches: " + result.getMismatches());
	}
}

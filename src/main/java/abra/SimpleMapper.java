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
	
	private static final int KMER_SIZE = 25;
	
	static final int UNMAPPED = -1;
	static final int HOMOLOGOUS_MAPPING = -2;
	
	static final char FORWARD_ORIENTATION = 0;
	static final char REVERSE_ORIENTATION = 1;
	
	enum Orientation {
		UNSET, FORWARD, REVERSE;
	}
	
	// Represents a single contig
	private String ref;
	private ReverseComplementor rc = new ReverseComplementor();
	private double maxMismatchRate;
	
	// kmer -> list of positions
	private Map<String, List<Integer>> kmerPositions = new HashMap<String, List<Integer>>();
	
	public SimpleMapper(String ref) {
		this(ref, .05);
	}
	
	public SimpleMapper(String ref, double maxMismatchRate) {
		this.ref = ref;
		this.maxMismatchRate = maxMismatchRate;
		
		for (int i=0; i<ref.length()-KMER_SIZE; i++) {
			String kmer = ref.substring(i, i+KMER_SIZE);
			if (!kmerPositions.containsKey(kmer)) {
				kmerPositions.put(kmer, new ArrayList<Integer>());
			}
			kmerPositions.get(kmer).add(i);
		}
	}
	
	private Map<Integer, Integer> getPositionMismatches(String bases) {
		// ref position -> # mismatches
		Map<Integer, Integer> posMismatches = new HashMap<Integer, Integer>();
		
		// Evaluate all reference positions where read and ref share a kmer
		for (int i=0; i<bases.length()-KMER_SIZE; i++) {
			String kmer = bases.substring(i, i+KMER_SIZE);
			
			if (kmerPositions.containsKey(kmer)) {
				for (int pos : kmerPositions.get(kmer)) {
					int refStartPos = pos - i;
					
					// Compare strings only if this position has not already been evaluated
					if (refStartPos >= 0 && !posMismatches.containsKey(refStartPos) && refStartPos <= ref.length() - bases.length()) {
						int mismatches = countMismatches(refStartPos, bases);
						posMismatches.put(refStartPos, mismatches);
					}
				}
			}
		}

		return posMismatches;
	}
	
	public SimpleMapperResult map(String read) {
		
		Map<Integer, Integer> forwardMismatches = getPositionMismatches(read);
		Map<Integer, Integer> reverseMismatches = getPositionMismatches(rc.reverseComplement(read));
		
		// Find position with fewest mismatches
		int bestMismatches = read.length() + 1;
		int bestPos = UNMAPPED;
		Orientation bestOrientation = Orientation.UNSET;
		
		// Search for matches to contig in forward orientation
		for (int pos : forwardMismatches.keySet()) {
			if (forwardMismatches.get(pos) < bestMismatches) {
				bestMismatches = forwardMismatches.get(pos);
				bestPos = pos;
				bestOrientation = Orientation.FORWARD;;
			} else if (forwardMismatches.get(pos) == bestMismatches) {
				bestPos = HOMOLOGOUS_MAPPING;
			}
		}
		
		// Search for matches to contig in reverse complement
		for (int pos : reverseMismatches.keySet()) {
			if (reverseMismatches.get(pos) < bestMismatches) {
				bestMismatches = reverseMismatches.get(pos);
				bestPos = pos;
				bestOrientation = Orientation.REVERSE;
			} else if (reverseMismatches.get(pos) == bestMismatches) {
				bestPos = HOMOLOGOUS_MAPPING;
			}
		}
		
		if (bestMismatches > read.length() * maxMismatchRate) {
			bestPos = UNMAPPED;
		}
		
		return new SimpleMapperResult(bestPos, bestMismatches, bestOrientation);
	}
	
	private int countMismatches(int refPosition, String read) {
		int mismatches = 0;
		for (int i=0; i<read.length(); i++) {
			if (read.charAt(i) != ref.charAt(refPosition+i)) {
				mismatches += 1;
				if (mismatches > read.length() * maxMismatchRate) {
					break;
				}
			}
		}
		
		return mismatches;
	}
	
	static class SimpleMapperResult {
		private int pos;
		private int mismatches;
		private Orientation orientation;
		
		SimpleMapperResult(int pos, int mismatches, Orientation orientation) {
			this.pos = pos;
			this.mismatches = mismatches;
			this.orientation = orientation;
		}
		
		public int getPos() {
			return pos;
		}
		public int getMismatches() {
			return mismatches;
		}
		public Orientation getOrientation() {
			return orientation;
		}
	}
	
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((ref == null) ? 0 : ref.hashCode());
		return result;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		SimpleMapper other = (SimpleMapper) obj;
		if (ref == null) {
			if (other.ref != null)
				return false;
		} else if (!ref.equals(other.ref))
			return false;
		return true;
	}
	
	public String toString() {
		return ref;
	}

	public static void main(String[] args) {
		String contig = "TTCAACTAGAGAGAGGTAAAAATTTTTCTAGAACATGAATTGCCCACTCCCCTCATTCCTTCTCAGAAACTAACTGAATTCCAGTGGGTGTGCCTGGCAAACCCAAAAGCAGTTTCTGTTCAGGATGCTGGTCTTACCTGTGAAGGCGTTCATGAACGTGGAGAGGGACCGGTTCAACATTTTGAAGAAAGGGTCTCTGCACGGATATTTCTGAGACCCACAAAGGACGGTATGCTCAAGAATGTGAGGAACACCAGTACTGTCCATGGGAGTGGTACGGAACTGCACGCTAGGGAAGAGAGAGGAATGGCACGCTAGGGAAGGCGAATGACCAGAACGCAAAAGGTTCAGCTTAGTGCTGCGGACACAGTTCCCAGATGCATCATCACCTCAGGCTACTAGAAATCATCATTCTGACACCACAATCCTCCAGCACAGGGTTTTCCAACTATA";
		String read = "TACTGTCCATGGGAGTGGTACGGAACTGCACGCTAGGGAAGAGAGAGGAATGGCACGCTAGGGAAGGCGAATGACCAGAACGCAAAAGGTTCAGCTTAGTG";
		
		SimpleMapper sm = new SimpleMapper(contig, .05);
		
		SimpleMapperResult result = sm.map(read);
		System.out.println("pos: " + result.getPos() + ". mismatches: " + result.getMismatches());
	}
}

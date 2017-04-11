package abra;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

/**
 * Simple ungapped alignment based upon initial exact seed match
 * @author lmose
 */
public class SimpleMapper {
	
	//TODO: May need to change for short reads or decreased maxMismatchRate
	//      Increasing when appropriate will speed things up.
	private static final int KMER_SIZE = 15;
	
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
	
	private Map<Integer, Integer> getPositionMismatches_old(String bases) {
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
	
	private Map<Integer, Integer> getPositionMismatches(String bases) {
		// ref position -> # mismatches
		Map<Integer, Integer> posMismatches = new HashMap<Integer, Integer>();
		
		// pos in ref -> kmer idx
		// i.e. For pos 100, kmers 0,2,3 (0,50,75 if KMER_SIZE=25) match
		
		Map<Integer, List<Integer>> refPosKmers = new HashMap<Integer, List<Integer>>();
		
		// Identify kmers in read that match reference
		int i=0;
		while (i<=bases.length()-KMER_SIZE) {		
			String kmer = bases.substring(i, i+KMER_SIZE);
			
			if (kmerPositions.containsKey(kmer)) {
				for (int pos : kmerPositions.get(kmer)) {
					int refStartPos = pos - i;
					if (refStartPos > 0 && refStartPos <= ref.length() - bases.length()) {
						if (!refPosKmers.containsKey(refStartPos)) {
							refPosKmers.put(refStartPos, new ArrayList<Integer>());
						}
						refPosKmers.get(refStartPos).add(i);
					}
				}
			}
			
			// Skip to next kmer
			if (i+KMER_SIZE <= bases.length()-KMER_SIZE) {
				i += KMER_SIZE;
			} else if (i < bases.length()-KMER_SIZE){
				// Last kmer may overlap previous kmer
				i = bases.length()-KMER_SIZE;
			} else {
				// Done.
				i = bases.length();
			}
		}
		
		// Count mismatches in read.  Skip over matching kmers for speed.
		for (int refStartPos : refPosKmers.keySet()) {
			int mismatches = 0;
			List<Integer> matchingKmers = refPosKmers.get(refStartPos);
			Iterator<Integer> kmerIter = matchingKmers.iterator();
			int currKmer = -1;
			if (kmerIter.hasNext()) {
				currKmer = kmerIter.next();
			}
			
			i=0;
			while (i<bases.length()) {
				
				if (i >= currKmer && i < currKmer+KMER_SIZE) {
					if (kmerIter.hasNext()) {
						currKmer = kmerIter.next();
					}
					
					// Skip over the matching kmer
					i += KMER_SIZE;
				} else {
					if (bases.charAt(i) != ref.charAt(refStartPos+i)) {
						mismatches += 1;
						if (mismatches > bases.length() * maxMismatchRate) {
							break;
						}
					}
					
					i += 1;
				}
			}

			posMismatches.put(refStartPos, mismatches);
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
	
	public String getSeq() {
		return this.ref;
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
//		String contig = "TTCAACTAGAGAGAGGTAAAAATTTTTCTAGAACATGAATTGCCCACTCCCCTCATTCCTTCTCAGAAACTAACTGAATTCCAGTGGGTGTGCCTGGCAAACCCAAAAGCAGTTTCTGTTCAGGATGCTGGTCTTACCTGTGAAGGCGTTCATGAACGTGGAGAGGGACCGGTTCAACATTTTGAAGAAAGGGTCTCTGCACGGATATTTCTGAGACCCACAAAGGACGGTATGCTCAAGAATGTGAGGAACACCAGTACTGTCCATGGGAGTGGTACGGAACTGCACGCTAGGGAAGAGAGAGGAATGGCACGCTAGGGAAGGCGAATGACCAGAACGCAAAAGGTTCAGCTTAGTGCTGCGGACACAGTTCCCAGATGCATCATCACCTCAGGCTACTAGAAATCATCATTCTGACACCACAATCCTCCAGCACAGGGTTTTCCAACTATA";
//		String read = "TACTGTCCATGGGAGTGGTACGGAACTGCACGCTAGGGAAGAGAGAGGAATGGCACGCTAGGGAAGGCGAATGACCAGAACGCAAAAGGTTCAGCTTAGTG";
		
//		SimpleMapper sm = new SimpleMapper(contig, .05);
		
		String contig1 = "TTCAACTAGAGAGAGGTAAAAATTTTTCTAGAACATGAATTGCCCACTCCCCTCATTCCTTCTCAGAAACTAACTGAATTCCAGTGGGTGTGCCTGGCAAACCCAAAAGCAGTTTCTGTTCAGGATGCTGGTCTTACCTGTGAAGGCGTTCATGAACGTGGAGAGGGACCGGTTCAACATTTTGAAGAAAGGGTCTCTGCACGGATATTTCTGAGACCCACAAAGGACGGTATGCTCAAGAATGTGAGGAACACCAGTACTGTCCATGGGAGTGGTACGGAACTGCACGCTAGGGAAGAGAGAGGAATGGCACGCTAGGGAAGGCGAATGACCAGAACGCAAAAGGTTCAGCTTAGTGCTGCGGACACAGTTCCCAGATGCATCATCACCTCAGGCTACTAGAAATCATCATTCTGACACCACAATCCTCCAGCACAGGGTTTTCCAACTATA";
		String read = "CCTGAACAGAAACTGCTTTTGGGAATGCCAGGCACACCCACTGGAATTCAGTTAGTTTCTGAGAAGGAATGAGGGGAGTGGGCAATTCATGTTCTAGAAA";

		SimpleMapper sm = new SimpleMapper(contig1, .05);
		
		SimpleMapperResult result = null;
		
		long start = System.currentTimeMillis();
		for (int i=0; i<100000; i++) {
			result = sm.map(read);
		}
		long end = System.currentTimeMillis();
		System.out.println("Elapsed: " + (end-start));
		
		System.out.println("pos: " + result.getPos() + ". mismatches: " + result.getMismatches());
	}
}

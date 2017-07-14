package abra;

import java.util.ArrayList;
import java.util.List;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.TextCigarCodec;

public class ContigAligner {
	
	private static final int MIN_ALIGNMENT_SCORE = 1;

	private static int MATCH;
	private static int MISMATCH;
	private static int GAP_OPEN_PENALTY;
	private static int GAP_EXTEND_PENALTY;
	
	private String refChr;
	private int refContextStart;
	String ref;
	private int minAnchorLength;
	private int maxAnchorMismatches;
	
	private List<Integer> junctionPositions = new ArrayList<Integer>();
	private List<Integer> junctionLengths = new ArrayList<Integer>();
	
	private CompareToReference2 localC2r;
	
	private NativeSemiGlobalAligner aligner = new NativeSemiGlobalAligner(MATCH, MISMATCH, GAP_OPEN_PENALTY, GAP_EXTEND_PENALTY);
	
	public static void init(int[] scoring) {
		
		for (int i=0; i<scoring.length; i++) {
			if (scoring[i] < 0) {
				String msg = "Please specify all Smith Waterman scores as positive values";
				Logger.error(msg);
				throw new RuntimeException(msg);
			}
		}
		
		MATCH = scoring[0];
		MISMATCH = -scoring[1];
		GAP_OPEN_PENALTY = -scoring[2];
		GAP_EXTEND_PENALTY = -scoring[3];
		
		Logger.info("SG match,mismatch,gap_open_penalty,gap_extend_penalty: " 
				+ MATCH + "," + MISMATCH + "," + GAP_OPEN_PENALTY + "," + GAP_EXTEND_PENALTY);
	}
	
	public ContigAligner(String ref, String refChr, int refStart, int minContigLength, int minAnchorLen, int maxAnchorMismatches) {
		this.ref = ref;
		this.refChr = refChr;
		this.refContextStart = refStart;
		this.minAnchorLength = minAnchorLen;
		this.maxAnchorMismatches = maxAnchorMismatches;
		this.localC2r = new CompareToReference2();
		localC2r.initLocal(refChr, ref);
	}
		
	public ContigAligner(String ref, String refChr, int refStart, int minContigLength, int minAnchorLen, int maxAnchorMismatches, 
			List<Integer> junctionPositions, List<Integer> junctionLengths) {
		this(ref, refChr, refStart, minContigLength, minAnchorLen, maxAnchorMismatches);
		this.junctionPositions = junctionPositions;
		this.junctionLengths = junctionLengths;
	}
	
	public ContigAlignerResult align(String seq) {
		
		ContigAlignerResult result = null;
		
		SemiGlobalAligner.Result sgResult = aligner.align(seq, ref);
		Logger.trace("SG Alignment [%s]:\t%s, possible: %d to: %s", seq, sgResult, seq.length()*MATCH, ref);
		if (sgResult.score > MIN_ALIGNMENT_SCORE && sgResult.score > sgResult.secondBest && sgResult.endPosition > 0) {
			Cigar cigar = TextCigarCodec.decode(sgResult.cigar);
			
			CigarElement first = cigar.getFirstCigarElement();
			CigarElement last = cigar.getLastCigarElement();
			
			// Do not allow indels at the edges of contigs.
			if (minAnchorLength > 0 && 
				(first.getOperator() != CigarOperator.M || first.getLength() < minAnchorLength || 
				last.getOperator() != CigarOperator.M || last.getLength() < minAnchorLength)) {

//				if ((first.getOperator() != CigarOperator.M || last.getOperator() != CigarOperator.M) &&
//						cigar.toString().contains("I")) {
//					Logger.trace("INDEL_NEAR_END: %s", cigar.toString());
//					return ContigAlignerResult.INDEL_NEAR_END;
//				}
//				
//				if ((first.getLength() < 5 || last.getLength() < 5) && 
//						cigar.toString().contains("I") &&
//						minAnchorLength >= 5) {
//					Logger.trace("INDEL_NEAR_END: %s", cigar.toString());
//					return ContigAlignerResult.INDEL_NEAR_END;						
//				}

				return null;
			}
				
			int endPos = sgResult.position + cigar.getReferenceLength();
			
			// Require first and last minAnchorLength bases of contig to be similar to ref
			int mismatches = 0;
			for (int i=0; i<minAnchorLength; i++) {
				if (seq.charAt(i) != ref.charAt(sgResult.position+i)) {
					mismatches += 1;
				}
			}
			
			if (mismatches > maxAnchorMismatches) {
				Logger.trace("Mismatches at beginning of: %s", seq);
			} else {
			
				mismatches = 0;
				for (int i=minAnchorLength; i>0; i--) {
					
					int seqIdx = seq.length()-i;
					int refIdx = endPos-i;
					
					if (seq.charAt(seqIdx) != ref.charAt(refIdx)) {
						mismatches += 1;
					}
				}
				
				if (mismatches > maxAnchorMismatches) {
					Logger.trace("Mismatches at end of: %s", seq);
				} else {
					result = finishAlignment(sgResult.position, endPos, sgResult.cigar, sgResult.score, seq);
				}
			}
		}
		
		return result;
	}
	
	ContigAlignerResult finishAlignment(int refStart, int refEnd, String alignedCigar, int score, String seq) {
		try {
			// Pad with remaining reference sequence
			String leftPad = ref.substring(0, refStart);
			String rightPad = "";
			if (refEnd < ref.length()-1) {
				rightPad = ref.substring(refEnd,ref.length());
			}
			String paddedSeq = leftPad + seq + rightPad;
			String cigar = CigarUtils.extendCigarWithMatches(alignedCigar, leftPad.length(), rightPad.length());
			Logger.trace("Padded contig: %s\t%s", cigar, paddedSeq);
			
			if (junctionPositions.size() > 0) {
				String oldCigar = cigar;
				cigar = CigarUtils.injectSplices(cigar, junctionPositions, junctionLengths);
				Logger.trace("Spliced Cigar.  old: %s, new: %s", oldCigar, cigar);
			}
			
			return new ContigAlignerResult(refStart-leftPad.length(), cigar, refChr, refContextStart, paddedSeq, score);
		} catch (StringIndexOutOfBoundsException e) {
			e.printStackTrace();
			System.err.println(String.format("index error: %d, %d, %s, %d, %s, %s", refStart, refEnd, alignedCigar, score, seq, ref));
			
			throw e;
		}
	}
	
	public static class ContigAlignerResult {
		
		// Used for testing
		static boolean PAD_CONTIG = true;
		
		private int localRefPos;
		private String cigar;
		
		private String chromosome;
		private int refContextStart;
		
		private String sequence;
		private int score;
		private boolean isSecondary = false;
		
		public static final ContigAlignerResult INDEL_NEAR_END = new ContigAlignerResult();

		private ContigAlignerResult() {
		}
		
		ContigAlignerResult(int refPos, String cigar, String chromosome, int refContextStart, String sequence, int score) {
			this.localRefPos = refPos;
			this.cigar = cigar;
			this.chromosome = chromosome;
			this.refContextStart = refContextStart;
			this.sequence = sequence;
			this.score = score;
		}
		
		public int getRefPos() {
			return localRefPos;
		}
		public String getCigar() {
			return cigar;
		}

		public String getChromosome() {
			return chromosome;
		}

		public int getRefContextStart() {
			return refContextStart;
		}
		
		// This is the actual genomic position 
		public int getGenomicPos() {
			return localRefPos + refContextStart;
		}
		
		public String getSequence() {
			return sequence;
		}
		
		public int getScore() {
			return score;
		}
		
		public boolean isSecondary() {
			return isSecondary;
		}

		public void setSecondary(boolean isSecondary) {
			this.isSecondary = isSecondary;
		}
	}
	
	public static void main(String[] args) {
		String ref = "AACAACAGATAATAACAAGTCCTAACCCTCTAGCTGCTTAGGCTGGCGGAGGCCCAGGGGCTCCCACGAGTTGGGTCCTTTCGCACCAGCACAGACTTACCTGATCTCGGTTGTTGATGTGAGAATAAGGAAGCTCCCCCGTCATCAGTTCATACAATACGATGCCATAGGAGTAGACATCCGACTGGAAACTGAATGGGTTGTTATCCTGCATTCGGATCACCTCTGGGGCCTACATGTATCACCATATGACAAAAGTGCATTTATCACCATATGACAGGCCTCACAGACATCTAGGGGCCAGGCTGTCCCTTTCATTAGTTATGAATGAG";
		String contig = "AACAACAGATAATAACAAGTCCTAACCCTCTAGCTGCTTAGGCTGGCGGAGGCCCAGGGGCTCCCACGAGTTGGGTCCTTTCGCACCAGCACAGACTTACCTGATCTCGGTTGTTGATGTGAGAATAAGGAAGCTCCCCCGTCATCAGTTCACAAAAGTGCATTTATCACCATATGACAGGCCTCACAGACATCTAGGGGCCAGGCTGTCCCTTTCATTAGTTATGAAT";
		
//		SSWAligner sw = new SSWAligner(ref, "chr3", 12626521, 50);
//		sw.align(contig);
	}
}

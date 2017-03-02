package abra;

import java.util.ArrayList;
import java.util.List;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.TextCigarCodec;
import ssw.Aligner;
import ssw.Alignment;

public class SSWAligner {
	
	private static final int MIN_ALIGNMENT_SCORE = 1;

//	private static final int MATCH = 20;
//	private static final int MISMATCH = -20;
//	private static final int GAP_OPEN_PENALTY = 41;
//	private static final int GAP_EXTEND_PENALTY = 1;
	
	// TODO: Optimize SW scoring
	// bwa scores
//	private static final int MATCH = 1;
//	private static final int MISMATCH = -4;
//	private static final int GAP_OPEN_PENALTY = 6;
//	private static final int GAP_EXTEND_PENALTY = 0;
	
	private static int MATCH;
	private static int MISMATCH;
	private static int GAP_OPEN_PENALTY;
	private static int GAP_EXTEND_PENALTY;
	
	private String refChr;
	private int refContextStart;
	private String ref;
	private int minContigLength;
	
	private List<Integer> junctionPositions = new ArrayList<Integer>();
	private List<Integer> junctionLengths = new ArrayList<Integer>();
	
	private static int [][] score;
	
	private boolean useSemiGlobal = true;
	private IndelShifter indelShifter = new IndelShifter();
	private CompareToReference2 c2r;
	
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
		
		Logger.info("SW match,mismatch,gap_open_penalty,gap_extend_penalty: " 
				+ MATCH + "," + MISMATCH + "," + GAP_OPEN_PENALTY + "," + GAP_EXTEND_PENALTY);
		
		score = new int[128][128];
		for (int i = 0; i < 128; i++) {
			for (int j = 0; j < 128; j++) {
				if (i == j) score[i][j] = MATCH;
				else score[i][j] = MISMATCH;
			}
		}
	}
	
	public SSWAligner(CompareToReference2 c2r, String ref, String refChr, int refStart, int minContigLength) {
		this.ref = ref;
		this.refChr = refChr;
		this.refContextStart = refStart;
		this.minContigLength = minContigLength;
		this.c2r = c2r;
	}
	
	public SSWAligner(CompareToReference2 c2r, String ref, String refChr, int refStart, int minContigLength, int junctionPos, int junctionLength) {
		this(c2r, ref, refChr, refStart, minContigLength);
		this.junctionPositions.add(junctionPos);
		this.junctionLengths.add(junctionLength);
	}
	
	public SSWAligner(CompareToReference2 c2r, String ref, String refChr, int refStart, int minContigLength, List<Integer> junctionPositions, List<Integer> junctionLengths) {
		this(c2r, ref, refChr, refStart, minContigLength);
		this.junctionPositions = junctionPositions;
		this.junctionLengths = junctionLengths;
	}
	
	public SSWAlignerResult align(String seq) {
		
		SSWAlignerResult result = null;
		
		if (useSemiGlobal) {
			SemiGlobalAligner aligner = new SemiGlobalAligner(MATCH, MISMATCH, GAP_OPEN_PENALTY, GAP_EXTEND_PENALTY);
			SemiGlobalAligner.Result sgResult = aligner.align(seq, ref);
			Logger.trace("SG Alignment [%s]:\t%s", seq, sgResult);
			if (sgResult.score > MIN_ALIGNMENT_SCORE && sgResult.score > sgResult.secondBest && sgResult.endPosition > 0) {
//			if (sgResult.score > MIN_ALIGNMENT_SCORE && sgResult.score > sgResult.secondBest) {
				Cigar cigar = TextCigarCodec.decode(sgResult.cigar);
				
				CigarElement first = cigar.getFirstCigarElement();
				CigarElement last = cigar.getLastCigarElement();
				
				// Do not allow indels at the edges of contigs.
				if (first.getOperator() == CigarOperator.M && first.getLength() >= 5 && 
					last.getOperator() == CigarOperator.M && last.getLength() >=5) {
					
					int endPos = sgResult.position + cigar.getReferenceLength();
					
//					cigar = indelShifter.shiftIndelsLeft(sgResult.position+this.refContextStart, endPos+this.refContextStart,
//							this.refChr, cigar, seq, c2r);
					
					cigar = indelShifter.shiftAllIndelsLeft(sgResult.position+this.refContextStart, endPos+this.refContextStart,
							this.refChr, cigar, seq, c2r);

					
					first = cigar.getFirstCigarElement();
					last = cigar.getLastCigarElement();
					
					String textCigar = TextCigarCodec.encode(cigar);
					Logger.trace("OLD_CIGAR: %s\tNEW_CIGAR%s", sgResult.cigar, textCigar);
					
					//TODO: Just do this once?
					if (first.getOperator() == CigarOperator.M && first.getLength() >= 5 && 
							last.getOperator() == CigarOperator.M && last.getLength() >=5) {
						
						result = finishAlignment(sgResult.position, endPos, textCigar, sgResult.score, seq);
					}
				}
			}
		} else {
			// Require minimum of minContigLength or 90% of the input sequence to align
			int minContigLen = Math.min(minContigLength, (int) (seq.length() * .9));
			Alignment aln = Aligner.align(seq.getBytes(), ref.getBytes(), score, GAP_OPEN_PENALTY, GAP_EXTEND_PENALTY, true);
			
			Logger.trace("Alignment [%s]:\t%s", seq, aln);
			
			// TODO: Optimize score requirements..
			if (aln != null && aln.score1 >= MIN_ALIGNMENT_SCORE && aln.score1 > aln.score2 && aln.read_end1 - aln.read_begin1 > minContigLen) {
							
				// Clip contig and remap if needed.
				// TODO: Trim sequence instead of incurring overhead of remapping
				
				int MAX_CLIP_BASES = Math.min(10, seq.length() / 10);
				if ((aln.read_begin1 > 0 || aln.read_end1 < seq.length()-1) &&
					(aln.read_begin1 < MAX_CLIP_BASES && aln.read_end1 > seq.length()-1-MAX_CLIP_BASES)) {
					
					seq = seq.substring(aln.read_begin1, aln.read_end1+1);
					aln = Aligner.align(seq.getBytes(), ref.getBytes(), score, GAP_OPEN_PENALTY, GAP_EXTEND_PENALTY, true);
				}
							
				// Requiring end to end alignment here...
				if (aln.read_begin1 == 0 && aln.read_end1 == seq.length()-1) {
					
					result = finishAlignment(aln.ref_begin1, aln.ref_end1, aln.cigar, aln.score1, seq);
				}
			}
		}
		
		return result;
	}
	
	SSWAlignerResult finishAlignment(int refStart, int refEnd, String alignedCigar, int score, String seq) {
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
			
			return new SSWAlignerResult(refStart-leftPad.length(), cigar, refChr, refContextStart, paddedSeq, score);
		} catch (StringIndexOutOfBoundsException e) {
			e.printStackTrace();
			System.err.println(String.format("index error: %d, %d, %s, %d, %s, %s", refStart, refEnd, alignedCigar, score, seq, ref));
			
			throw e;
		}
	}
	
	public static class SSWAlignerResult {
		
		// Used for testing
		static boolean PAD_CONTIG = true;
		
		private int localRefPos;
		private String cigar;
		
		private String chromosome;
		private int refContextStart;
		
		private String sequence;
		private int score;
		
		SSWAlignerResult(int refPos, String cigar, String chromosome, int refContextStart, String sequence, int score) {
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
	}
	
	public static void main(String[] args) {
		String ref = "AACAACAGATAATAACAAGTCCTAACCCTCTAGCTGCTTAGGCTGGCGGAGGCCCAGGGGCTCCCACGAGTTGGGTCCTTTCGCACCAGCACAGACTTACCTGATCTCGGTTGTTGATGTGAGAATAAGGAAGCTCCCCCGTCATCAGTTCATACAATACGATGCCATAGGAGTAGACATCCGACTGGAAACTGAATGGGTTGTTATCCTGCATTCGGATCACCTCTGGGGCCTACATGTATCACCATATGACAAAAGTGCATTTATCACCATATGACAGGCCTCACAGACATCTAGGGGCCAGGCTGTCCCTTTCATTAGTTATGAATGAG";
		String contig = "AACAACAGATAATAACAAGTCCTAACCCTCTAGCTGCTTAGGCTGGCGGAGGCCCAGGGGCTCCCACGAGTTGGGTCCTTTCGCACCAGCACAGACTTACCTGATCTCGGTTGTTGATGTGAGAATAAGGAAGCTCCCCCGTCATCAGTTCACAAAAGTGCATTTATCACCATATGACAGGCCTCACAGACATCTAGGGGCCAGGCTGTCCCTTTCATTAGTTATGAAT";
		
//		SSWAligner sw = new SSWAligner(ref, "chr3", 12626521, 50);
//		sw.align(contig);
	}
}

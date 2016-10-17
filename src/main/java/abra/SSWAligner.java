package abra;

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
	private static final int MATCH = 1;
	private static final int MISMATCH = -4;
	private static final int GAP_OPEN_PENALTY = 6;
	private static final int GAP_EXTEND_PENALTY = 0;
	
	private String refChr;
	private int refStart;
	private String ref;
	private int junctionPos;
	private int junctionLength;
	
	private static int [][] score;
	static {
		try {
			System.loadLibrary("sswjni");
		} catch (java.lang.UnsatisfiedLinkError e) {
			System.out.println(String.format("Cannot find libsswjni.so. Has the library been built and LD_LIBRARY_PATH or -Djava.library.path set appropriately?\n%s", e));
			throw e;
		}
		
		score = new int[128][128];
		for (int i = 0; i < 128; i++) {
			for (int j = 0; j < 128; j++) {
				if (i == j) score[i][j] = MATCH;
				else score[i][j] = MISMATCH;
			}
		}
	}
	
	public SSWAligner(String ref, String refChr, int refStart) {
		this(ref, refChr, refStart, -1, -1);
	}
	
	public SSWAligner(String ref, String refChr, int refStart, int junctionPos, int junctionLength) {
		this.ref = ref;
		this.refChr = refChr;
		this.refStart = refStart;
		this.junctionPos = junctionPos;
		this.junctionLength = junctionLength;
	}
	
	public SSWAlignerResult align(String seq) {
		
		SSWAlignerResult result = null;
		
		Alignment aln = Aligner.align(seq.getBytes(), ref.getBytes(), score, GAP_OPEN_PENALTY, GAP_EXTEND_PENALTY, true);
		
		System.err.println("Alignment [" +  seq + "] :\t" + aln);
		
		// TODO: Optimize score requirements..
		if (aln != null && aln.score1 >= MIN_ALIGNMENT_SCORE && aln.score1 > aln.score2) {
			
			// TODO: Is this best number?
			int MAX_CLIP_BASES = 10;
						
			// Clip contig and remap if needed.
			// TODO: Trim sequence and Cigar instead of incurring overhead of remapping
			if ((aln.read_begin1 > 0 && aln.read_begin1 < MAX_CLIP_BASES) ||
				(aln.read_end1 < seq.length()-1 && aln.read_end1 > seq.length()-1-MAX_CLIP_BASES)) {
				
				seq = seq.substring(aln.read_begin1, aln.read_end1+1);
				aln = Aligner.align(seq.getBytes(), ref.getBytes(), score, GAP_OPEN_PENALTY, GAP_EXTEND_PENALTY, true);
				System.err.println("Trimmed Alignment [" + seq + "] :\t" + aln);
			}
			
			// Requiring end to end alignment here...
			if (aln.read_begin1 == 0 && aln.read_end1 == seq.length()-1) {
				
				// Pad with remaining reference sequence
				String leftPad = ref.substring(0, aln.ref_begin1);
				String rightPad = ref.substring(aln.ref_end1+1,ref.length());
				String paddedSeq = leftPad + seq + rightPad;
				String cigar = CigarUtils.extendCigarWithMatches(aln.cigar, leftPad.length(), rightPad.length());
				
				System.err.println("padded seq: " + paddedSeq);
				
				if (junctionPos > 0) {
					String oldCigar = cigar;
					cigar = CigarUtils.injectSplice(cigar, junctionPos, junctionLength);
					System.err.println("Spliced Cigar.  old: " + oldCigar + ", new: " + cigar);
				}
				
				result = new SSWAlignerResult(aln.ref_begin1-leftPad.length(), cigar, refChr, refStart, paddedSeq, aln.score1);
			}
		}
		
		return result;
	}
	
	public static class SSWAlignerResult {
		
		// Used for testing
		static boolean PAD_CONTIG = true;
		
		private int localRefPos;
		private String cigar;
		
		private String chromosome;
		private int refContextStart;
		
		private String sequence;
		private short score;
		
		SSWAlignerResult(int refPos, String cigar, String chromosome, int refContextStart, String sequence, short score) {
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
		
		public short getScore() {
			return score;
		}
	}
	
	public static void main(String[] args) {
		String ref = "AACAACAGATAATAACAAGTCCTAACCCTCTAGCTGCTTAGGCTGGCGGAGGCCCAGGGGCTCCCACGAGTTGGGTCCTTTCGCACCAGCACAGACTTACCTGATCTCGGTTGTTGATGTGAGAATAAGGAAGCTCCCCCGTCATCAGTTCATACAATACGATGCCATAGGAGTAGACATCCGACTGGAAACTGAATGGGTTGTTATCCTGCATTCGGATCACCTCTGGGGCCTACATGTATCACCATATGACAAAAGTGCATTTATCACCATATGACAGGCCTCACAGACATCTAGGGGCCAGGCTGTCCCTTTCATTAGTTATGAATGAG";
		String contig = "AACAACAGATAATAACAAGTCCTAACCCTCTAGCTGCTTAGGCTGGCGGAGGCCCAGGGGCTCCCACGAGTTGGGTCCTTTCGCACCAGCACAGACTTACCTGATCTCGGTTGTTGATGTGAGAATAAGGAAGCTCCCCCGTCATCAGTTCACAAAAGTGCATTTATCACCATATGACAGGCCTCACAGACATCTAGGGGCCAGGCTGTCCCTTTCATTAGTTATGAAT";
		
		SSWAligner sw = new SSWAligner(ref, "chr3", 12626521);
		sw.align(contig);
	}
}

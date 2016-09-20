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
	private static final int GAP_EXTEND_PENALTY = 1;
	
	private String refChr;
	private int refStart;
	private String ref;
	
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
		this.ref = ref;
		this.refChr = refChr;
		this.refStart = refStart;
	}
	
	public SSWAlignerResult align(String seq) {
		SSWAlignerResult result = null;
		
		Alignment aln = Aligner.align(seq.getBytes(), ref.getBytes(), score, GAP_OPEN_PENALTY, GAP_EXTEND_PENALTY, true);
		
		System.out.println("Alignment [" +  seq + "] :\t" + aln);
		
		// TODO: Optimize score requirements..
		if (aln != null && aln.score1 >= MIN_ALIGNMENT_SCORE && aln.score1 > aln.score2) {
			
			// Requiring end to end alignment here...
			if (aln.read_begin1 == 0 && aln.read_end1 == seq.length()-1) {
				result = new SSWAlignerResult(aln.ref_begin1, aln.cigar, refChr, refStart);
			}
		}
		
		return result;
	}
	
	public static class SSWAlignerResult {
		private int localRefPos;
		private String cigar;
		
		private String chromosome;
		private int refContextStart;
		
		SSWAlignerResult(int refPos, String cigar, String chromosome, int refContextStart) {
			this.localRefPos = refPos;
			this.cigar = cigar;
			this.chromosome = chromosome;
			this.refContextStart = refContextStart;
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
	}
	
	public static void main(String[] args) {
		String ref = "AACAACAGATAATAACAAGTCCTAACCCTCTAGCTGCTTAGGCTGGCGGAGGCCCAGGGGCTCCCACGAGTTGGGTCCTTTCGCACCAGCACAGACTTACCTGATCTCGGTTGTTGATGTGAGAATAAGGAAGCTCCCCCGTCATCAGTTCATACAATACGATGCCATAGGAGTAGACATCCGACTGGAAACTGAATGGGTTGTTATCCTGCATTCGGATCACCTCTGGGGCCTACATGTATCACCATATGACAAAAGTGCATTTATCACCATATGACAGGCCTCACAGACATCTAGGGGCCAGGCTGTCCCTTTCATTAGTTATGAATGAG";
		String contig = "AACAACAGATAATAACAAGTCCTAACCCTCTAGCTGCTTAGGCTGGCGGAGGCCCAGGGGCTCCCACGAGTTGGGTCCTTTCGCACCAGCACAGACTTACCTGATCTCGGTTGTTGATGTGAGAATAAGGAAGCTCCCCCGTCATCAGTTCACAAAAGTGCATTTATCACCATATGACAGGCCTCACAGACATCTAGGGGCCAGGCTGTCCCTTTCATTAGTTATGAAT";
		
		SSWAligner sw = new SSWAligner(ref, "chr3", 12626521);
		sw.align(contig);
	}
}

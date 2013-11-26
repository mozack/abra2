package abra.utils.trio;

public class LocusGenotype  {

	private String chromosome;
	private int pos;
	private int end;
	private String ref;
	private String alt1;
	private String alt2;
//	private boolean isPass;
	private double altAlleleFreq;
	private double qual;
	
	public enum Genotype { REF_REF, REF_ALT1, ALT1_ALT1, ALT1_ALT2, UNK };
	
	public Genotype gt = Genotype.REF_REF;
	
	public static final LocusGenotype REFERENCE = new LocusGenotype(TrioVcfReader.NO_CHROMOSOME, Integer.MAX_VALUE, Integer.MAX_VALUE, null, null, null, Genotype.REF_REF, 1000, 0);
	
//	private String filter;
	
	public LocusGenotype(String chromosome, int pos, int end, String ref, String alt1, String alt2, Genotype gt, double qual, double altAlleleFreq) {
		this.chromosome = chromosome;
		this.pos = pos;
		this.end = end;
		this.ref = ref;
		this.alt1 = alt1;
		this.alt2 = alt2;
		this.gt = gt;
		this.altAlleleFreq = altAlleleFreq;
		this.qual = qual;
		
//		this.isPass = filter.equals("PASS") ? true : false;
	}

	public String getRef() {
		return ref;
	}

	public String getAlt1() {
		return alt1;
	}

	public String getAlt2() {
		return alt2;
	}

	public Genotype getGt() {
		return gt;
	}

	public String getChromosome() {
		return chromosome;
	}

	public int getPos() {
		return pos;
	}
	
	public int getEnd() {
		return end;
	}
	
	public boolean isSnp() {
		boolean isSnp = false;
		
		if ((gt == Genotype.REF_ALT1) || (gt == Genotype.ALT1_ALT1)) {
			if (ref.length() == 1 && alt1.length() == 1) {
				isSnp = true;
			}
		}
		
		if (gt == Genotype.ALT1_ALT2) {
			if (alt1.length() == 1 && alt2.length() == 1) {
				isSnp = true;
			}
		}
		
		if (gt == Genotype.UNK) {
			if (ref.length() == 1 && alt1.length() == 1) {
				isSnp = true;
			}
		}
		
		return isSnp;
	}
	
	public boolean isSnpOrReference() {
		return gt == Genotype.REF_REF || isSnp();
	}
	
	public boolean isIndel() {
		return !isSnpOrReference();
	}
	
	public int getLenth() {
		int length = 1;
		
		if (isIndel()) {
			length = Math.max(alt1.length()-1, ref.length()-1);
		}
		
		return length;
	}
	
	private String genotypeString() {
		switch (gt) {
			case REF_REF:
				return "0/0";
			case REF_ALT1:
				return "0/1";
			case ALT1_ALT1:
				return "1/1";
			case ALT1_ALT2:
				return "1/2";
			default:
				return "-/-";
		}
	}
	
	public double getAltAlleleFreq() {
		return this.altAlleleFreq;
	}
	
	public String summary() {
		if (this == REFERENCE) {
			return "REFERENCE";
		}
		
		return chromosome + ":" + pos + "-" + end + ":" + genotypeString() + ":" + ref + ":" + alt1 + ":" + (alt2 != null ? alt2 : ""); 
				
	}
	
	public String toString() {
		return summary();
	}
	
	public boolean isVariant() {
//		return gt == Genotype.REF_ALT1 || gt == Genotype.ALT1_ALT1 || gt == Genotype.ALT1_ALT2;
		return gt == Genotype.REF_ALT1 || gt == Genotype.ALT1_ALT1;
	}
		
	public boolean isAbraIndelConflict() {
		if (this == REFERENCE) return false;
		
		// Due to multiple indels at single locus (not handled well upstream)
		if ((ref.length() > 1) && (alt2 != null)) {
			return true;
		}
		
		// Due to multiple indels at single locus (not handled well upstream)
		if ((alt2 != null) && (alt1.length() > 1) && (alt2.length() > 1)) {
			return true;
		}
		
		// filter complex calls (for now...)
		if ((ref.length() > 1) && (alt1.length() > 1)) {
			return true;
		}
	
		return false;
	} 
	
	public boolean isFiltered() {
		
		/*
		return qual < 10 ||              // Did not pass caller's filters
//		return qual < 150.0 ||              // Did not pass caller's filters
				gt == Genotype.UNK;    // Input gt not in N/N format
		*/

		return gt == Genotype.UNK;
	}
}

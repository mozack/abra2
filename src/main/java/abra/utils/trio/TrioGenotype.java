package abra.utils.trio;

import java.util.ArrayList;
import java.util.List;

import abra.utils.trio.LocusGenotype.Genotype;

public class TrioGenotype {

	private LocusGenotype father;
	private LocusGenotype mother;
	private LocusGenotype child;
	
	private int otherEvents = 0;
	
	public enum GenotypeStatus { 
		VALID("VALID"), CONFLICT("CONFLICT"), FILTERED("FILTERED"), INDEL_CONFLICT_ISAAC("FILTERED_INDEL_CONFLICT_ISAAC"), INDEL_CONFLICT_ABRA("FILTERED_INDEL_CFLICT_ABRA"), INDEL_LOW_FREQ("FILTERED_INDEL_LOW_FREQ");
	
		private String str;
		private GenotypeStatus(String str) {
			this.str = str;
		}
		
		public String toString() {
			return str;
		}
	};
	
	
	public TrioGenotype(LocusGenotype father, LocusGenotype mother, LocusGenotype child) {
		
		if (father == LocusGenotype.REFERENCE && mother == LocusGenotype.REFERENCE && child == LocusGenotype.REFERENCE) {
			throw new IllegalArgumentException("Entire trio should not be reference at single locus");
		}
		
		this.father = father;
		this.mother = mother;
		this.child = child;
	}

	public LocusGenotype getFather() {
		return father;
	}

	public LocusGenotype getMother() {
		return mother;
	}

	public LocusGenotype getChild() {
		return child;
	}
	
	public String summary() {
		String varType = "REF";
		if (hasVariant()) {
			if (isSnpOrReference()) {
				varType = "SNP";
			} else {
				varType = "INDEL";
			}
		}
		
//		String chr = "";
//		if (father != LocusGenotype.REFERENCE) {
//			chr = father.getChromosome();
//		} else if (mother != LocusGenotype.REFERENCE) {
//			chr = mother.getChromosome();
//		} else if (child != LocusGenotype.REFERENCE) {
//			chr = child.getChromosome();
//		}
//		
//		int startPos = 0;
		
		String range = getChromosome() + "\t" + min(father.getPos(), mother.getPos(), child.getPos()) + "\t" +
				max(father.getEnd(), mother.getEnd(), child.getEnd());
		
		return range + "\t" + getStatus() + "\tf:" + father.summary() + "\tm:" + mother.summary() + "\tc:" + child.summary() + "\t" + varType + "\t" + getLength();  
	}
	
	private String getChromosome() {
		String chromosome = "";
		
		if (father != LocusGenotype.REFERENCE) {
			chromosome = father.getChromosome();
		} else if (mother != LocusGenotype.REFERENCE) {
			chromosome = mother.getChromosome();
		} else if (child != LocusGenotype.REFERENCE) {
			chromosome = child.getChromosome();
		}
		
		return chromosome;
	}
	
	private int min(int i1, int i2, int i3) {
		int min = Math.min(i1, i2);
		min = Math.min(min, i3);
		
		return min;
	}
	
//	private int maxEnd(LocusGenotype lg1, LocusGenotype lg2, LocusGenotype lg3) {
//		int max = -1;
//		if (lg1 != LocusGenotype.REFERENCE) {
//			max = lg1.getEnd();
//		}
//		
//		if (lg2 != LocusGenotype.REFERENCE) {
//			max = Math.max(max,  lg2.getEnd());
//		}
//		
//		if (lg3 != LocusGenotype.REFERENCE) {
//			max = Math.max(max,  lg3.getEnd());
//		}
//	}
	
	private int max(int i1, int i2, int i3) {
		int max = Math.max(i1, i2);
		max = Math.max(max, i3);
		
		return max;
	}
	
	public String toString() {
		return summary();
	}
	
	public boolean hasVariant() {
		return father.isVariant() || mother.isVariant() || child.isVariant();
	}
	
	public boolean isSnpOrReference() {
		//return hasVariant() && father.isSnpOrReference() && mother.isSnpOrReference() && child.isSnpOrReference();
		return father.isSnpOrReference() && mother.isSnpOrReference() && child.isSnpOrReference();
	}
	
	public int getLength() {
		return max(father.getLenth(), child.getLenth(), mother.getLenth());
	}
	
	public GenotypeStatus getStatus() {

		if (father.isFiltered() || mother.isFiltered() || child.isFiltered()) {
//			if (father.isFiltered() && mother.isFiltered() && child.isFiltered()) {
				return GenotypeStatus.FILTERED;
			}
		
		if (father.isAbraIndelConflict() || mother.isAbraIndelConflict() || child.isAbraIndelConflict()) {
			return GenotypeStatus.INDEL_CONFLICT_ABRA;
		}
		
		GenotypeStatus status = GenotypeStatus.CONFLICT;
		
		String[] childAlleles = getAlleles(child);
		String[] motherAlleles = getAlleles(mother);
		String[] fatherAlleles = getAlleles(father);
		
		String c1 = childAlleles[0];
		String c2 = childAlleles[1];
		String m1 = motherAlleles[0];
		String m2 = motherAlleles[1];
		String f1 = fatherAlleles[0];
		String f2 = fatherAlleles[1];
		
		if ((c1.equals(m1) || c1.equals(m2)) &&
			(c2.equals(f1) || c2.equals(f2))) {
			
			status = GenotypeStatus.VALID;
		}
		
		if ((c2.equals(m1) || c2.equals(m2)) &&
			(c1.equals(f1) || c1.equals(f2))) {
				
			status = GenotypeStatus.VALID;
		}
		
		boolean isThresholdExceeded = false;
		double MIN_ALT_ALLELE_FREQ = 0.2;
//		if ((status == GenotypeStatus.CONFLICT) && !(isSnpOrReference())) {
		
		if (!(isSnpOrReference())) {
			if (child.isIndel() && child.getAltAlleleFreq() >  MIN_ALT_ALLELE_FREQ) {
				isThresholdExceeded = true;
			} else if (mother.isIndel() && mother.getAltAlleleFreq() >  MIN_ALT_ALLELE_FREQ) {
				isThresholdExceeded = true;
			} else if (father.isIndel() && father.getAltAlleleFreq() >  MIN_ALT_ALLELE_FREQ) {
				isThresholdExceeded = true;
			}
			
			// Filter indels @ < 20% as conflicts or valid
//			if ((child.isIndel() || mother.isIndel() || father.isIndel()) && !isThresholdExceeded) {
//				status = GenotypeStatus.INDEL_LOW_FREQ;
//			}
		}
				
		
		return status;
	}
	
	private String[] getAlleles(LocusGenotype lg) {
		if (lg.getGt() == Genotype.REF_REF) {
			return new String[] { "REF", "REF" };
		} else if (lg.getGt() == Genotype.REF_ALT1) {
//			if (lg.isIndel() && lg.getAltAlleleFreq() < .2) {
//				return new String[] { "REF", "REF" };
//			}
			return new String[] { "REF", lg.getRef() + ":" + lg.getAlt1() };
		} else if (lg.getGt() == Genotype.ALT1_ALT1) {
			String alt1 = lg.getRef() + ":" + lg.getAlt1();
			return new String[] { alt1, alt1 };
		} else if (lg.getGt() == Genotype.ALT1_ALT2) {
			
			// TODO: Need factor this down to simplest representation for comparison
			String alt1 = lg.getRef() + ":" + lg.getAlt1();
			String alt2 = lg.getRef() + ":" + lg.getAlt2();
			
			return new String[] { alt1, alt2 };
		} else {
			return null;
		}
	}
}

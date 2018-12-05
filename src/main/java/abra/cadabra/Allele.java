package abra.cadabra;

public class Allele {
	public enum Type {
		A, T, C, G, INS, DEL, MNP, UNK
	}

	private Type type;
	private int length;
	
	private String seq;
	
	public static final Allele UNK = new Allele(Type.UNK, 1);
	public static final Allele A = new Allele(Type.A, 1);
	public static final Allele T = new Allele(Type.T, 1);
	public static final Allele C = new Allele(Type.C, 1);
	public static final Allele G = new Allele(Type.G, 1);
	
	public Allele(Type type, int length) {
		this.type = type;
		this.length = length;
	}
	
	public static Allele getMnpAllele(String seq) {
		Allele allele = new Allele(Type.MNP, seq.length());
		allele.seq = seq;
		
		return allele;
	}
		
	public static Allele getAllele (char base) {
		switch (base) {
			case 'A':
				return A;
			case 'T':
				return T;
			case 'C':
				return C;
			case 'G':
				return G;
			default:
				return UNK;
		}
	}
	
	public Type getType() {
		return type;
	}
	
	public int getLength() {
		return length;
	}
	
	public String toString() {
		String str = "";
		switch (type) {
			case A:
				str = "A";
				break;
			case C:
				str = "C";
				break;				
			case T:
				str = "T";
				break;
			case G:
				str = "G";
				break;
			case INS:
				str = "INS";
				break;
			case DEL:
				str = "DEL";
				break;
			case MNP:
				str = "MNP";
				break;
			case UNK:
				str = "UNK";
				break;
		}
		
		return str;
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + length;
		result = prime * result + ((seq == null) ? 0 : seq.hashCode());
		result = prime * result + ((type == null) ? 0 : type.hashCode());
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
		Allele other = (Allele) obj;
		if (length != other.length)
			return false;
		if (seq == null) {
			if (other.seq != null)
				return false;
		} else if (!seq.equals(other.seq))
			return false;
		if (type != other.type)
			return false;
		return true;
	}
}

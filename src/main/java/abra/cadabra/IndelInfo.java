package abra.cadabra;

import htsjdk.samtools.CigarElement;

public class IndelInfo {

	private CigarElement cigarElement;
	private int readIndex = -1;
	private String insertBases;

	public IndelInfo(CigarElement cigarElement, int readIndex) {
		this.cigarElement = cigarElement;
		this.readIndex = readIndex;
	}

	public CigarElement getCigarElement() {
		return cigarElement;
	}

	public int getReadIndex() {
		return readIndex;
	}
	
	public void setReadIndex(int readIndex) {
		this.readIndex = readIndex;
	}
	
	public String getInsertBases() {
		return insertBases;
	}

	public void setInsertBases(String insertBases) {
		this.insertBases = insertBases;
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result
				+ ((cigarElement == null) ? 0 : cigarElement.hashCode());
		result = prime * result
				+ ((insertBases == null) ? 0 : insertBases.hashCode());
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
		IndelInfo other = (IndelInfo) obj;
		if (cigarElement == null) {
			if (other.cigarElement != null)
				return false;
		} else if (!cigarElement.equals(other.cigarElement))
			return false;
		if (insertBases == null) {
			if (other.insertBases != null)
				return false;
		} else if (!insertBases.equals(other.insertBases))
			return false;
		return true;
	}
}

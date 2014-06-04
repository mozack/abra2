package abra.cadabra;

import net.sf.samtools.CigarElement;

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
}

package abra;

/**
 * Representation of a Feature (i.e. a line in a GTF file)
 * 
 * @author Lisle Mose (lmose at unc dot edu)
 */
public class Feature {
	private String seqname;
	private long start;  // 1 based
	private long end;    // inclusive
	
	public Feature(String seqname, long start, long end) {
		this.seqname      = seqname;
		this.start        = start;
		this.end          = end;
	}
	
	public String getSeqname() {
		return seqname;
	}
	
	public long getStart() {
		return start;
	}
	
	public long getEnd() {
		return end;
	}
	
	public String getDescriptor() {
		return seqname + "_" + start + "_" + end;
	}
	
	public long getLength() {
		return end-start;
	}
	
	public String toString() {
		return getDescriptor();
	}
}

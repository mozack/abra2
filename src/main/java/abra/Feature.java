/* Copyright 2013 University of North Carolina at Chapel Hill.  All rights reserved. */
package abra;

import net.sf.samtools.SAMRecord;

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
	
	public void pad(long padding) {
		start -= padding;
		end += padding;
	}
	
	public void setEnd(long end) {
		this.end = end;
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
	
	private boolean spansCoordinate(int coord) {
		return (coord >= start) && (coord <= end); 
	}
	
	private boolean spansEitherCoordinate(int coord1, int coord2) {
		return spansCoordinate(coord1) || spansCoordinate(coord2);
	}
	
	public boolean overlapsRead(SAMRecord read) {
		return ((this.seqname.equals(read.getReferenceName())) &&
			(spansEitherCoordinate(read.getAlignmentStart(), read.getAlignmentEnd())));
	}
	
	public boolean overlaps(String chromosome, int startPos, int stopPos) {
		return ((this.seqname.equals(chromosome)) &&
				(spansEitherCoordinate(startPos, stopPos)));
	}
}

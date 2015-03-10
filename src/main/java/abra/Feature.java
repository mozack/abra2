/* Copyright 2013 University of North Carolina at Chapel Hill.  All rights reserved. */
package abra;

import htsjdk.samtools.SAMRecord;

/**
 * Representation of a Feature (i.e. a line in a GTF file)
 * 
 * @author Lisle E. Mose (lmose at unc dot edu)
 */
public class Feature {
	private String seqname;
	private long start;  // 1 based
	private long end;    // inclusive
	private String additionalInfo;
	
	// Optional kmerSize value specific to ABRA assembly.
	private int kmerSize;

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
	
	public void setAdditionalInfo(String info) {
		this.additionalInfo = info;
	}
	
	public String getAdditionalInfo() {
		return additionalInfo;
	}
	
	
	private boolean isWithin(long coord, long start, long stop) {
		return coord >= start && coord <= stop;
	}
	
	private boolean overlaps(long start1, long stop1, long start2, long stop2) {
		return 
				isWithin(start1, start2, stop2) ||
				isWithin(stop1, start2, stop2) ||
				isWithin(start2, start1, stop1) ||
				isWithin(stop2, start1, stop1);
	}
	
	public boolean overlapsRead(SAMRecord read) {
		return overlaps(read.getReferenceName(), read.getAlignmentStart(), read.getAlignmentEnd());
	}
	
	public boolean overlaps(String chromosome, int startPos, int stopPos) {
		return ((this.seqname.equals(chromosome)) && overlaps(start, end, startPos, stopPos));
	}
	
	public int getKmer() {
		return kmerSize;
	}

	public void setKmer(int kmer) {
		this.kmerSize = kmer;
	}
}

/* Copyright 2013 University of North Carolina at Chapel Hill.  All rights reserved. */
package abra;

import java.util.List;

import htsjdk.samtools.SAMFileHeader;
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
		int alignmentEnd = Math.max(read.getAlignmentEnd(), read.getAlignmentStart() + read.getReadLength());
		return overlaps(read.getReferenceName(), read.getAlignmentStart(), alignmentEnd);
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
	
	public static int findFirstOverlappingRegion(SAMFileHeader samHeader, SAMRecord read, List<Feature> regions, int start) {
		if (start < 0) {
			start = 0;
		}
		
		// Adjust alignment end for soft clipping / unmapped reads
		int alignmentEnd = Math.max(read.getAlignmentEnd(), read.getAlignmentStart() + read.getReadLength());

		for (int idx=start; idx<regions.size(); idx++) {
			Feature region = regions.get(idx);
			if ( (read.getReferenceIndex() < samHeader.getSequenceDictionary().getSequenceIndex(region.getSeqname())) ||
				 (read.getReferenceName().equals(region.getSeqname()) && alignmentEnd < region.getStart()) ) {
				
				// This read is in between regions
				// TODO: adjust start region here
				return -1;
			} else if (region.overlaps(read.getReferenceName(), read.getAlignmentStart(), alignmentEnd)) {
				return idx;
			}
		}
		
		// This read is beyond all regions
		return -1;
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + (int) (end ^ (end >>> 32));
		result = prime * result + ((seqname == null) ? 0 : seqname.hashCode());
		result = prime * result + (int) (start ^ (start >>> 32));
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
		Feature other = (Feature) obj;
		if (end != other.end)
			return false;
		if (seqname == null) {
			if (other.seqname != null)
				return false;
		} else if (!seqname.equals(other.seqname))
			return false;
		if (start != other.start)
			return false;
		return true;
	}
}

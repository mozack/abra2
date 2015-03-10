/* Copyright 2013 University of North Carolina at Chapel Hill.  All rights reserved. */
package abra;

import htsjdk.samtools.SAMRecord;

/**
 * A Pair of SAMRecord's
 * 
 * @author Lisle E. Mose (lmose at unc dot edu)
 */
public class ReadPair {
    private SAMRecord read1;
    private SAMRecord read2;
    
    private String hashString = null;
    
    ReadPair(SAMRecord read1, SAMRecord read2) {
        this.read1 = read1;
        this.read2 = read2;
    }
    
    public SAMRecord getRead1() {
        return read1;
    }
    
    public SAMRecord getRead2() {
        return read2;
    }
    
    public String toString() {
        String r1 = read1 != null ? read1.getReadName() : "null";
        String r2 = read2 != null ? read2.getReadName() : "null";
        return "read1: " + r1 + ", read2: " + r2;
    }
    
    private synchronized String getHashString() {
    	if (hashString == null) {
    		hashString = read1 != null ? read1.getSAMString() : "null" +
    					 read2 != null ? read2.getSAMString() : "null";
    	}
    	
    	return hashString;
    }
    
    @Override
    public int hashCode() {
    	return getHashString().hashCode();
    }
    
    @Override
    public boolean equals(Object obj) {
    	ReadPair that = (ReadPair) obj;
    	return this.getHashString().equals(that.getHashString());
    }
}

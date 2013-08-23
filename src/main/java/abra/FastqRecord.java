/* Copyright 2013 University of North Carolina at Chapel Hill.  All rights reserved. */
package abra;

import java.util.Arrays;


/**
 * Representation of a single Fastq record.
 * Some code here may be specific to paired end processing.
 * 
 * @author Lisle Mose (lmose at unc dot edu)
 */
public class FastqRecord {
    
    public static final int NUM_LINES = 4;
    
    private String[] lines = new String[NUM_LINES];
    
    private QualityConverter qualityConverter;
    
    public FastqRecord(String[] lines) {
        if (lines.length != NUM_LINES) {
            throw new IllegalArgumentException("Invalid number of lines for FastqRecord: [" + Arrays.toString(lines) + "]");
        }
        this.lines = lines;
    }
    
    FastqRecord(String id, String sequence, String quality) {
    	lines[0] = id;
    	lines[1] = sequence;
    	lines[2] = "+";
    	lines[3] = quality;
    }
    
    public String getId() {
        return lines[0];
    }
    
    public String[] getLines() {
        return lines;
    }
    
    public void setQualityConverter(QualityConverter qualityConverter) {
    	this.qualityConverter = qualityConverter;
    }
    
    /**
     * Returns the portion of the id string leading up to "/"
     */
    public String getBaseId() {
        int slashIdx = getId().indexOf("/");
        int spaceIdx = getId().indexOf(" ");
        
        if ((slashIdx == -1) && (spaceIdx == -1)) {
            return getId();
        }
        
        int idx = -1;
        if (slashIdx == -1) {
            idx = spaceIdx;
        } else if (spaceIdx == -1) {
            idx = slashIdx;
        } else {
            idx = spaceIdx < slashIdx ? spaceIdx : slashIdx;
        }
        
        return getId().substring(0, idx);
    }
    
    /**
     * Returns true if this FastqRecord has the same base id as the input FastqRecord 
     */
    public boolean hasSameBaseId(FastqRecord rec) {
        return rec != null && this.getBaseId().equals(rec.getBaseId());
    }
    
    public String toString() {
        return Arrays.toString(lines);
    }
    
    public int hashcode() {
        return Arrays.hashCode(lines);
    }
    
    public boolean equals(Object obj) {
        FastqRecord that = (FastqRecord) obj;
        return Arrays.equals(this.lines, that.lines);
    }
    
    public void stripNonReadInfoInId() {
        int idx = lines[0].indexOf(" ");
        if (idx > 0) {
            lines[0] = lines[0].substring(0, idx);
        }
    }
    
    public void appendToId(String suffix) {
        if (!lines[0].endsWith(suffix)) {
            lines[0] = lines[0] + suffix;
        }
    }
    
    String getQuality() {
    	return lines[3];
    }
    
    public String getSequence() {
    	return lines[1];
    }
    
    private void setQuality(String quality) {
    	lines[3] = quality;
    }
    
    public void phred33To64() {
    	
    	if (qualityConverter == null) {
    		throw new RuntimeException("Please set QualityConverter.");
    	}
    	
    	String phred33 = getQuality();
    	    	
    	setQuality(qualityConverter.phred33ToPhred64(phred33));
    }
}

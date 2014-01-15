/* Copyright 2013 University of North Carolina at Chapel Hill.  All rights reserved. */
package abra;

/**
 * Utility class that can be used to convert quality scores.
 * 
 * @author Lisle E. Mose (lmose at unc dot edu)
 */
public class QualityConverter {
	
	private static final int PHRED33_TO_PHRED64_DIFF = 31;

	/**
	 * Convert from phred33 to phred64
	 */
    public String phred33ToPhred64(String phred33) {
    	StringBuffer phred64 = new StringBuffer();
    	
    	for (int i=0; i<phred33.length(); i++) {
    		phred64.append((char) (phred33.charAt(i) + PHRED33_TO_PHRED64_DIFF));
    	}
    	
    	return phred64.toString();
    }

	/**
	 * Convert from phred64 to phred33
	 */
    public String phred64ToPhred33(String phred64) {
    	StringBuffer phred33 = new StringBuffer();
    	
    	for (int i=0; i<phred64.length(); i++) {
    		phred33.append((char) (phred64.charAt(i) - PHRED33_TO_PHRED64_DIFF));
    	}
    	
    	return phred33.toString();
    }    
}

/* Copyright 2013 University of North Carolina at Chapel Hill.  All rights reserved. */
package abra;

import htsjdk.samtools.SAMRecord;

/**
 * Concrete implementations of this interface are responsible for outputting
 * realigned reads to the output BAM file.
 * 
 * @author Lisle E. Mose (lmose at unc dot edu)
 */
public interface RealignmentWriter {

	public void addAlignment(SAMRecord updatedRead, SAMRecord origRead);
	
	public int flush();
}

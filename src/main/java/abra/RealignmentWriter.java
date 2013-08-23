/* Copyright 2013 University of North Carolina at Chapel Hill.  All rights reserved. */
package abra;

import net.sf.samtools.SAMRecord;

public interface RealignmentWriter {

	public void addAlignment(SAMRecord updatedRead, SAMRecord origRead);
	
	public int flush();
	
}

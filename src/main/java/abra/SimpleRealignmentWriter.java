/* Copyright 2013 University of North Carolina at Chapel Hill.  All rights reserved. */
package abra;

import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMRecord;

/**
 * Simple realignment writer.  Does not consider read pair information when writing to output BAM.
 * 
 * @author Lisle E. Mose (lmose at unc dot edu)
 */
public class SimpleRealignmentWriter implements RealignmentWriter {

	private SAMFileWriter writer;
	private int realignCount = 0;
	private IndelShifter indelShifter = new IndelShifter();
	private boolean isTightAlignment = false;
	private CompareToReference2 c2r;
	
	public SimpleRealignmentWriter(CompareToReference2 c2r, SAMFileWriter writer, boolean isTightAlignment) {
		this.writer = writer;
		this.c2r = c2r;
		this.isTightAlignment = isTightAlignment;
	}
	
	@Override
	public void addAlignment(SAMRecord updatedRead, SAMRecord origRead) {
		
		if (updatedRead != null) {
			// Output realigned read
			addAlignment(updatedRead);
			if (updatedRead.getAttribute("YO") != null) {
				realignCount += 1;
			}
		} else {
			// Output original read
//			realigner.adjustForStrand(contigAlignedRead.getReadNegativeStrandFlag(), origRead);
			addAlignment(origRead);
		}
	}
	
	private void addAlignment(SAMRecord read) {
		if (isTightAlignment) {
			writer.addAlignment(indelShifter.shiftIndelsLeft(read, c2r));
		} else {
			writer.addAlignment(read);
		}
	}

	@Override
	public int flush() {
		return realignCount;
	}

}

/* Copyright 2013 University of North Carolina at Chapel Hill.  All rights reserved. */
package abra;

import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMRecord;

public class SimpleRealignmentWriter implements RealignmentWriter {

	private SAMFileWriter writer;
	private ReAligner realigner;
	private int realignCount = 0;
	private IndelShifter indelShifter = new IndelShifter();
	private boolean isTightAlignment = false;
	
	public SimpleRealignmentWriter(ReAligner realigner, SAMFileWriter writer, boolean isTightAlignment) {
		this.writer = writer;
		this.realigner = realigner;
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
			writer.addAlignment(indelShifter.shiftIndelsLeft(read, realigner.getC2r()));
		} else {
			writer.addAlignment(read);
		}
	}

	@Override
	public int flush() {
		return realignCount;
	}

}

/* Copyright 2013 University of North Carolina at Chapel Hill.  All rights reserved. */
package abra;

import java.io.IOException;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileWriter;

/**
 * Thread runnable class for read adjustment.
 * 
 * @author Lisle E. Mose (lmose at unc dot edu)
 */
public class AdjustReadsRunnable implements Runnable {
	
	private ReadAdjuster readAdjuster;
	private String sortedAlignedToContig;
	private SAMFileWriter outputSam;
	private boolean isTightAlignment;
	private String tempDir;
	private SAMFileHeader samHeader;
	
	public AdjustReadsRunnable(ReadAdjuster readAdjuster, String sortedAlignedToContig, SAMFileWriter outputSam,
			boolean isTightAlignment, String tempDir, SAMFileHeader samHeader) {

		this.readAdjuster = readAdjuster;
		this.sortedAlignedToContig = sortedAlignedToContig;
		this.outputSam = outputSam;
		this.isTightAlignment = isTightAlignment;
		this.tempDir = tempDir;
		this.samHeader = samHeader;
	}

	@Override
	public void run() {
		try {
			readAdjuster.adjustReads(sortedAlignedToContig, outputSam, isTightAlignment, tempDir, samHeader);
		} catch (IOException e) {
			e.printStackTrace();
			throw new RuntimeException(e);
		}
	}
}

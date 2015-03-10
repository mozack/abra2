/* Copyright 2013 University of North Carolina at Chapel Hill.  All rights reserved. */
package abra;

import java.io.IOException;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;

/**
 * Thread runnable class for read adjustment.
 * 
 * @author Lisle E. Mose (lmose at unc dot edu)
 */
public class AdjustReadsRunnable extends AbraRunnable {
	
	private ReadAdjuster readAdjuster;
	private String sortedAlignedToContig;
	private SAMFileWriter outputSam;
	private boolean isTightAlignment;
	private String tempDir;
	private SAMFileHeader samHeader;
	
	public AdjustReadsRunnable(ThreadManager threadManager, ReadAdjuster readAdjuster, String sortedAlignedToContig, SAMFileWriter outputSam,
			boolean isTightAlignment, String tempDir, SAMFileHeader samHeader) {

		super(threadManager);
		this.readAdjuster = readAdjuster;
		this.sortedAlignedToContig = sortedAlignedToContig;
		this.outputSam = outputSam;
		this.isTightAlignment = isTightAlignment;
		this.tempDir = tempDir;
		this.samHeader = samHeader;
	}

	@Override
	public void go() throws Exception {
		readAdjuster.adjustReads(sortedAlignedToContig, outputSam, isTightAlignment, tempDir, samHeader);
	}
}

/* Copyright 2013 University of North Carolina at Chapel Hill.  All rights reserved. */
package abra;

import java.io.IOException;
import java.io.InputStream;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;

/**
 * Thread runnable class for read adjustment.
 * 
 * @author Lisle E. Mose (lmose at unc dot edu)
 */
public class AdjustReadsStreamRunnable extends AbraRunnable {
	
	private ReadAdjuster readAdjuster;
	private SAMFileWriter outputSam;
	private boolean isTightAlignment;
	private String tempDir;
	private SAMFileHeader samHeader;
	private InputStream inputStream;

	public AdjustReadsStreamRunnable(ThreadManager threadManager, ReadAdjuster readAdjuster, SAMFileWriter outputSam,
			boolean isTightAlignment, String tempDir, SAMFileHeader samHeader) {

		super(threadManager);
		this.readAdjuster = readAdjuster;
		this.outputSam = outputSam;
		this.isTightAlignment = isTightAlignment;
		this.tempDir = tempDir;
		this.samHeader = samHeader;
	}

	@Override
	public void go() throws Exception {
		readAdjuster.adjustReads(inputStream, outputSam, isTightAlignment, tempDir, samHeader);
	}
	
	public void setInputStream(InputStream inputStream) {
		this.inputStream = inputStream;
	}
}

/* Copyright 2013 University of North Carolina at Chapel Hill.  All rights reserved. */
package abra;

import java.io.IOException;
import java.io.InputStream;
import java.util.Queue;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;

/**
 * Thread runnable class for read adjustment.
 * 
 * @author Lisle E. Mose (lmose at unc dot edu)
 */
public class AdjustReadsQueueRunnable extends AbraRunnable {
	
	private ReadAdjuster readAdjuster;
	private SAMFileWriter outputSam;
	private boolean isTightAlignment;
	private String tempDir;
	private SAMFileHeader samHeader;
	private Queue<SAMRecord> queue;
	private MutableBoolean isDone;

	public AdjustReadsQueueRunnable(ThreadManager threadManager, ReadAdjuster readAdjuster, SAMFileWriter outputSam,
			boolean isTightAlignment, String tempDir, SAMFileHeader samHeader, MutableBoolean isDone) {

		super(threadManager);
		this.readAdjuster = readAdjuster;
		this.outputSam = outputSam;
		this.isTightAlignment = isTightAlignment;
		this.tempDir = tempDir;
		this.samHeader = samHeader;
		this.isDone = isDone;
	}

	@Override
	public void go() throws Exception {
		readAdjuster.adjustReads(queue, outputSam, isTightAlignment, tempDir, samHeader, isDone);
	}
	
	public void setReadQueue(Queue<SAMRecord> queue) {
		this.queue = queue;
	}
	
	public void setDone() {
		isDone.setValue(true);
	}
}

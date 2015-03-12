/* Copyright 2013 University of North Carolina at Chapel Hill.  All rights reserved. */
package abra;

import java.io.IOException;

import htsjdk.samtools.SAMFileWriter;

/**
 * Thread entry point for read alignment.
 * 
 * @author Lisle E. Mose (lmose at unc dot edu)
 */
@Deprecated
public class AlignReadsRunnable extends AbraRunnable {
		
	private String tempDir;
	private String inputSam;
	private String cleanContigsFasta;
	private CompareToReference2 c2r;
	private SAMFileWriter finalOutputSam;
	private String alignedToContigSam;
	private ReAligner reAligner;
	
	public AlignReadsRunnable(ThreadManager threadManager, ReAligner realigner, String tempDir, String inputSam, String cleanContigsFasta,
			CompareToReference2 c2r, SAMFileWriter finalOutputSam, String alignedToContigSam) {

		super(threadManager);
		this.reAligner = realigner;
		this.tempDir = tempDir;
		this.inputSam = inputSam;
		this.cleanContigsFasta = cleanContigsFasta;
		this.c2r = c2r;
		this.finalOutputSam = finalOutputSam;
		this.alignedToContigSam = alignedToContigSam;
	}

	@Override
	public void go() throws Exception{
//		reAligner.alignReads(tempDir, inputSam, cleanContigsFasta, c2r, finalOutputSam, alignedToContigSam);
	}
}

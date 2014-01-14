/* Copyright 2013 University of North Carolina at Chapel Hill.  All rights reserved. */
package abra;

import java.io.IOException;

import net.sf.samtools.SAMFileWriter;

/**
 * Thread entry point for read alignment.
 * 
 * @author Lisle E. Mose (lmose at unc dot edu)
 */
public class AlignReadsRunnable implements Runnable {
		
	private ReAligner realigner;
	private String tempDir;
	private String inputSam;
	private String cleanContigsFasta;
	private CompareToReference2 c2r;
	private SAMFileWriter finalOutputSam;
	private String alignedToContigSam;
	
	public AlignReadsRunnable(ReAligner realigner, String tempDir, String inputSam, String cleanContigsFasta,
			CompareToReference2 c2r, SAMFileWriter finalOutputSam, String alignedToContigSam) {

		this.realigner = realigner;
		this.tempDir = tempDir;
		this.inputSam = inputSam;
		this.cleanContigsFasta = cleanContigsFasta;
		this.c2r = c2r;
		this.finalOutputSam = finalOutputSam;
		this.alignedToContigSam = alignedToContigSam;
	}

	@Override
	public void run() {
		try {
			realigner.alignReads(tempDir, inputSam, cleanContigsFasta, c2r, finalOutputSam, alignedToContigSam);
		} catch (IOException e) {
			e.printStackTrace();
			throw new RuntimeException(e);
		} catch (InterruptedException e) {
			e.printStackTrace();
			throw new RuntimeException(e);
		}
	}
}

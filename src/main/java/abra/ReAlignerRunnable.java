/* Copyright 2013 University of North Carolina at Chapel Hill.  All rights reserved. */
package abra;

/**
 * Thread entry point for chromsome / reference sequence specific processing.
 * 
 * @author Lisle E. Mose (lmose at unc dot edu)
 */
public class ReAlignerRunnable extends AbraRunnable {
	private int chromosomeChunkIdx;
	private ReAligner reAligner;
	
	public ReAlignerRunnable(ThreadManager threadManager, ReAligner reAligner, int chromosomeChunkIdx) {
		super(threadManager);
		this.chromosomeChunkIdx = chromosomeChunkIdx;
		this.reAligner = reAligner;
	}
	
	@Override
	public void go() throws Exception {		
		reAligner.processChromosomeChunk(chromosomeChunkIdx);
	}
}

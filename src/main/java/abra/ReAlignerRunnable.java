/* Copyright 2013 University of North Carolina at Chapel Hill.  All rights reserved. */
package abra;

/**
 * Thread entry point for chromsome / reference sequence specific processing.
 * 
 * @author Lisle E. Mose (lmose at unc dot edu)
 */
public class ReAlignerRunnable extends AbraRunnable {
	private String chromosome;
	private ReAligner reAligner;
	
	public ReAlignerRunnable(ThreadManager threadManager, ReAligner reAligner, String chromosome) {
		super(threadManager);
		this.chromosome = chromosome;
		this.reAligner = reAligner;
	}
	
	@Override
	public void go() throws Exception {		
		reAligner.processChromosome(chromosome);
	}
}

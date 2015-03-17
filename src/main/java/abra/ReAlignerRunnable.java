/* Copyright 2013 University of North Carolina at Chapel Hill.  All rights reserved. */
package abra;

/**
 * Thread entry point for region specific processing.
 * 
 * @author Lisle E. Mose (lmose at unc dot edu)
 */
public class ReAlignerRunnable extends AbraRunnable {
	private Feature region;
	private ReAligner reAligner;
	
	public ReAlignerRunnable(ThreadManager threadManager, ReAligner reAligner, Feature region) {
		super(threadManager);
		this.region = region;
		this.reAligner = reAligner;
	}
	
	@Override
	public void go() throws Exception {		
		reAligner.processRegion(region);
	}
}

/* Copyright 2013 University of North Carolina at Chapel Hill.  All rights reserved. */
package abra;

/**
 * Thread entry point for region specific processing.
 * 
 * @author Lisle E. Mose (lmose at unc dot edu)
 */
public class ReAlignerRunnable implements Runnable {

	private ReAligner reAligner;
	private Feature region;
	private boolean isAlive = true;
	
	public ReAlignerRunnable(ReAligner reAligner, Feature region) {
		this.reAligner = reAligner;
		this.region = region;
	}
	
	@Override
	public void run() {
		
		try {
			reAligner.processRegion(region);
		} catch (Exception e) {
			e.printStackTrace();
			throw new RuntimeException(e);
		} finally {
			reAligner.removeThread(this);
			isAlive = false;
		}
	}
	
	public boolean isAlive() {
		return isAlive;
	}
}

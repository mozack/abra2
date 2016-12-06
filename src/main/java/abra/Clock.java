/* Copyright 2013 University of North Carolina at Chapel Hill.  All rights reserved. */
package abra;

/**
 * Simple class used to log elapsed wall clock times.
 * 
 * @author Lisle E. Mose (lmose at unc dot edu)
 */
public class Clock {

	private String descriptor;
	private long startMsecs;
	private long stopMsecs;
	
	public Clock(String descriptor) {
		this.descriptor = descriptor;
	}
	
	public void start() {
		this.startMsecs = System.currentTimeMillis();
	}
	
	public long elapsedSeconds() {
		return (stopMsecs - startMsecs) / 1000;
	}
	
	public void stopAndPrint() {
		this.stopMsecs = System.currentTimeMillis();
		
		Logger.info("Clock time in " + descriptor + ": " + elapsedSeconds());
	}
}

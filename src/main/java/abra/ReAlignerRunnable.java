/* Copyright 2013 University of North Carolina at Chapel Hill.  All rights reserved. */
package abra;

import java.util.List;

/**
 * Thread entry point for region specific processing.
 * 
 * @author Lisle E. Mose (lmose at unc dot edu)
 */
public class ReAlignerRunnable extends AbraRunnable {
	private List<Feature> regions;
	private ReAligner reAligner;
	private long createTime;
	
	public ReAlignerRunnable(ThreadManager threadManager, ReAligner reAligner, List<Feature> regions) {
		super(threadManager);
		createTime = System.nanoTime();
		this.regions = regions;
		this.reAligner = reAligner;
	}
	
	@Override
	public void go() throws Exception {
		long s = System.nanoTime();
		synchronized(AbraRunnable.class) {
			System.out.println("ReAlignerRunnable time_to_go: " + (s-createTime)/1000000);
			System.out.println("ReAlignerRunnable spawn_time: " + (s-spawnStartTime)/1000000);
		}
		
		for (Feature region : regions) {
			reAligner.processRegion(region);
		}
		long e = System.nanoTime();
		
		System.out.println("ReAlignerRunnable elapsed: " + (e-s)/1000000);
	}
}

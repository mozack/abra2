/* Copyright 2013 University of North Carolina at Chapel Hill.  All rights reserved. */
package abra;

public class SortBamRunnable implements Runnable {

	private ReAligner realigner;
	private String input;
	private String output;
	private String sortOrder;
	
	public SortBamRunnable(ReAligner realigner, String input, String output, String sortOrder) {
		this.realigner = realigner;
		this.input = input;
		this.output = output;
		this.sortOrder = sortOrder;
	}
	
	@Override
	public void run() {
		realigner.sortBam(input, output, sortOrder);		
	}	
}

package abra;

import java.io.IOException;


public class ReAlignerRunnable implements Runnable {

	private ReAligner reAligner;
	private Feature region;
	private String inputSam;
	private boolean isAlive = true;
	
	public ReAlignerRunnable(ReAligner reAligner, Feature region, String inputSam) {
		this.reAligner = reAligner;
		this.region = region;
		this.inputSam = inputSam;
	}
	
	@Override
	public void run() {
		
		try {
			reAligner.processRegion(region, inputSam);
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

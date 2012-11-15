package abra;

import java.io.IOException;

public class UpdateMismatchAndEditDistanceRunnable implements Runnable {
	
	private String inputBam;
	private String outputBam;
	private ReAligner realigner;

	public UpdateMismatchAndEditDistanceRunnable(String inputBam, String outputBam, ReAligner realigner) {
		this.inputBam = inputBam;
		this.outputBam = outputBam;
		this.realigner = realigner;
	}
	
	@Override
	public void run() {
		try {
			realigner.updateMismatchAndEditDistance(inputBam, outputBam);
		} catch (IOException e) {
			e.printStackTrace();
			throw new RuntimeException(e);
		}
	}
}

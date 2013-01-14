package abra;

import java.io.IOException;

@Deprecated
public class UpdateMismatchAndEditDistanceRunnable implements Runnable {
	
	private String inputBam;
	private String outputBam;
	private ReAligner realigner;

	@Override
	public void run() {
		throw new RuntimeException ("Deprecated.");
	}
}

package abra;

import java.io.IOException;

public class AdjustReadsRunnable implements Runnable {
	
	private ReAligner realigner;
	private String sortedOriginalReads;
	private String sortedAlignedToContig;
	private String outputSam;
	private boolean isTightAlignment;
	
	public AdjustReadsRunnable(ReAligner realigner, String sortedOriginalReads, String sortedAlignedToContig, String outputSam,
			boolean isTightAlignment) {
		this.realigner = realigner;
		this.sortedOriginalReads = sortedOriginalReads;
		this.sortedAlignedToContig = sortedAlignedToContig;
		this.outputSam = outputSam;
		this.isTightAlignment = isTightAlignment;
	}

	@Override
	public void run() {
		try {
			realigner.adjustReads(sortedOriginalReads, sortedAlignedToContig, outputSam, isTightAlignment);
		} catch (IOException e) {
			e.printStackTrace();
			throw new RuntimeException(e);
		}
	}
}

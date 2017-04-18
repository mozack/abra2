package abra;

public class SortedSAMWriterRunnable extends AbraRunnable {
	
	private SortedSAMWriter writer;
	private int sampleIdx;
	private String inputBam;

	public SortedSAMWriterRunnable(ThreadManager threadManager, SortedSAMWriter writer, int sampleIdx, String inputBam) {
		super(threadManager);
		this.writer = writer;
		this.sampleIdx = sampleIdx;
		this.inputBam = inputBam;
	}

	@Override
	public void go() throws Exception {
		writer.outputFinal(sampleIdx, inputBam);
	}

}

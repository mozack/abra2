package abra;

public class SortedSAMWriterRunnable extends AbraRunnable {
	
	private SortedSAMWriter writer;
	private int sampleIdx;

	public SortedSAMWriterRunnable(ThreadManager threadManager, SortedSAMWriter writer, int sampleIdx) {
		super(threadManager);
		this.writer = writer;
		this.sampleIdx = sampleIdx;
	}

	@Override
	public void go() throws Exception {
		writer.outputFinal(sampleIdx);
	}

}

package abra;

import htsjdk.samtools.SAMRecord;

import java.io.IOException;
import java.io.PipedInputStream;
import java.io.PipedOutputStream;
import java.util.Queue;

public class AlignerStdoutHandler implements StdoutHandler {
	
	private	Thread queueThread;
	private Thread adjustThread;
	private Thread draino;
	private AdjustReadsQueueRunnable adjustReadsRunnable;
	
	public AlignerStdoutHandler(AdjustReadsQueueRunnable adjustReadsRunnable) {
		this.adjustReadsRunnable = adjustReadsRunnable;
	}

	public void process(Process proc) throws IOException {
		
		PipedInputStream pis = new PipedInputStream(MAX_BYTES_TO_BUFFER);
		PipedOutputStream pos = new PipedOutputStream();
		pos.connect(pis);
		
		// Drain stdout and write to the piped output stream
		draino = new Thread(new Draino(proc.getInputStream(), pos));
		draino.start();
		
		// Read piped input stream and update read queue
		Queue<SAMRecord> queue = new ConcurrentQueue<SAMRecord>();
		queueThread = new Thread(new ReadInputStreamRunnable(pis, queue));
	 
		// Process read queue content
		adjustReadsRunnable.setReadQueue(queue);
		adjustThread = new Thread(adjustReadsRunnable);
		adjustThread.start();
		
		queueThread.start();
	}
	
	public void postProcess() throws InterruptedException {
		draino.join();
		queueThread.join();
		adjustReadsRunnable.setDone();
		adjustThread.join();
	}
}

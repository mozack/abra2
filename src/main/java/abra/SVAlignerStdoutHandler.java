package abra;

import htsjdk.samtools.DefaultSAMRecordFactory;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

import java.io.IOException;
import java.io.PipedInputStream;
import java.io.PipedOutputStream;

public class SVAlignerStdoutHandler implements StdoutHandler{
	
	private Thread draino;
	private int readLength;
	private SAMFileHeader header;
	private Thread counterThread;
	private SVReadCounterRunnable counter; 
	
	public SVAlignerStdoutHandler(int readLength, SAMFileHeader header) {
		this.readLength = readLength;
		this.header = header;
	}

	@Override
	public void process(Process proc) throws IOException {
		
		PipedInputStream pis = new PipedInputStream(MAX_BYTES_TO_BUFFER);
		PipedOutputStream pos = new PipedOutputStream();
		pos.connect(pis);
		
		// Drain stdout and write to the piped output stream
		draino = new Thread(new Draino(proc.getInputStream(), pos));
		draino.start();
		
		final SamReader reader =
		        SamReaderFactory.make()
		                .validationStringency(ValidationStringency.SILENT)
		                .samRecordFactory(DefaultSAMRecordFactory.getInstance())
		                .open(SamInputResource.of(pis));
		
		counter = new SVReadCounterRunnable(reader, readLength, header);
		
		counterThread = new Thread(counter);
		counterThread.start();
	}

	@Override
	public void postProcess() throws InterruptedException {
		
		draino.join();
		counterThread.join();
	}

	public SVReadCounter getCounter() {
		return counter.getCounter();
	}
}

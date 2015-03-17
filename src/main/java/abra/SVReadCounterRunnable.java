package abra;


import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SamReader;

public class SVReadCounterRunnable implements Runnable {

	private int readLength;
	private SAMFileHeader header;
	private SamReader reader;
	private SVReadCounter counter;
	
	public SVReadCounterRunnable(SamReader reader, int readLength, SAMFileHeader header) {
		this.readLength = readLength;
		this.header = header;
		this.reader = reader;
	}
	
	@Override
	public void run() {
		counter = new SVReadCounter();
		counter.countReadsSupportingBreakpoints(reader, readLength, header);
	}

	public SVReadCounter getCounter() {
		return counter;
	}
}

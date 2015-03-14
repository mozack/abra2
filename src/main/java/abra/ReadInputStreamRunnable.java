package abra;

import htsjdk.samtools.DefaultSAMRecordFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

import java.io.InputStream;
import java.util.Queue;

public class ReadInputStreamRunnable extends AbraRunnable {

	private static final int MAX_QUEUE_SIZE = 600000;
	
	private InputStream is;
	private Queue<SAMRecord> queue;
	
	public ReadInputStreamRunnable(ThreadManager threadManager, InputStream is, Queue<SAMRecord> queue) {
		super(threadManager);
		this.is = is;
		this.queue = queue;
	}

	@Override
	public void go() throws Exception {
		
		final SamReader reader =
		        SamReaderFactory.make()
		                .validationStringency(ValidationStringency.SILENT)
		                .samRecordFactory(DefaultSAMRecordFactory.getInstance())
		                .open(SamInputResource.of(is));
	
		
		for (SAMRecord read : reader) {
			while (queue.size() > MAX_QUEUE_SIZE) {
				System.out.println("Queue too big");
				Thread.sleep(100);
			}
			queue.add(read);
		}
		
		reader.close();
	}	
}

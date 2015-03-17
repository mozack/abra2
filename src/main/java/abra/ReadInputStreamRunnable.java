package abra;

import htsjdk.samtools.DefaultSAMRecordFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

import java.io.IOException;
import java.io.InputStream;
import java.util.Queue;

public class ReadInputStreamRunnable implements Runnable {

	private static final int MAX_QUEUE_SIZE = 600000;
	
	private InputStream is;
	private Queue<SAMRecord> queue;
	
	public ReadInputStreamRunnable(InputStream is, Queue<SAMRecord> queue) {
		this.is = is;
		this.queue = queue;
	}

	public void run() {
		
		final SamReader reader =
		        SamReaderFactory.make()
		                .validationStringency(ValidationStringency.SILENT)
		                .samRecordFactory(DefaultSAMRecordFactory.getInstance())
		                .open(SamInputResource.of(is));
	
		
		for (SAMRecord read : reader) {
			while (queue.size() > MAX_QUEUE_SIZE) {
//				System.out.println("Queue too big");
				try {
					Thread.sleep(100);
				} catch (InterruptedException e) { }
			}
			queue.add(read);
		}
		
		try {
			reader.close();
		} catch (IOException e) {
			e.printStackTrace();
			throw new RuntimeException(e);
		}
	}	
}

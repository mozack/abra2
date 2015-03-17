package abra;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.ListIterator;
import java.util.Random;

import htsjdk.samtools.SAMRecord;

/**
 * Downsamples reads on the fly.  Useful for when depth gets too high.
 */
public class DownsampledReadList {
	
	private Random random = new Random(1);
	private List<SAMRecord> reads = new ArrayList<SAMRecord>();
	private int maxReads;
	private int totalReads = 0;

	public DownsampledReadList(int maxReads) {
		this.maxReads = maxReads;
	}
	
	public void add(SAMRecord read) {
		totalReads += 1;
		if (reads.size() < maxReads) {
			reads.add(read);
		} else {
			int slot = random.nextInt(totalReads);
			if (slot < maxReads) {
				reads.set(slot, read);
			}
		}
	}
	
	public List<SAMRecord> getReads() {
		return reads;
	}
	
	public int getTotalReadCount() {
		return totalReads;
	}
}

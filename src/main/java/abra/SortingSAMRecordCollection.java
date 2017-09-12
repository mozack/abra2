package abra;

import java.io.File;
import java.util.Iterator;

import htsjdk.samtools.BAMRecordCodec;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordCoordinateComparator;
import htsjdk.samtools.SAMRecordQueryNameComparator;
import htsjdk.samtools.util.SortingCollection2;

/**
 * Wrapper class for HTSJDK's SortingCollection specific to SAMRecord and adds size tracking
 *
 * @author lmose
 */
public class SortingSAMRecordCollection implements Iterable<SAMRecord> {
	
	private SortingCollection2<SAMRecord> reads;
	private int size = 0;
	
	// Array backing the SortingCollection.
	// Re-using this avoids the cost of reallocating the large array each time the SortingCollection is flushed
//	private SAMRecord[] records;
	
	public static SortingSAMRecordCollection newSortByCoordinateInstance(SAMRecord[] recordArray, SAMFileHeader header, int maxRecordsInRAM, String tempDir) {
		return new SortingSAMRecordCollection(recordArray, header, new SAMRecordCoordinateComparator(), maxRecordsInRAM, tempDir);
	}
	
	public static SortingSAMRecordCollection newSortByNameInstance(SAMRecord[] recordArray, SAMFileHeader header, int maxRecordsInRAM, String tempDir) {
		return new SortingSAMRecordCollection(recordArray, header, new SAMRecordQueryNameComparator(), maxRecordsInRAM, tempDir);
	}
	
	private SortingSAMRecordCollection(SAMRecord[] recordArray, SAMFileHeader header, java.util.Comparator<SAMRecord> comparator, int maxRecordsInRAM, String tempDir) {
		reads = SortingCollection2.newInstance(recordArray, SAMRecord.class, new BAMRecordCodec(header), comparator, maxRecordsInRAM, new File(tempDir));
	}

	@Override
	public Iterator<SAMRecord> iterator() {
		return reads.iterator();
	}
	
	public void add(SAMRecord read) {
		reads.add(read);
		size += 1;
	}
	
	public int size() {
		return size;
	}
	
	public void cleanup() {
		reads.cleanup();
	}
}

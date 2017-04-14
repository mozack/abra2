package abra;

import static org.testng.Assert.assertEquals;

import java.util.ArrayList;
import java.util.List;

import org.testng.annotations.Test;

// TODO: Test more thoroughly including multiple samples
public class MultiSamReaderTest {

	@Test (groups = "unit")
	public void testReadTwoBams() {
		String[] bams = new String[] { "test-data/sample1.bam", "test-data/sample2.bam" };
		int minMapq = 20;
		boolean isPairedEnd = false;
		String chromosome = "10";
		Feature chromosomeChunk = new Feature(chromosome, 3000000, 4000000);
		
		MultiSamReader rdr = new MultiSamReader(bams, minMapq, isPairedEnd, chromosomeChunk);
		
		List<SAMRecordWrapper> reads = new ArrayList<SAMRecordWrapper>();
		for (SAMRecordWrapper read : rdr) {
			reads.add(read);
		}
		
		assertEquals(reads.size(), 9);
		assertEquals(reads.get(0).getSamRecord().getReadName(), "ERR194161.458962555");
		assertEquals(reads.get(0).getSampleIdx(), 0);
		
		assertEquals(reads.get(1).getSamRecord().getReadName(), "ERR194161.458962561");
		assertEquals(reads.get(1).getSampleIdx(), 1);
		
		assertEquals(reads.get(2).getSamRecord().getReadName(), "ERR194161.458962547");
		assertEquals(reads.get(2).getSampleIdx(), 1);
		
		assertEquals(reads.get(3).getSamRecord().getReadName(), "ERR194161.458962586");
		assertEquals(reads.get(3).getSampleIdx(), 0);
		
		assertEquals(reads.get(4).getSamRecord().getReadName(), "ERR194161.458962577");
		assertEquals(reads.get(4).getSampleIdx(), 0);
		
		assertEquals(reads.get(5).getSamRecord().getReadName(), "ERR194161.458962575");
		assertEquals(reads.get(5).getSampleIdx(), 1);
		
		assertEquals(reads.get(6).getSamRecord().getReadName(), "ERR194161.458962567");
		assertEquals(reads.get(6).getSampleIdx(), 0);
		
		// Same read in both samples here.  Order doesn't matter.
		assertEquals(reads.get(7).getSamRecord().getReadName(), "ERR194161.458962591");
		assertEquals(reads.get(8).getSamRecord().getReadName(), "ERR194161.458962591");
	}
}

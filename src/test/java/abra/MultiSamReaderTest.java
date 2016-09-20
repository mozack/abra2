package abra;

import static org.testng.Assert.assertEquals;

import org.testng.annotations.Test;

// TODO: Test more thoroughly including multiple samples
public class MultiSamReaderTest {

	@Test (groups = "unit")
	public void testReadSingleBam() {
		String[] bam = new String[] { "demo/abra_demo.bam" };
		int minMapq = 20;
		boolean isPairedEnd = false;
		
		MultiSamReader rdr = new MultiSamReader(bam, minMapq, isPairedEnd);
		
		int numReads = 0;
		for (SAMRecordWrapper read : rdr) {
			numReads += 1;
		}
		
		assertEquals(numReads, 202);
	}
}

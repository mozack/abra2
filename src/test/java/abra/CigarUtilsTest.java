package abra;

import static org.testng.Assert.assertEquals;

import org.testng.annotations.Test;

public class CigarUtilsTest {

	@Test (groups = "unit")
	public void testMatchSubset() {
		String contigCigar = "500M";
		String readContig = CigarUtils.subsetCigarString(50, 100, contigCigar);
		assertEquals(readContig, "100M");
	}
	
	@Test (groups = "unit")
	public void testExactMatch() {
		String contigCigar = "100M";
		String readContig = CigarUtils.subsetCigarString(0, 100, contigCigar);
		assertEquals(readContig, "100M");
	}
	
	@Test (groups = "unit")
	public void testDeletion() {
		String contigCigar = "10M50D490M";
		String readContig = CigarUtils.subsetCigarString(5, 100, contigCigar);
		assertEquals(readContig, "5M50D95M");
	}
	
	@Test (groups = "unit")
	public void testSkipDeletionAtStart() {
		String contigCigar = "10M50D490M";
		String readContig = CigarUtils.subsetCigarString(10, 100, contigCigar);
		assertEquals(readContig, "100M");
	}
	
	@Test (groups = "unit")
	public void testSkipDeletionWithinStart() {
		String contigCigar = "10M50D490M";
		String readContig = CigarUtils.subsetCigarString(30, 100, contigCigar);
		assertEquals(readContig, "100M");
	}

	@Test (groups = "unit")
	public void testInsertion() {
		String contigCigar = "10M50I490M";
		String readContig = CigarUtils.subsetCigarString(5, 100, contigCigar);
		assertEquals(readContig, "5M50I45M");
	}
	
	@Test (groups = "unit")
	public void testInsertionAtStart() {
		String contigCigar = "10M50I490M";
		String readContig = CigarUtils.subsetCigarString(10, 100, contigCigar);
		assertEquals(readContig, "50I50M");
	}
	
	@Test (groups = "unit")
	public void testInsertionPartialOverlapAtStart() {
		String contigCigar = "10M50I490M";
		String readContig = CigarUtils.subsetCigarString(30, 100, contigCigar);
		assertEquals(readContig, "30I70M");
	}
	
	@Test (groups = "unit")
	public void testSkipDeletionAtEnd() {
		String contigCigar = "300M50D200M";
		String readContig = CigarUtils.subsetCigarString(200, 100, contigCigar);
		assertEquals(readContig, "100M");
	}
	
	@Test (groups = "unit")
	public void testSkipInsertionAtEnd() {
		String contigCigar = "300M50I200M";
		String readContig = CigarUtils.subsetCigarString(200, 100, contigCigar);
		assertEquals(readContig, "100M");
	}
	
	@Test (groups = "unit")
	public void testInsertionPartialOverlapAtEnd() {
		String contigCigar = "300M50I200M";
		String readContig = CigarUtils.subsetCigarString(201, 100, contigCigar);
		assertEquals(readContig, "99M1I");
	}
}

/* Copyright 2013 University of North Carolina at Chapel Hill.  All rights reserved. */
package abra;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMRecord;

import org.testng.Assert;
import org.testng.annotations.Test;

public class SAMRecordUtilsTest {

	@Test (groups = "unit" )
	public void testRemoveSoftClips_withDeletionAndSoftClipAtEnd() {
		SAMRecord read = new SAMRecord(null);
		
		read.setReadName("TEST1");
		read.setCigarString("5M2D3M2S");
		read.setReadString("CCCCCCAGCC");
		
		SAMRecordUtils.removeSoftClips(read);
		
		Assert.assertEquals(read.getCigarString(), "5M2D3M");
		Assert.assertEquals(read.getReadString(), "CCCCCCAG");
	}
	
	@Test (groups = "unit" )
	public void testRemoveSoftClips_withDeletionAndSoftClipAtStart() {
		SAMRecord read = new SAMRecord(null);
		
		read.setReadName("TEST1");
		read.setCigarString("2S5M2D3M");
		read.setReadString("CCCCCCAGCC");
		
		SAMRecordUtils.removeSoftClips(read);
		
		Assert.assertEquals(read.getCigarString(), "5M2D3M");
		Assert.assertEquals(read.getReadString(), "CCCCAGCC");
	}
	
	@Test (groups = "unit")
	public void testSubsetCigar_exact() {
		Cigar cigar = getCigar("100M");
		Cigar cigar2 = SAMRecordUtils.subset(cigar, 0, 99);
		Assert.assertEquals(cigar2.toString(), "100M");
	}
	
	@Test (groups = "unit")
	public void testSubsetCigar_lessThan1Elem() {
		Cigar cigar = getCigar("100M");
		Cigar cigar2 = SAMRecordUtils.subset(cigar, 0, 90);
		Assert.assertEquals(cigar2.toString(), "91M");
	}
	
	@Test (groups = "unit")
	public void testSubsetCigar_lessThan1Elem2() {
		Cigar cigar = getCigar("100M");
		Cigar cigar2 = SAMRecordUtils.subset(cigar, 10, 99);
		Assert.assertEquals(cigar2.toString(), "90M");
	}
	
	@Test (groups = "unit")
	public void testSubsetCigar_lessThan1Elem3() {
		Cigar cigar = getCigar("100M");
		Cigar cigar2 = SAMRecordUtils.subset(cigar, 10, 90);
		Assert.assertEquals(cigar2.toString(), "81M");
	}
	
	@Test (groups = "unit")
	public void testSubsetCigar_lessThanFirstElem() {
		Cigar cigar = getCigar("100M100I100M");
		Cigar cigar2 = SAMRecordUtils.subset(cigar, 10, 90);
		Assert.assertEquals(cigar2.toString(), "81M");
	}
	
	@Test (groups = "unit")
	public void testSubsetCigar_spanElems1() {
		Cigar cigar = getCigar("100M100I100M");
		Cigar cigar2 = SAMRecordUtils.subset(cigar, 0, 100);
		Assert.assertEquals(cigar2.toString(), "100M1I");
	}
	
	@Test (groups = "unit")
	public void testSubsetCigar_spanElems2() {
		Cigar cigar = getCigar("100M100I100M");
		Cigar cigar2 = SAMRecordUtils.subset(cigar, 1, 100);
		Assert.assertEquals(cigar2.toString(), "99M1I");
	}
	
	@Test (groups = "unit")
	public void testSubsetCigar_spanElems3() {
		Cigar cigar = getCigar("100M100I100M");
		Cigar cigar2 = SAMRecordUtils.subset(cigar, 1, 199);
		Assert.assertEquals(cigar2.toString(), "99M100I");
	}
	
	@Test (groups = "unit")
	public void testSubsetCigar_skipElem() {
		Cigar cigar = getCigar("100M100I100M");
		Cigar cigar2 = SAMRecordUtils.subset(cigar, 100, 199);
		Assert.assertEquals(cigar2.toString(), "100I");
	}
	
	@Test (groups = "unit")
	public void testSubsetCigar_skipElem2() {
		Cigar cigar = getCigar("100M100I100M");
		Cigar cigar2 = SAMRecordUtils.subset(cigar, 200, 299);
		Assert.assertEquals(cigar2.toString(), "100M");
	}

	@Test (groups = "unit")
	public void testSubsetCigar_skipAndSpan1() {
		Cigar cigar = getCigar("100M100I100M");
		Cigar cigar2 = SAMRecordUtils.subset(cigar, 100, 249);
		Assert.assertEquals(cigar2.toString(), "100I50M");
	}
	
	@Test (groups = "unit")
	public void testSubsetCigar_skipAndSpan2() {
		Cigar cigar = getCigar("100M100I100M");
		Cigar cigar2 = SAMRecordUtils.subset(cigar, 150, 249);
		Assert.assertEquals(cigar2.toString(), "50I50M");
	}
	
	@Test (groups = "unit")
	public void testSubsetCigar_spanAll1() {
		Cigar cigar = getCigar("100M100I100M");
		Cigar cigar2 = SAMRecordUtils.subset(cigar, 0, 299);
		Assert.assertEquals(cigar2.toString(), "100M100I100M");
	}
	
	@Test (groups = "unit")
	public void testSubsetCigar_spanAll2() {
		Cigar cigar = getCigar("100M100I100M");
		Cigar cigar2 = SAMRecordUtils.subset(cigar, 0, 298);
		Assert.assertEquals(cigar2.toString(), "100M100I99M");
	}
	
	@Test (groups = "unit")
	public void testSubsetCigar_spanAll3() {
		Cigar cigar = getCigar("100M100I100M");
		Cigar cigar2 = SAMRecordUtils.subset(cigar, 1, 299);
		Assert.assertEquals(cigar2.toString(), "99M100I100M");
	}
	
	@Test (groups = "unit")
	public void testSubsetCigar_spanAll4() {
		Cigar cigar = getCigar("100M100I100M");
		Cigar cigar2 = SAMRecordUtils.subset(cigar, 1, 298);
		Assert.assertEquals(cigar2.toString(), "99M100I99M");
	}
	
	@Test (groups = "unit")
	public void testSubsetCigar_long() {
		Cigar cigar = getCigar("100M100I100M100I100M");
		Cigar cigar2 = SAMRecordUtils.subset(cigar, 1, 498);
		Assert.assertEquals(cigar2.toString(), "99M100I100M100I99M");
	}
	
	@Test (groups = "unit")
	public void testSubsetCigar_deletion() {
		Cigar cigar = getCigar("100M100D100M");
		Cigar cigar2 = SAMRecordUtils.subset(cigar, 1, 198);
		Assert.assertEquals(cigar2.toString(), "99M100D99M");
	}
	
	private Cigar getCigar(String str) {
		SAMRecord read = new SAMRecord(null);
		read.setCigarString(str);
		return read.getCigar();
	}
}

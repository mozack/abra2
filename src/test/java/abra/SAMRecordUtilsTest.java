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
	
	@Test (groups = "unit")
	public void testGetLeadingClips() {
		SAMRecord read = new SAMRecord(null);
		Cigar cigar = getCigar("10S15M8S");
		read.setCigar(cigar);
		String leading = SAMRecordUtils.getLeadingClips(read);
		Assert.assertEquals(leading, "10S");
	}
	
	@Test (groups = "unit")
	public void testGetLeadingClips_empty() {
		SAMRecord read = new SAMRecord(null);
		Cigar cigar = getCigar("15M8S");
		read.setCigar(cigar);
		String leading = SAMRecordUtils.getLeadingClips(read);
		Assert.assertEquals(leading, "");
	}
	
	@Test (groups = "unit")
	public void testGetTrailingClips() {
		SAMRecord read = new SAMRecord(null);
		Cigar cigar = getCigar("10S15M8S");
		read.setCigar(cigar);
		String leading = SAMRecordUtils.getTrailingClips(read);
		Assert.assertEquals(leading, "8S");
	}

	@Test (groups = "unit")
	public void testGetTrailingClips_empty() {
		SAMRecord read = new SAMRecord(null);
		Cigar cigar = getCigar("10S15M");
		read.setCigar(cigar);
		String leading = SAMRecordUtils.getTrailingClips(read);
		Assert.assertEquals(leading, "");
	}
	
	@Test (groups = "unit")
	public void testGetMappedReadPortion() {
		SAMRecord read = new SAMRecord(null);
		Cigar cigar = getCigar("4S5M7S");
		String seq = "ATCGCGATACCCCCCC";
		read.setCigar(cigar);
		read.setReadString(seq);
		
		String unclipped = SAMRecordUtils.getMappedReadPortion(read);
		Assert.assertEquals(unclipped, "CGATA");
	}
	
	private Cigar getCigar(String str) {
		SAMRecord read = new SAMRecord(null);
		read.setCigarString(str);
		return read.getCigar();
	}
	
	@Test (groups = "unit")
	public void testMergeSequences() {
		String s1 = "ATCGATTACGACTTGGGCAA";
		String s2 =    "GATTACGACTTGGGCAACGC";
		String q1 = "01234567890123456789";
		String q2 =    "ABCDEF!ABCDEF!FF!!!!";
		
		Pair<String, String> merged = SAMRecordUtils.mergeSequences(s1, s2, q1, q2);
		Assert.assertEquals(merged.getFirst(),  "ATCGATTACGACTTGGGCAACGC");
		Assert.assertEquals(merged.getSecond(), "012ABCDEF9ABCDEF6FF9!!!");
	}
	
	@Test (groups = "unit")
	public void testMergeSequences_discordantMismatch() {
		String s1 = "ATCGATTACGACTTGGGCTTTAA";
		String s2 =    "GATTACGACTTGGGCCTTAACGC";
		String q1 = "AAAAAAAAAAAAAAAAAAAAAAA";
		String q2 =    "BBBBBBBBBBBBBBBBBBBBBBB";
		
		Pair<String, String> merged = SAMRecordUtils.mergeSequences(s1, s2, q1, q2);
		Assert.assertEquals(merged.getFirst(),  "ATCGATTACGACTTGGGCNTTAACGC");
		Assert.assertEquals(merged.getSecond(), "AAABBBBBBBBBBBBBBB!BBBBBBB");
	}
	
	@Test (groups = "unit")
	public void testMergeSequences_preferredMismatch() {
		String s1 = "ATCGATTACGACTTGGGCTTTAA";
		String s2 =    "GATTACGACTTGGGCCTTAACGC";
		String q1 = "#######################";
		String q2 =    "BBBBBBBBBBBBBBBBBBBBBBB";
		
		Pair<String, String> merged = SAMRecordUtils.mergeSequences(s1, s2, q1, q2);
		Assert.assertEquals(merged.getFirst(),  "ATCGATTACGACTTGGGCCTTAACGC");
		Assert.assertEquals(merged.getSecond(), "###BBBBBBBBBBBBBBBBBBBBBBB");
	}
	
	@Test (groups = "unit")
	public void testMergeSequences_tooManyMismatches() {
		String s1 = "ATCGATTACGACTTGGGCTTTAA";
		String s2 =    "GATTACGACTTGGGCCCCTAACGC";
		String q1 = "#######################";
		String q2 =    "BBBBBBBBBBBBBBBBBBBBBBB";
		
		Pair<String, String> merged = SAMRecordUtils.mergeSequences(s1, s2, q1, q2);
		Assert.assertEquals(merged,  null);
	}
	
//	@Test (groups = "unit")
//	public void testMergeSequences_mismatchNearStart() {
//		String s1 = "ATCGAATACGACTTGGGCAA";
//		String s2 =    "GATTACGACTTGGGCAACGC";
//		String q1 = "01234567890123456789";
//		String q2 =    "ABCDEF!ABCDEF!FF!!!!";
//		
//		Pair<String, String> merged = SAMRecordUtils.mergeSequences(s1, s2, q1, q2);
//		Assert.assertEquals(merged.getFirst(),  "ATCGATTACGACTTGGGCAACGC");
//		Assert.assertEquals(merged.getSecond(), "012ABCDEF9ABCDEF6FF9!!!");
//	}
	
	@Test (groups = "unit")
	public void testMergeSequences_multipleHeadHits() {
		String s1 = "AAATCAGGCTAGGCTAAGCTATGATGTTCCTTAGATTAGGTGTATTAAATCCATTTTCAACTTACAATATTTTCAACTTACGACGAGTTTATCAGGAAGT";
		String s2 = "TTCAACTTACAATATTTTCAACTTACGACGAGTTTATCAGGAAGTAACACCATCGTAAGTCAAGTAGCATCTGTATCAGGCAAAGTCATAGAACCATTTT";
		String q1 = "<<AAAFFFFFFFFFFFFFFFAFFFFFFFFFFF<FAFFAFFFFFFFFAA<FFF.FFFFAF7FFAFFFFFFFFFFFFFFFFFFFFF<FFFFFFFFAFFFFFF";
		String q2 = "AFFFFF<AF<AFAFFFFFFFFFF<FFF7FFFFFFFFFFFFFFAFFFAFFFFFFFFFFFF.FFFFFFFFFAFFFFFFF<FFFFFFFAFFFFFFFFFAA<AA";
		
		Pair<String, String> merged = SAMRecordUtils.mergeSequences(s1, s2, q1, q2);
		
		Assert.assertEquals(merged.getFirst(), "AAATCAGGCTAGGCTAAGCTATGATGTTCCTTAGATTAGGTGTATTAAATCCATTTTCAACTTACAATATTTTCAACTTACGACGAGTTTATCAGGAAGTAACACCATCGTAAGTCAAGTAGCATCTGTATCAGGCAAAGTCATAGAACCATTTT");
		Assert.assertEquals(merged.getSecond(), "<<AAAFFFFFFFFFFFFFFFAFFFFFFFFFFF<FAFFAFFFFFFFFAA<FFF.FFFFFFFFFAFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFAFFFFFFFFFFFF.FFFFFFFFFAFFFFFFF<FFFFFFFAFFFFFFFFFAA<AA");
	}
}

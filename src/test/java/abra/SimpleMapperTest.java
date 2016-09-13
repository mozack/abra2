package abra;

import static org.testng.Assert.assertEquals;
import org.testng.annotations.Test;
import abra.SimpleMapper.SimpleMapperResult;


public class SimpleMapperTest {
	
	private String contig1 = "TTCAACTAGAGAGAGGTAAAAATTTTTCTAGAACATGAATTGCCCACTCCCCTCATTCCTTCTCAGAAACTAACTGAATTCCAGTGGGTGTGCCTGGCAAACCCAAAAGCAGTTTCTGTTCAGGATGCTGGTCTTACCTGTGAAGGCGTTCATGAACGTGGAGAGGGACCGGTTCAACATTTTGAAGAAAGGGTCTCTGCACGGATATTTCTGAGACCCACAAAGGACGGTATGCTCAAGAATGTGAGGAACACCAGTACTGTCCATGGGAGTGGTACGGAACTGCACGCTAGGGAAGAGAGAGGAATGGCACGCTAGGGAAGGCGAATGACCAGAACGCAAAAGGTTCAGCTTAGTGCTGCGGACACAGTTCCCAGATGCATCATCACCTCAGGCTACTAGAAATCATCATTCTGACACCACAATCCTCCAGCACAGGGTTTTCCAACTATA";
	private String contig2 = "ATCGATCGATATCGATCGATATCGATCGATATCGATCGATATCGATCGATATCGATCGATATCGATCGATATCGATCGATATCGATCGATATCGATCGATATCGATCGATATCGATCGATATCGATCGATATCGATCGATATCGATCGATATCGATCGATATCGATCGATATCGATCGATATCGATCGAT";

	@Test (groups = "unit" )
	public void testMapExact() {
		SimpleMapper sm = new SimpleMapper(contig1);
		String read = "TACTGTCCATGGGAGTGGTACGGAACTGCACGCTAGGGAAGAGAGAGGAATGGCACGCTAGGGAAGGCGAATGACCAGAACGCAAAAGGTTCAGCTTAGTG";
		SimpleMapperResult smr = sm.map(read);
		assertEquals(0, smr.getMismatches());
		assertEquals(257, smr.getPos());
	}
	
	@Test (groups = "unit" )
	public void testMapOneMismatch() {
		SimpleMapper sm = new SimpleMapper(contig1);
		String read = "TACTGTCCATGGGAGTGCTACGGAACTGCACGCTAGGGAAGAGAGAGGAATGGCACGCTAGGGAAGGCGAATGACCAGAACGCAAAAGGTTCAGCTTAGTG";
		SimpleMapperResult smr = sm.map(read);
		assertEquals(1, smr.getMismatches());
		assertEquals(257, smr.getPos());
	}
	
	@Test (groups = "unit" )
	public void testMapFiveMismatches() {
		SimpleMapper sm = new SimpleMapper(contig1);
		String read = "AACTGTCCATGGGAGTGGTACGTTTCTGCACGCTAGGGAAGAGAGAGGAATGGCACGCTAGGGAAGGCGAATGACCAGAACGCAAAAGGTTCAGCTTAGTC";
		SimpleMapperResult smr = sm.map(read);
		assertEquals(5, smr.getMismatches());
		assertEquals(257, smr.getPos());
	}
	
	@Test (groups = "unit" )
	public void testMapSixMismatches() {
		SimpleMapper sm = new SimpleMapper(contig1);
		String read = "AACTGTCCATGGGAGTGGTACGTTTCTGCACGCTAGGGAAGAGAGAGGAAAGGCACGCTAGGGAAGGCGAATGACCAGAACGCAAAAGGTTCAGCTTAGTC";
		SimpleMapperResult smr = sm.map(read);
		assertEquals(SimpleMapper.UNMAPPED, smr.getPos());
	}
	
	@Test (groups = "unit" )
	public void testMapNoSeedMatch() {
		SimpleMapper sm = new SimpleMapper(contig1);
		String read = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
		SimpleMapperResult smr = sm.map(read);
		assertEquals(SimpleMapper.UNMAPPED, smr.getPos());
	}
	
	@Test (groups = "unit" )
	public void testMapAmbiguousMatch() {
		SimpleMapper sm = new SimpleMapper(contig2);
		String read = "CGATCGATATCGATCGATATCGATCGATATCGATCGATATCGATCGATATCGATCGATATCGATCGATATCGATCGATATCGATCGATATCGATCGATA";
		SimpleMapperResult smr = sm.map(read);
		assertEquals(SimpleMapper.HOMOLOGOUS_MAPPING, smr.getPos());
		assertEquals(0, smr.getMismatches());
	}
	
	@Test (groups = "unit" )
	public void testMapAmbiguousMatchWithMismatches() {
		SimpleMapper sm = new SimpleMapper(contig2);
		String read = "CGATCGATATCGATCGATAACGATCGATATCGATCGATATCGATCGATATCGATCGATATCGATCGATATCGATCGATATCGATCGATTACGATCGATA";
		SimpleMapperResult smr = sm.map(read);
		assertEquals(SimpleMapper.HOMOLOGOUS_MAPPING, smr.getPos());
		assertEquals(3, smr.getMismatches());
	}
	
	@Test (groups = "unit" )
	public void testSimple1Mismatch() {
		SimpleMapper sm = new SimpleMapper("ATCGAAAAAATTTTTTCCCCCCGGGGGGATCGGCTAATCG");
		String read =                          "ATAAAATTTTTTCCCCCCGGGGGGATCG";
	
		SimpleMapperResult smr = sm.map(read);
		assertEquals(4, smr.getPos());
		assertEquals(1, smr.getMismatches());
		
	}
}

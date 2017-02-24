package abra;

import static org.testng.Assert.assertEquals;
import org.testng.annotations.Test;

import abra.SimpleMapper.Orientation;
import abra.SimpleMapper.SimpleMapperResult;


public class SimpleMapperTest {
	
	private String contig1 = "TTCAACTAGAGAGAGGTAAAAATTTTTCTAGAACATGAATTGCCCACTCCCCTCATTCCTTCTCAGAAACTAACTGAATTCCAGTGGGTGTGCCTGGCAAACCCAAAAGCAGTTTCTGTTCAGGATGCTGGTCTTACCTGTGAAGGCGTTCATGAACGTGGAGAGGGACCGGTTCAACATTTTGAAGAAAGGGTCTCTGCACGGATATTTCTGAGACCCACAAAGGACGGTATGCTCAAGAATGTGAGGAACACCAGTACTGTCCATGGGAGTGGTACGGAACTGCACGCTAGGGAAGAGAGAGGAATGGCACGCTAGGGAAGGCGAATGACCAGAACGCAAAAGGTTCAGCTTAGTGCTGCGGACACAGTTCCCAGATGCATCATCACCTCAGGCTACTAGAAATCATCATTCTGACACCACAATCCTCCAGCACAGGGTTTTCCAACTATA";
	private String contig2 = "ATCGATCGATATCGATCGATATCGATCGATATCGATCGATATCGATCGATATCGATCGATATCGATCGATATCGATCGATATCGATCGATATCGATCGATATCGATCGATATCGATCGATATCGATCGATATCGATCGATATCGATCGATATCGATCGATATCGATCGATATCGATCGATATCGATCGAT";
	
	private static final double DEFAULT_MISMATCH_RATE = .05;

	@Test (groups = "unit" )
	public void testMapExact() {
		SimpleMapper sm = new SimpleMapper(contig1, DEFAULT_MISMATCH_RATE);
		String read = "TACTGTCCATGGGAGTGGTACGGAACTGCACGCTAGGGAAGAGAGAGGAATGGCACGCTAGGGAAGGCGAATGACCAGAACGCAAAAGGTTCAGCTTAGTG";
		SimpleMapperResult smr = sm.map(read);
		assertEquals(0, smr.getMismatches());
		assertEquals(257, smr.getPos());
		assertEquals(Orientation.FORWARD, smr.getOrientation());
	}
	
	@Test (groups = "unit" )
	public void testMapOneMismatch() {
		SimpleMapper sm = new SimpleMapper(contig1, DEFAULT_MISMATCH_RATE);
		String read = "TACTGTCCATGGGAGTGCTACGGAACTGCACGCTAGGGAAGAGAGAGGAATGGCACGCTAGGGAAGGCGAATGACCAGAACGCAAAAGGTTCAGCTTAGTG";
		SimpleMapperResult smr = sm.map(read);
		assertEquals(1, smr.getMismatches());
		assertEquals(257, smr.getPos());
		assertEquals(Orientation.FORWARD, smr.getOrientation());
	}
	
	@Test (groups = "unit" )
	public void testMapFiveMismatches() {
		SimpleMapper sm = new SimpleMapper(contig1, DEFAULT_MISMATCH_RATE);
		String read = "AACTGTCCATGGGAGTGGTACGTTTCTGCACGCTAGGGAAGAGAGAGGAATGGCACGCTAGGGAAGGCGAATGACCAGAACGCAAAAGGTTCAGCTTAGTC";
		SimpleMapperResult smr = sm.map(read);
		assertEquals(5, smr.getMismatches());
		assertEquals(257, smr.getPos());
		assertEquals(Orientation.FORWARD, smr.getOrientation());
	}
	
	@Test (groups = "unit" )
	public void testMapSixMismatches() {
		SimpleMapper sm = new SimpleMapper(contig1, DEFAULT_MISMATCH_RATE);
		String read = "AACTGTCCATGGGAGTGGTACGTTTCTGCACGCTAGGGAAGAGAGAGGAAAGGCACGCTAGGGAAGGCGAATGACCAGAACGCAAAAGGTTCAGCTTAGTC";
		SimpleMapperResult smr = sm.map(read);
		assertEquals(SimpleMapper.UNMAPPED, smr.getPos());
	}
	
	@Test (groups = "unit" )
	public void testMapSixMismatchesIncreasedMismatchRate() {
		SimpleMapper sm = new SimpleMapper(contig1, .06);
		String read = "AACTGTCCATGGGAGTGGTACGTTTCTGCACGCTAGGGAAGAGAGAGGAAAGGCACGCTAGGGAAGGCGAATGACCAGAACGCAAAAGGTTCAGCTTAGTC";
		SimpleMapperResult smr = sm.map(read);
		assertEquals(6, smr.getMismatches());
		assertEquals(257, smr.getPos());
	}

	
	@Test (groups = "unit" )
	public void testMapNoSeedMatch() {
		SimpleMapper sm = new SimpleMapper(contig1, DEFAULT_MISMATCH_RATE);
		String read = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
		SimpleMapperResult smr = sm.map(read);
		assertEquals(SimpleMapper.UNMAPPED, smr.getPos());
	}
	
	@Test (groups = "unit" )
	public void testMapAmbiguousMatch() {
		SimpleMapper sm = new SimpleMapper(contig2, DEFAULT_MISMATCH_RATE);
		String read = "CGATCGATATCGATCGATATCGATCGATATCGATCGATATCGATCGATATCGATCGATATCGATCGATATCGATCGATATCGATCGATATCGATCGATA";
		SimpleMapperResult smr = sm.map(read);
		assertEquals(SimpleMapper.HOMOLOGOUS_MAPPING, smr.getPos());
		assertEquals(0, smr.getMismatches());
	}
	
	@Test (groups = "unit" )
	public void testMapAmbiguousMatchWithMismatches() {
		SimpleMapper sm = new SimpleMapper(contig2, DEFAULT_MISMATCH_RATE);
		String read = "CGATCGATATCGATCGATAACGATCGATATCGATCGATATCGATCGATATCGATCGATATCGATCGATATCGATCGATATCGATCGATTACGATCGATA";
		SimpleMapperResult smr = sm.map(read);
		assertEquals(SimpleMapper.HOMOLOGOUS_MAPPING, smr.getPos());
		assertEquals(3, smr.getMismatches());
	}
	
	@Test (groups = "unit" )
	public void testSimple1Mismatch() {
		SimpleMapper sm = new SimpleMapper("ATCGAAAAAATTTTTTCCCCCCGGGGGGATCGGCTAATCG", DEFAULT_MISMATCH_RATE);
		String read =                          "ATAAAATTTTTTCCCCCCGGGGGGATCG";
	
		SimpleMapperResult smr = sm.map(read);
		assertEquals(4, smr.getPos());
		assertEquals(1, smr.getMismatches());
		assertEquals(Orientation.FORWARD, smr.getOrientation());
	}
	
	@Test (groups = "unit")
	public void testShortAmbiguousMatch() {
		String contig1 = "ATCGATCGATCGATCGATCGATCGATCGATCGATCG";
		String read    = "ACCGATCGATCGATCGATCGATCGATCGATCG";
		
		SimpleMapper sm = new SimpleMapper(contig1, DEFAULT_MISMATCH_RATE);
		SimpleMapperResult smr = sm.map(read);
		assertEquals(smr.getPos(), SimpleMapper.HOMOLOGOUS_MAPPING);
		assertEquals(1, smr.getMismatches());		
	}
	
	@Test (groups = "unit" )
	public void testReverseComplementExact() {
		SimpleMapper sm = new SimpleMapper(contig1, DEFAULT_MISMATCH_RATE);
		String read = "CCTGAACAGAAACTGCTTTTGGGTTTGCCAGGCACACCCACTGGAATTCAGTTAGTTTCTGAGAAGGAATGAGGGGAGTGGGCAATTCATGTTCTAGAAA";
		SimpleMapperResult smr = sm.map(read);
		assertEquals(0, smr.getMismatches());
		assertEquals(24, smr.getPos());
		assertEquals(Orientation.REVERSE, smr.getOrientation());
	}
	
	@Test (groups = "unit" )
	public void testReverseComplement2Mismatches() {
		SimpleMapper sm = new SimpleMapper(contig1, DEFAULT_MISMATCH_RATE);
		String read = "CCTGAACAGAAACTGCTTTTGGGAATGCCAGGCACACCCACTGGAATTCAGTTAGTTTCTGAGAAGGAATGAGGGGAGTGGGCAATTCATGTTCTAGAAA";
		SimpleMapperResult smr = sm.map(read);
		assertEquals(2, smr.getMismatches());
		assertEquals(24, smr.getPos());
		assertEquals(Orientation.REVERSE, smr.getOrientation());
	}
}

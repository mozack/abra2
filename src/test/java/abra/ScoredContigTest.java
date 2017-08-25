package abra;

import java.util.List;

import org.testng.annotations.Test;

import static org.testng.Assert.assertEquals;

public class ScoredContigTest {

	@Test (groups = "unit")
	public void testConvertAndFilter() {
		String contigStrings = 
				">chr7_55242125_55242525_53_-3.583754\n" +
				"contig1\n" +
				">chr7_55242125_55242525_54_-3.333877\n" +
				"contig2\n" +
				">chr7_55242125_55242525_55_-3.525407\n" +
				"contig3\n" +
				">chr7_55242125_55242525_56_-3.370506\n" +
				"contig4\n" +
				">chr7_55242125_55242525_57_-3.008790\n" +
				"contig5\n" +
				">chr7_55242125_55242525_58_-2.607389\n" +
				"contig6\n" +
				">chr7_55242125_55242525_59_-2.357512\n" +
				"contig7\n" +
				">chr7_55242125_55242525_60_-2.549043\n" +
				"contig8\n" +
				">chr7_55242125_55242525_61_-2.394141\n" +
				"contig9\n" +
				">chr7_55242125_55242525_62_-2.911880\n" +
				"contig10\n" +
				">chr7_55242125_55242525_63_-2.707760\n" +
				"contig11\n" +
				">chr7_55242125_55242525_64_-5.462782\n" +
				"contig12\n";
		
		int maxContigs = 3;
		List<ScoredContig> contigs = ScoredContig.convertAndFilter(contigStrings, maxContigs, new StringBuffer());
		assertEquals(contigs.size(), 3);
		assertEquals(contigs.get(0).getScore(), -2.357512);
		assertEquals(contigs.get(0).getContig(), "contig7");
		assertEquals(contigs.get(1).getScore(), -2.394141);
		assertEquals(contigs.get(1).getContig(), "contig9");
		assertEquals(contigs.get(2).getScore(), -2.549043);
		assertEquals(contigs.get(2).getContig(), "contig8");
	}
}

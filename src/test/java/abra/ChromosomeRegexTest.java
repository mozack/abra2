package abra;

import static org.testng.Assert.assertFalse;
import static org.testng.Assert.assertTrue;

import org.testng.annotations.Test;

public class ChromosomeRegexTest {
	
	private ChromosomeRegex cr = new ChromosomeRegex(ChromosomeRegex.DEFAULT_SKIP_REGEX);

	@Test (groups = "unit")
	public void testChr1() {
		assertFalse(cr.matches("chr1"));
	}
	
	@Test (groups = "unit")
	public void testChrX() {
		assertFalse(cr.matches("chrX"));
	}
	
	@Test (groups = "unit")
	public void testHg19Decoy() {
		assertTrue(cr.matches("hs37d5"));
	}
	
	@Test (groups = "unit")
	public void testHg19Gl() {
		assertTrue(cr.matches("GL000238.1"));
	}
	
	@Test (groups = "unit")
	public void testRandom() {
		assertTrue(cr.matches("chr1_KI270706v1_random"));
	}
	
	@Test (groups = "unit")
	public void testUnplaced() {
		assertTrue(cr.matches("chrUn_KI270320v1"));
	}
}

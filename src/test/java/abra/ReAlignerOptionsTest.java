package abra;

import static org.testng.Assert.assertEquals;

import org.testng.annotations.Test;

/**
 * Unit tests for ReAlignerOptions.
 * 
 * @author Lisle E. Mose (lmose at unc dot edu)
 */
public class ReAlignerOptionsTest {

	@Test (groups = "unit")
	public void testNoParams() {
		ReAlignerOptions options = new ReAlignerOptions();
		options.parseOptions(new String[0]);
	}
	
	@Test (groups = "unit")
	public void testBamParams() {
		ReAlignerOptions options = new ReAlignerOptions();
		options.parseOptions("--in in1.bam,in2.bam --out out1.bam,out2.bam".split("\\s"));
		String[] input = options.getInputFiles();
		String[] output = options.getOutputFiles();
		
		assertEquals(input.length, 2);
		assertEquals(input[0], "in1.bam");
		assertEquals(input[1], "in2.bam");
		assertEquals(output.length, 2);
		assertEquals(output[0], "out1.bam");
		assertEquals(output[1], "out2.bam");
	}
}

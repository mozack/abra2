package abra.cadabra;

import static org.testng.Assert.assertEquals;
import org.testng.annotations.Test;

public class HomopolymerRunTest {

	@Test (groups = "unit")
	public void testBasic() {
		String seq = "ATCGATCGTATTTTTTTTTT";
		
		HomopolymerRun hrun = HomopolymerRun.find(seq);
		assertEquals(10, hrun.getLength());
		assertEquals('T', hrun.getBase());
		assertEquals(10, hrun.getPos());
	}
	
	@Test (groups = "unit")
	public void testLeading() {
		String seq = "AAAAAAATAAAATACGACTA";
		
		HomopolymerRun hrun = HomopolymerRun.find(seq);
		assertEquals(7, hrun.getLength());
		assertEquals('A', hrun.getBase());
		assertEquals(0, hrun.getPos());
	}
	
	@Test (groups = "unit")
	public void testTrailing() {
		String seq = "AAATAAATAAAGTACGCCCC";
		
		HomopolymerRun hrun = HomopolymerRun.find(seq);
		assertEquals(4, hrun.getLength());
		assertEquals('C', hrun.getBase());
		assertEquals(16, hrun.getPos());
	}
	
	@Test (groups = "unit")
	public void testNone() {
		String seq = "AAATAAATAAAGTACGCGCC";
		
		HomopolymerRun hrun = HomopolymerRun.find(seq);
		assertEquals(null, hrun);
	}
}

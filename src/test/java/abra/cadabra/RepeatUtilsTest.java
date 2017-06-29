package abra.cadabra;

import org.testng.Assert;
import org.testng.annotations.Test;

public class RepeatUtilsTest {
	
	@Test (groups = "unit")
	public void testGetRepeatUnit_Str1() {
		String rp  = RepeatUtils.getRepeatUnit("ATATAT");
		Assert.assertEquals(rp, "AT");
	}
	
	@Test (groups = "unit")
	public void testGetRepeatUnit_Str2() {
		String rp  = RepeatUtils.getRepeatUnit("ATCATC");
		Assert.assertEquals(rp, "ATC");
	}
	
	@Test (groups = "unit")
	public void testGetRepeatUnit_Str3() {
		String rp  = RepeatUtils.getRepeatUnit("ATCGATCGATCGATCG");
		Assert.assertEquals(rp, "ATCG");
	}

	@Test (groups = "unit")
	public void testGetRepeatUnit_NearStr1() {
		String rp  = RepeatUtils.getRepeatUnit("ATATA");
		Assert.assertEquals(rp, "ATATA");
	}

	@Test (groups = "unit")
	public void testGetRepeatUnit_Homopolymer() {
		String rp = RepeatUtils.getRepeatUnit("GGG");
		Assert.assertEquals(rp, "G");
	}
	
	@Test (groups = "unit")
	public void testGetRepeatUnit_NearHomopolymer1() {
		String rp = RepeatUtils.getRepeatUnit("GGGGC");
		Assert.assertEquals(rp, "GGGGC");
	}
	
	@Test (groups = "unit")
	public void testGetRepeatUnit_NearHomopolymer2() {
		String rp = RepeatUtils.getRepeatUnit("CGGGG");
		Assert.assertEquals(rp, "CGGGG");
	}
	
	@Test (groups = "unit")
	public void testGetRepeatUnit_NearHomopolymer3() {
		String rp = RepeatUtils.getRepeatUnit("GGCGG");
		Assert.assertEquals(rp, "GGCGG");
	}
	
	@Test (groups = "unit")
	public void testGetRepeatUnit_SingleNt() {
		String rp  = RepeatUtils.getRepeatUnit("T");
		Assert.assertEquals(rp, "T");
	}
	
	@Test (groups = "unit")
	public void testGetRepeatPeriod_HpRun() {
		int period  = RepeatUtils.getRepeatPeriod("T", "TTTTTATATAT");
		Assert.assertEquals(period, 5);
	}
	
	@Test (groups = "unit")
	public void testGetRepeatPeriod_NoRepeat1() {
		int period  = RepeatUtils.getRepeatPeriod("T", "ATTTTATATAT");
		Assert.assertEquals(period, 0);
	}
	
	@Test (groups = "unit")
	public void testGetRepeatPeriod_Str1() {
		int period  = RepeatUtils.getRepeatPeriod("TA", "TATATCTA");
		Assert.assertEquals(period, 2);
	}
	
	@Test (groups = "unit")
	public void testGetRepeatPeriod_Str2() {
		int period  = RepeatUtils.getRepeatPeriod("TAC", "TACTAC");
		Assert.assertEquals(period, 2);
	}
	
	@Test (groups = "unit")
	public void testGetRepeatPeriod_Str3() {
		int period  = RepeatUtils.getRepeatPeriod("ATT", "ATTATTATTATTATTATTATTATT");
		Assert.assertEquals(period, 8);
	}
	
	@Test (groups = "unit")
	public void testGetRepeatPeriod_NoRepeat() {
		int period  = RepeatUtils.getRepeatPeriod("A", "GGGG");
		Assert.assertEquals(period, 0);
	}
}

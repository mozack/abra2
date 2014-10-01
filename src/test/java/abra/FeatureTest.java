package abra;

import org.testng.Assert;
import org.testng.annotations.Test;

public class FeatureTest {

	@Test (groups = "unit")
	public void testOverlaps() {
	
		Feature feature = new Feature("chr1", 100, 200);
		
		Assert.assertTrue(feature.overlaps("chr1", 50, 150));
		Assert.assertTrue(feature.overlaps("chr1", 1, 100));
		Assert.assertTrue(feature.overlaps("chr1", 200, 300));
		Assert.assertTrue(feature.overlaps("chr1", 150, 250));
		Assert.assertTrue(feature.overlaps("chr1", 1, 300));
		Assert.assertTrue(feature.overlaps("chr1", 100, 200));
		Assert.assertTrue(feature.overlaps("chr1", 100, 201));
		Assert.assertTrue(feature.overlaps("chr1", 99, 200));
		Assert.assertTrue(feature.overlaps("chr1", 2, 101));
		Assert.assertTrue(feature.overlaps("chr1", 199, 299));
		
		Assert.assertFalse(feature.overlaps("chr2", 50, 150));
		Assert.assertFalse(feature.overlaps("chr1", 1, 99));
		Assert.assertFalse(feature.overlaps("chr2", 101, 200));
		Assert.assertFalse(feature.overlaps("chr1", 1000, 1001));
	}
}

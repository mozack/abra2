package abra;

import java.util.Arrays;
import java.util.List;
import java.util.Map;

import static org.testng.Assert.assertEquals;
import org.testng.annotations.Test;

public class JunctionUtilsTest {
	
	private void assertJunctionEquals(Feature actual, String chr, int start, int end) {
		assertEquals(actual.getSeqname(), chr, "Junction: " + chr + ":" + start + "-" + end);
		assertEquals(actual.getStart(), start, "Junction: " + chr + ":" + start + "-" + end);
		assertEquals(actual.getEnd(), end, "Junction: " + chr + ":" + start + "-" + end);
	}

	@Test (groups = "unit")
	public void testGetRegionJunctions() throws Exception {
		RegionLoader loader = new RegionLoader();
		List<Feature> junctions = loader.load("test-data/junctions1.tab", false);
		assertEquals(junctions.size(), 29);
		
		Feature region = new Feature("chr4", 1803001, 1803401);
		List<Feature> regions = Arrays.asList(region);
		int readLength = 48;
		int maxRegionLength = 400;
		Map<Feature, List<Feature>> regionJunctionMap = JunctionUtils.getRegionJunctions(regions, junctions, readLength, maxRegionLength);
		
		List<Feature> regionJunctions = regionJunctionMap.get(region);
		
		assertEquals(regionJunctions.size(), 12);
		assertJunctionEquals(regionJunctions.get(0), "chr4", 1800000, 1801530);
		assertJunctionEquals(regionJunctions.get(1), "chr4", 1801540, 1803093);
		assertJunctionEquals(regionJunctions.get(2), "chr4", 1803264, 1803346);
		assertJunctionEquals(regionJunctions.get(3), "chr4", 1803471, 1803561);
		assertJunctionEquals(regionJunctions.get(4), "chr4", 1803471, 1803590);
		assertJunctionEquals(regionJunctions.get(5), "chr4", 1803651, 1808025);
		assertJunctionEquals(regionJunctions.get(6), "chr4", 1803714, 1805418);
		assertJunctionEquals(regionJunctions.get(7), "chr4", 1803753, 1804640);
		assertJunctionEquals(regionJunctions.get(8), "chr4", 1803753, 1805418);
		assertJunctionEquals(regionJunctions.get(9), "chr4", 1803753, 1806056);
		assertJunctionEquals(regionJunctions.get(10), "chr4", 1803753, 1806550);
		assertJunctionEquals(regionJunctions.get(11), "chr4", 1808055, 1808272);
				 
	}
}

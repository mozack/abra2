package abra;

import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import static org.testng.Assert.assertEquals;
import static org.testng.Assert.assertTrue;
import org.testng.annotations.Test;

public class JunctionUtilsTest {
	
	private void assertJunctionEquals(Feature actual, String chr, int start, int end) {
		assertEquals(actual.getSeqname(), chr, "Junction: " + chr + ":" + start + "-" + end);
		assertEquals(actual.getStart(), start, "Junction: " + chr + ":" + start + "-" + end);
		assertEquals(actual.getEnd(), end, "Junction: " + chr + ":" + start + "-" + end);
	}
	
	@Test (groups = "unit")
	public void testLoadJunctionsFromGtf() throws Exception {
		Set<Feature> junctions = JunctionUtils.loadJunctionsFromGtf("test-data/annotation1.gtf");
		assertEquals(junctions.size(), 4);
		assertTrue(junctions.contains(new Feature("chr7", 55087059, 55209978)));
		assertTrue(junctions.contains(new Feature("chr7", 55210131, 55210997)));
		assertTrue(junctions.contains(new Feature("chr7", 55211182, 55218986)));
		assertTrue(junctions.contains(new Feature("chr7", 55211182, 55214298)));
	}

	@Test (groups = "unit")
	public void testGetRegionJunctions() throws Exception {
		
		RegionLoader loader = new RegionLoader();
		List<Feature> junctions = loader.load("test-data/junctions1.tab", false);
		assertEquals(junctions.size(), 29);
		
		Feature region = new Feature("chr4", 1803001, 1803401);
		List<Feature> regions = Arrays.asList(region);
		int readLength = 100;
		int maxRegionLength = 400;
		Map<Feature, List<Feature>> regionJunctionMap = JunctionUtils.getRegionJunctions(regions, junctions, readLength, maxRegionLength);
		
		List<Feature> regionJunctions = regionJunctionMap.get(region);
		
		// TODO: Revisit this.
		//assertEquals(regionJunctions.size(), 12);
		
		assertEquals(regionJunctions.size(), 8);
		assertJunctionEquals(regionJunctions.get(0), "chr4", 1800000, 1801530);
		assertJunctionEquals(regionJunctions.get(1), "chr4", 1801251, 1801473);
		assertJunctionEquals(regionJunctions.get(2), "chr4", 1801540, 1803093);
		assertJunctionEquals(regionJunctions.get(3), "chr4", 1803264, 1803346);		
		assertJunctionEquals(regionJunctions.get(4), "chr4", 1803471, 1803561);
		assertJunctionEquals(regionJunctions.get(5), "chr4", 1803471, 1803590);
		assertJunctionEquals(regionJunctions.get(6), "chr4", 1803651, 1808025);

		//		assertJunctionEquals(regionJunctions.get(6), "chr4", 1803714, 1805418);
//		assertJunctionEquals(regionJunctions.get(7), "chr4", 1803753, 1804640);
//		assertJunctionEquals(regionJunctions.get(8), "chr4", 1803753, 1805418);
//		assertJunctionEquals(regionJunctions.get(9), "chr4", 1803753, 1806056);
//		assertJunctionEquals(regionJunctions.get(10), "chr4", 1803753, 1806550);

		assertJunctionEquals(regionJunctions.get(7), "chr4", 1808055, 1808272);
				 
		
//		List<List<Feature>> junctionPerms = JunctionUtils.combineJunctions(regionJunctions, readLength);
//		System.out.println(junctionPerms.size());
	}
	
	@Test (groups = "unit")
	public void testCombineJunctions() throws Exception {
		
		int readLength = 50;
		
		Feature j1 = new Feature("chr1", 10000, 10100);
		Feature j2 = new Feature("chr1", 10110, 10200);
		Feature j3 = new Feature("chr1", 10110, 10300);
		Feature j4 = new Feature("chr1", 10330, 10500);
		
		List<Feature> inputJunctions = Arrays.asList(j1, j2, j3, j4);
		List<List<Feature>> junctionPerms = JunctionUtils.combineJunctions(new Feature("chr1", 10000, 10400), inputJunctions, new HashSet<Feature>(), readLength, readLength);
		assertEquals(junctionPerms.size(), 8);
		
		// Expected permutations
		List<Feature> p1 = Arrays.asList(j1);
		List<Feature> p2 = Arrays.asList(j2);
		List<Feature> p3 = Arrays.asList(j3);
		List<Feature> p4 = Arrays.asList(j4);
		List<Feature> p5 = Arrays.asList(j1,j2);
		List<Feature> p6 = Arrays.asList(j1,j3);
		List<Feature> p7 = Arrays.asList(j3,j4);
		List<Feature> p8 = Arrays.asList(j1,j3,j4);
		
		assertTrue(junctionPerms.contains(p1));
		assertTrue(junctionPerms.contains(p2));
		assertTrue(junctionPerms.contains(p3));
		assertTrue(junctionPerms.contains(p4));
		assertTrue(junctionPerms.contains(p5));
		assertTrue(junctionPerms.contains(p6));
		assertTrue(junctionPerms.contains(p7));
		assertTrue(junctionPerms.contains(p8));
	}
}

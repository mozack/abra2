/* Copyright 2013 University of North Carolina at Chapel Hill.  All rights reserved. */
package abra;

import static org.testng.Assert.assertEquals;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import htsjdk.samtools.SAMRecord;

import org.testng.Assert;
import org.testng.annotations.Test;

public class RealignerTest {

	/*
	@Test (groups = "unit")
	public void testUpdateReadAlignment() {
		SAMRecord contigRead = new SAMRecord(null);
		contigRead.setReferenceName("chr21");
		contigRead.setCigarString("201M50I1105M2D563M");
		contigRead.setAlignmentStart(36205060);
		
		SAMRecord origRead = new SAMRecord(null);
		
		List<ReadBlock> blocks = ReadBlock.getReadBlocks(contigRead);
		Assert.assertEquals(blocks.size(), 5);
		
		ReadPosition readPosition = new ReadPosition(origRead, 1266, 0);
	}
	*/
	
	// TODO: Move to RegionLoaderTest
	@Test (groups = "unit")
	public void testCollapseRegions() {
		List<Feature> input = new ArrayList<Feature>();
		input.add(new Feature("chr20", 1, 2000));
		input.add(new Feature("chr20", 2050, 10000));
		input.add(new Feature("chr20", 10020, 20000));
		input.add(new Feature("chr20", 20100, 20200));
		input.add(new Feature("chr21", 20201, 20300));
		
		List<Feature> features = RegionLoader.collapseRegions(input, 70);
		assertEquals(features.size(), 3);
		validateFeature(features.get(0), 1, 20000);
		validateFeature(features.get(1), 20100, 20200);
		validateFeature(features.get(2), 20201, 20300);
	}
	
	@Test (groups = "unit")
	public void testPairJunctions() {
		Feature j1 = new Feature("chr8", 27303437, 27308265);
		Feature j2 = new Feature("chr8", 27308413, 27308559);
		Feature j3 = new Feature("chr8", 27308596, 27309001);
		Feature j4 = new Feature("chr8", 27309027, 27310630);
		Feature j5 = new Feature("chr8", 27309227, 27314630);

		List<Feature> junctions = Arrays.asList(j1, j2, j3, j4, j5);
		
		ReAligner r = new ReAligner();
		List<Pair<Feature, Feature>> junctionPairs = r.pairJunctions(junctions, 50);
		assertEquals(junctionPairs.size(), 2);
		
		// 1st pair
		assertEquals(j2, junctionPairs.get(0).getFirst());
		assertEquals(j3, junctionPairs.get(0).getSecond());
		
		// 2nd pair
		assertEquals(j3, junctionPairs.get(1).getFirst());
		assertEquals(j4, junctionPairs.get(1).getSecond());
	}
	
	@Test (groups = "unit")
	public void testPairJunctions_cannotAppearInSameContig() {
		Feature j1 = new Feature("chr1", 755395, 755674);
		Feature j2 = new Feature("chr1", 755430, 755674);
		Feature j3 = new Feature("chr1", 755535, 755674);
		
		List<Feature> junctions = Arrays.asList(j1, j2, j3);
		
		ReAligner r = new ReAligner();
		List<Pair<Feature, Feature>> junctionPairs = r.pairJunctions(junctions, 5000);
		assertEquals(junctionPairs.size(), 0);
	}
	
	/*
	@Test (groups = "unit")
	public void testSplitRegions() {
		List<Feature> input = new ArrayList<Feature>();
		input.add(new Feature("chr20", 1, 2000));
		input.add(new Feature("chr20", 2001, 10000));
		input.add(new Feature("chr20", 9523233, 9523367));
		input.add(new Feature("chr21", 36205258, 36206898));
		input.add(new Feature("chr22", 21303999, 21308037));
		
		ReAligner realigner = new ReAligner();
		List<Feature> features = realigner.splitRegions(input, 2000, 500, 200);
		
		assertEquals(features.size(), 9);
		// first
		validateFeature(features.get(0), 1, 2000);
		// second
		validateFeature(features.get(1), 2001, 4201);
		validateFeature(features.get(2), 4001, 6201);
		validateFeature(features.get(3), 6001, 8201);
		validateFeature(features.get(4), 8001, 10000);
		// third
		validateFeature(features.get(5), 9523233, 9523367);
		// fourth
		validateFeature(features.get(6), 36205258, 36206898);
		// fifth
		validateFeature(features.get(7), 21303999, 21306199);
		validateFeature(features.get(8), 21305999, 21308037);		
	}
	*/
	
	private void validateFeature(Feature feature, int expectedStart, int expectedEnd) {
		assertEquals(feature.getStart(), expectedStart);
		assertEquals(feature.getEnd(), expectedEnd);
	}
}

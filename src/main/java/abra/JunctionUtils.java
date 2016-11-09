package abra;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

public class JunctionUtils {
	
	/**
	 * Return Map with key=Region, value = Sorted list of junctions that may be relevant to the region
	 */
	static Map<Feature, List<Feature>> getRegionJunctions(List<Feature> chromosomeRegions, List<Feature> chromosomeJunctions,
			int readLength, int maxRegionLength) {
		System.err.println("Assigning junctions to regions");
		// TODO: Brute force matching of regions / junctions
		// TODO: Match up more efficiently
		
		// key = region, value = junction list
		Map<Feature, List<Feature>> regionJunctions = new HashMap<Feature, List<Feature>>();
		
		Map<Integer, List<Feature>> chromosomeJunctionsByStart = new HashMap<Integer, List<Feature>>();
		for (Feature junction : chromosomeJunctions) {
			addToChromosomePositionMap(junction, (int) junction.getStart(), chromosomeJunctionsByStart);
		}

		Map<Integer, List<Feature>> chromosomeJunctionsByEnd = new HashMap<Integer, List<Feature>>();
		for (Feature junction : chromosomeJunctions) {
			addToChromosomePositionMap(junction, (int) junction.getEnd(), chromosomeJunctionsByEnd);
		}
		
		for (Feature region : chromosomeRegions) {
		
			// Junctions for current region
			Set<Feature> localJunctions = new HashSet<Feature>();;
			
			for (int pos=(int) region.getStart()-maxRegionLength; pos<region.getEnd()+maxRegionLength; pos++) {
				if (chromosomeJunctionsByStart.containsKey(pos)) {
					localJunctions.addAll(chromosomeJunctionsByStart.get(pos));
				}
				
				if (chromosomeJunctionsByEnd.containsKey(pos)) {
					localJunctions.addAll(chromosomeJunctionsByEnd.get(pos));
				}
			}
			
			// Add neighboring junctions (up to 2 additional splices)
			addNeighboringJunctions(localJunctions, chromosomeJunctionsByStart, chromosomeJunctionsByEnd, readLength);			
			addNeighboringJunctions(localJunctions, chromosomeJunctionsByStart, chromosomeJunctionsByEnd, readLength);
			
			List<Feature> localJunctionList = new ArrayList<Feature>(localJunctions);
			Collections.sort(localJunctionList, new JunctionComparator());
			
			regionJunctions.put(region, localJunctionList);
		}
		
		System.err.println("Done assigning junctions to regions");
		
		return regionJunctions;
	}
	
	private static void addToChromosomePositionMap(Feature junction, int position, Map<Integer, List<Feature>> chromosomeJunctionsByPosition) {
		if (!chromosomeJunctionsByPosition.containsKey(position)) {
			chromosomeJunctionsByPosition.put(position, new ArrayList<Feature>());
		}
		
		chromosomeJunctionsByPosition.get(position).add(junction);
	}

	// Given the set of current junctions, add any other junctions that may be within a read length distance
	private static void addNeighboringJunctions(Set<Feature> currJunctions, Map<Integer, List<Feature>> chromosomeJunctionsByStart,
			Map<Integer, List<Feature>> chromosomeJunctionsByEnd, int readLength) {
		List<Feature> toAdd = new ArrayList<Feature>();
		
		for (Feature junction : currJunctions) {
			// Look for junctions with endpoint within read length of current junction start
			for (int i=0; i<readLength; i++) {
				int idx = (int) junction.getStart() - i;
				if (chromosomeJunctionsByEnd.containsKey(idx)) {
					toAdd.addAll(chromosomeJunctionsByEnd.get(idx));
				}
			}
			
			// Look for junctions with start within read length of current junction end
			for (int i=0; i<readLength; i++) {
				int idx = (int) junction.getEnd() + i;
				if (chromosomeJunctionsByStart.containsKey(idx)) {
					toAdd.addAll(chromosomeJunctionsByStart.get(idx));
				}
			}
		}
		
		currJunctions.addAll(toAdd);
	}
	
	public static List<List<Feature>> combineJunctions(List<Feature> junctions, int readLength) {
		List<List<Feature>> combinedJunctions = new ArrayList<List<Feature>>();
		
		// Get all possible permutations of junctions regardless of validity
		List<List<Feature>> junctionLists = combineAllJunctions(junctions, readLength);
		
		for (List<Feature> currJunctions : junctionLists) {
			if (isJunctionCombinationValid(currJunctions, readLength)) {
				combinedJunctions.add(currJunctions);
			}
		}
		
		return combinedJunctions;
	}
	
	// Produce all possible junction permutations from the input list.
	private static List<List<Feature>> combineAllJunctions(List<Feature> junctions, int readLength) {
		List<List<Feature>> junctionLists = null;
		
		if (junctions.size() == 1) {
			junctionLists = Arrays.asList((List<Feature>) new ArrayList<Feature>(), (List<Feature>) new ArrayList<Feature>());
			// Return 2 lists, one with the junction and one without.
			junctionLists.get(1).add(junctions.get(0));
		} else if (junctions.size() > 1) {
			junctionLists = new ArrayList<List<Feature>>();
			Feature currentJunction = junctions.get(0);
			List<List<Feature>> subJuncs = combineAllJunctions(junctions.subList(1, junctions.size()), readLength);
			// For each returned list, create a new list with and without the current junction
			for (List<Feature> subJuncList : subJuncs) {
				// Pass along sub list without current junction
				junctionLists.add(subJuncList);
				List<Feature> newList = new ArrayList<Feature>();
				// Add new sublist with current junction
				newList.add(currentJunction);
				newList.addAll(subJuncList);
				
				if (isJunctionCombinationValid(newList, readLength)) {
					junctionLists.add(newList);
				}
			}
		} else {
			junctionLists = new ArrayList<List<Feature>>();
		}
		
		return junctionLists;
	}
	
	// Assuming all inputs on same chromosome
	protected static boolean isJunctionCombinationValid(List<Feature> junctions, int readLength) {
//		if (1==1) {
//			return true;
//		}
		
		for (int i=0; i<junctions.size()-1; i++) {
			
			// End of left junction must be less than start of right junction 
			if (junctions.get(i).getEnd() >= junctions.get(i+1).getStart()) {
				return false;
			}
			
			// Distance between junctions must be less than readLength*2
			if (junctions.get(i+1).getStart() - junctions.get(i).getEnd() > readLength*2) {
				return false;
			}
		}
		
		return junctions.size() > 0;
	}

	// Sort strictly based upon start and end pos.  Chromosome ignored.
	static class JunctionComparator implements Comparator<Feature> {

		@Override
		public int compare(Feature j1, Feature j2) {
			int ret = 0;
			
			if (j1.getStart() < j2.getStart()) {
				ret = -1;
			} else if (j1.getStart() > j2.getStart()) {
				ret = 1;
			} else if (j1.getEnd() < j2.getEnd()) {
				ret = -1;
			} else if (j1.getEnd() > j2.getEnd()) {
				ret = 1;
			}
			
			return ret;
		}
	}
}

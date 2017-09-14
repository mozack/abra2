package abra;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
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

	public static final int MAX_JUNCTION_PERMUTATIONS = 1024;
	public static final int MAX_POTENTIAL_PERMUTATIONS = 65536;
	/**
	 * Load junctions from GTF using exons grouped by transcript_id
	 * Sort order is unspecified 
	 */
	public static Set<Feature> loadJunctionsFromGtf(String filename) throws FileNotFoundException, IOException {
		Logger.debug("Loading annotated junctions from %s", filename);
		Set<Exon> exonSet = new HashSet<Exon>();
		
		BufferedReader reader = new BufferedReader(new FileReader(filename));
		
		try {
			String line = reader.readLine();
			while (line != null) {
				if (!line.startsWith("#")) {
					String[] fields = line.split("\\t");
					if (fields.length >= 9 && fields[2].equals("exon")) {
						String chrom = fields[0];
						int start = Integer.parseInt(fields[3]);
						int stop = Integer.parseInt(fields[4]);
						String attributes = fields[8];
						String[] attrFields = attributes.split(";");
						for (String attr : attrFields) {
							attr = attr.trim();
							if (attr.startsWith("transcript_id")) {
								int idx = attr.indexOf("transcript_id") + "transcript_id".length();
								String transcriptId = attr.substring(idx, attr.length());
								exonSet.add(new Exon(chrom, start, stop, transcriptId));
							}
						}
					}
				}
				
				line = reader.readLine();
			}
		} finally {
			reader.close();
		}
		
		List<Exon> exons = new ArrayList<Exon>(exonSet);
		Collections.sort(exons);
		
		Set<Feature> junctions = new HashSet<Feature>();
		Exon prevExon = null;
		
		for (Exon exon : exons) {
			if (prevExon != null && exon.getTranscriptId().equals(prevExon.getTranscriptId())) {
				// Create junction, adjusting coordinates to match first / last position in intron
				// similar to STAR's junction output
				Feature junction = new Feature(exon.getChrom(), prevExon.getStop()+1, exon.getStart()-1);
				junctions.add(junction);
			}
			
			prevExon = exon;
		}
		
		Logger.info("Loaded " + junctions.size() + " annotated junctions");
		
		return junctions;
	}
	
	/**
	 * Return Map with key=Region, value = Sorted list of junctions that may be relevant to the region
	 */
	static Map<Feature, List<Feature>> getRegionJunctions(List<Feature> chromosomeRegions, List<Feature> chromosomeJunctions,
			int readLength, int maxRegionLength) {
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
			
			int regionOverlap = readLength;
//			for (int pos=(int) region.getStart()-maxRegionLength; pos<region.getEnd()+maxRegionLength; pos++) {
//			for (int pos=(int) region.getStart()-readLength; pos<region.getEnd()+readLength; pos++) {
			for (int pos=(int) region.getStart()-regionOverlap; pos<region.getEnd()+regionOverlap; pos++) {
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
	
	public static List<List<Feature>> combineJunctions(Feature region, List<Feature> junctions, Set<Feature> preferredJunctions, int maxJuncDist, int readLength) throws TooManyJunctionPermutationsException {
		List<List<Feature>> combinedJunctions = new ArrayList<List<Feature>>();
		
		// Get all possible permutations of junctions regardless of validity
		List<List<Feature>> junctionLists = combineAllJunctions(region, junctions, maxJuncDist, readLength);
		
		for (List<Feature> currJunctions : junctionLists) {
			if (isJunctionCombinationValid(region, currJunctions, maxJuncDist, readLength) &&
				containsInRegionJunction(region, currJunctions, preferredJunctions, readLength)) {
				combinedJunctions.add(currJunctions);
			}
		}
		
		return combinedJunctions;
	}
	
	// Produce all possible junction permutations from the input list.
	private static List<List<Feature>> combineAllJunctions(Feature region, List<Feature> junctions, int maxJuncDist, int readLength) throws TooManyJunctionPermutationsException {
		List<List<Feature>> junctionLists = null;
		
		if (junctions.size() == 1) {
			junctionLists = Arrays.asList((List<Feature>) new ArrayList<Feature>(), (List<Feature>) new ArrayList<Feature>());
			// Return 2 lists, one with the junction and one without.
			junctionLists.get(1).add(junctions.get(0));
		} else if (junctions.size() > 1) {
			junctionLists = new ArrayList<List<Feature>>();
			Feature currentJunction = junctions.get(0);
			List<List<Feature>> subJuncs = combineAllJunctions(region, junctions.subList(1, junctions.size()), maxJuncDist, readLength);
			// For each returned list, create a new list with and without the current junction
			for (List<Feature> subJuncList : subJuncs) {
				// Pass along sub list without current junction
				junctionLists.add(subJuncList);
				List<Feature> newList = new ArrayList<Feature>();
				// Add new sublist with current junction
				newList.add(currentJunction);
				newList.addAll(subJuncList);
				
				if (isJunctionCombinationValid(region, newList, maxJuncDist, readLength)) {
					junctionLists.add(newList);
				}
				
				if (junctionLists.size() > MAX_POTENTIAL_PERMUTATIONS) {
					throw new TooManyJunctionPermutationsException();
				}
			}
		} else {
			junctionLists = new ArrayList<List<Feature>>();
		}
		
		return junctionLists;
	}
	
	private static boolean containsInRegionJunction(Feature region, List<Feature> junctions, Set<Feature> preferredJunctions, int readLength) {
		
		// Require at least one junction endpoint to be within the region.
		long start = region.getStart() - readLength;
		long end = region.getEnd() + readLength;

		for (Feature junction : junctions) {
			
			if (preferredJunctions.contains(junction)) {
				return true;
			}
		
			if (junction.getStart() > start && junction.getStart() < end) {
				return true;
			}
			
			if (junction.getEnd() > start && junction.getEnd() < end) {
				return true;
			}
		}
		
		return false;
	}
	
	// Assuming all inputs on same chromosome
	protected static boolean isJunctionCombinationValid(Feature region, List<Feature> junctions, int maxRegionSize, int readLength) {

		Feature paddedRegion = new Feature(region.getSeqname(), region.getStart()-maxRegionSize/2, region.getEnd()+maxRegionSize/2);
		
		for (int i=0; i<junctions.size()-1; i++) {
			
			int maxDist = maxRegionSize;
			
			Feature left = junctions.get(i);
			Feature right = junctions.get(i+1);
			
			// End of left junction must be less than start of right junction 
			if (left.getEnd() >= right.getStart()) {
				return false;
			}
			
			// Include all junctions within neighboring regions (padded region)
			// Junctions outside of padded regions should be within a read length of
			// another junction start/end point.
			
			if (!left.overlaps(paddedRegion) || !right.overlaps(paddedRegion)) {
				// If one of left or right junction is out of region,
				// limit junc distance to half region size

				// TODO: Should this be the window overlap size?
//				maxDist = maxRegionSize / 2;
				maxDist = readLength;
			}
			
			// Distance between junctions must be less than maxDist
			if (right.getStart() - left.getEnd() > maxDist) {				
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
	
	static class Exon implements Comparable<Exon> {
		int start;
		int stop;
		String chrom;
		String transcriptId;
		
		Exon(String chrom, int start, int stop, String transcriptId) {
			this.start = start;
			this.stop = stop;
			this.chrom = chrom;
			this.transcriptId = transcriptId; 
		}
		
		public int getStart() {
			return start;
		}

		public int getStop() {
			return stop;
		}

		public String getChrom() {
			return chrom;
		}

		public String getTranscriptId() {
			return transcriptId;
		}

		@Override
		public int compareTo(Exon that) {
			int cmp = this.transcriptId.compareTo(that.transcriptId);
			if (cmp == 0) {
				cmp = this.start - that.start;
			}
			
			return cmp;
		}

		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result + ((chrom == null) ? 0 : chrom.hashCode());
			result = prime * result + start;
			result = prime * result + stop;
			result = prime * result + ((transcriptId == null) ? 0 : transcriptId.hashCode());
			return result;
		}

		@Override
		public boolean equals(Object obj) {
			if (this == obj)
				return true;
			if (obj == null)
				return false;
			if (getClass() != obj.getClass())
				return false;
			Exon other = (Exon) obj;
			if (chrom == null) {
				if (other.chrom != null)
					return false;
			} else if (!chrom.equals(other.chrom))
				return false;
			if (start != other.start)
				return false;
			if (stop != other.stop)
				return false;
			if (transcriptId == null) {
				if (other.transcriptId != null)
					return false;
			} else if (!transcriptId.equals(other.transcriptId))
				return false;
			return true;
		}
	}
	
	static class TooManyJunctionPermutationsException extends Exception {
		
	}
}

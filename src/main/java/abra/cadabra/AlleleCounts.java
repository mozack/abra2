package abra.cadabra;

import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import htsjdk.samtools.SAMRecord;

public class AlleleCounts {
	
	private int count;
	private int totalCount;  // Includes overlapping pairs
	private int fwd;
	private int rev;
	private int minReadIdx = Integer.MAX_VALUE;
	private int maxReadIdx = Integer.MIN_VALUE;
	private Map<String, Integer> insertBaseCounts = new HashMap<String, Integer>();
	private Set<String> readIds = new HashSet<String>();
	
	
	public static final AlleleCounts EMPTY_COUNTS;
	
	static {
		EMPTY_COUNTS = new AlleleCounts();
		EMPTY_COUNTS.count = 0;
		EMPTY_COUNTS.fwd = 0;
		EMPTY_COUNTS.rev = 0;
		EMPTY_COUNTS.minReadIdx = 0;
		EMPTY_COUNTS.maxReadIdx = 0;
		EMPTY_COUNTS.totalCount = 0;
	}
	
	public int getCount() {
		return count;
	}
	
	public int getFwd() {
		return fwd;
	}
	
	public int getRev() {
		return rev;
	}
	
	public int getMinReadIdx() {
		return minReadIdx;
	}
	
	public int getMaxReadIdx() {
		return maxReadIdx;
	}
	
	public void incrementCount(SAMRecord read) {
		if (!readIds.contains(read.getReadName())) {
			// Don't allow multiple ends of fragment to be double counted
			count += 1;
		}
		
		totalCount += 1;
		
		if (read.getReadNegativeStrandFlag()) {
			incrementRev();
		} else {
			incrementFwd();
		}
		
		readIds.add(read.getReadName());
	}
	
	public int getTotalCount() {
		return totalCount;
	}

	private void incrementFwd() {
		fwd += 1;
	}
	
	private void incrementRev() {
		rev += 1;
	}
	
	public void clearReadIds() {
		readIds.clear();
	}
	
	public void updateReadIdx(int idx) {
		if (idx < minReadIdx) {
			minReadIdx = idx;
		}
		
		if (idx > maxReadIdx) {
			maxReadIdx = idx;
		}
	}
	
	public void updateInsertBases(String bases) {
		if (bases != null) {
			if (insertBaseCounts.containsKey(bases)) {
				insertBaseCounts.put(bases, insertBaseCounts.get(bases)+1);
			} else {
				insertBaseCounts.put(bases, 1);
			}
		}
	}
	
	public String getPreferredInsertBases() {
		int max = 0;
		String maxBases = "";
		for (String bases : insertBaseCounts.keySet()) {
			int count = insertBaseCounts.get(bases);
			if (count > max) {
				max = count;
				maxBases = bases;
			}
		}
		
		return maxBases;
	}
	
	public static int sum(Collection<AlleleCounts> counts) {
		int sum = 0;
		
		for (AlleleCounts ac : counts) {
			sum += ac.getCount();
		}
		
		return sum;
	}
	
	public String toString() {
		return "count: " + count;
	}
}

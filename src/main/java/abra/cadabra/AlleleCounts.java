package abra.cadabra;

import java.util.HashMap;
import java.util.Map;

public class AlleleCounts {
	
	private int count;
	private int fwd;
	private int rev;
	private int minReadIdx = Integer.MAX_VALUE;
	private int maxReadIdx = Integer.MIN_VALUE;
	private Map<String, Integer> insertBaseCounts = new HashMap<String, Integer>();
	
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
	
	public void incrementCount() {
		count += 1;
	}
	
	public void incrementFwd() {
		fwd += 1;
	}
	
	public void incrementRev() {
		rev += 1;
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
	
	public String toString() {
		return "count: " + count;
	}
}

package abra.cadabra;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
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
	private Map<Integer, List<String>> alignmentEnds = new HashMap<Integer, List<String>>();
	private int spanEnd = -1;

	private int spanningCount = -1;
	
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
	
	public void setSpanEnd(int spanEnd) {
		this.spanEnd = spanEnd;
	}
	
	public int getCount() {
		if (spanEnd <= 0) {
			return count;
		}
		
		if (spanningCount < 0) {
			// Count number of distinct reads with an alignment end that reaches the spanEnd.
			// Distinct read ids are used to not double count overlapping fragments
			Set<String> readIds = new HashSet<String>();
			
			for (Integer alignmentEnd : alignmentEnds.keySet()) {
				if (alignmentEnd > spanEnd) {
					readIds.addAll(alignmentEnds.get(alignmentEnd));
				}
			}
			
			spanningCount = readIds.size();
		}
		
		return spanningCount;
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
		if (!alignmentEnds.containsKey(read.getAlignmentEnd())) {
			alignmentEnds.put(read.getAlignmentEnd(), new ArrayList<String>());
		}
		
		alignmentEnds.get(read.getAlignmentEnd()).add(read.getReadName());
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
	
	public static void setSpanEnd(int spanEnd, Map<Allele, AlleleCounts> counts) {
		for (Allele allele : counts.keySet()) {
			if (allele.getType() != Allele.Type.DEL && allele.getType() != Allele.Type.INS) {
				AlleleCounts ac = counts.get(allele);
				ac.setSpanEnd(spanEnd);
			}
		}
	}
	
	public String toString() {
		return "count: " + count;
	}
}

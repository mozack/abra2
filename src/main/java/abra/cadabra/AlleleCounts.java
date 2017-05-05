package abra.cadabra;

public class AlleleCounts {
	
	private int count;
	private int fwd;
	private int rev;
	private int minReadIdx = Integer.MAX_VALUE;
	private int maxReadIdx = Integer.MIN_VALUE;
	
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
}

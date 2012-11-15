package abra;

//TODO: Get rid of this class
public class Counts implements Cloneable {

	private boolean isTerminatedAtRepeat = false;
	
	public void setTerminatedAtRepeat(boolean isTerminatedAtRepeat) {
		this.isTerminatedAtRepeat = isTerminatedAtRepeat;
	}
	
	public boolean isTerminatedAtRepeat() {
		return isTerminatedAtRepeat;
	}
	
	public Object clone() {
		Counts clone = new Counts();
//		clone.edgeCounts.addAll(this.edgeCounts);
//		clone.totalEdgeCount = this.totalEdgeCount;
		
		return clone;
	}
	
	public String toString() {
		/*
		return 
			"_numedges:" + getNumEdges() + 
			"_totaledgecounts:" + getTotalEdgeCount() +
			"_medianedgecount:" + getMedianEdgeCount() +
			"_minedgecount:" + getMinEdgeCount() + 
			"_terminatedatrepeat:" + isTerminatedAtRepeat;
		*/
		return "_terminatedatrepeat:" + isTerminatedAtRepeat;
	}
}
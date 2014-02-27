package abra;


public class BreakpointCandidate {
	private Feature region;
	private int spanningReadPairCount;
	
	public BreakpointCandidate(Feature region, int spanningReadPairCount) {
		this.region = region;
		this.spanningReadPairCount = spanningReadPairCount;
	}
	
	public Feature getRegion() {
		return region;
	}
	
	public int getSpanningReadPairCount() {
		return spanningReadPairCount;
	}
}

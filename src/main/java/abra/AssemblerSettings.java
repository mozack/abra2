package abra;

public class AssemblerSettings {

	private int kmerSize;
	private int minEdgeFrequency;
	private int minNodeFrequncy;
	private int minContigLength;
	private double minEdgeRatio;
	private int maxPotentialContigs;
	private double minContigRatio;
	private int minUniqueReads;
		
	public int getKmerSize() {
		return kmerSize;
	}
	
	public void setKmerSize(int kmerSize) {
		this.kmerSize = kmerSize;
	}
	
	public int getMinEdgeFrequency() {
		return minEdgeFrequency;
	}
	
	public void setMinEdgeFrequency(int minEdgeFrequency) {
		this.minEdgeFrequency = minEdgeFrequency;
	}
	
	public int getMinNodeFrequncy() {
		return minNodeFrequncy;
	}
	
	public void setMinNodeFrequncy(int minNodeFrequncy) {
		this.minNodeFrequncy = minNodeFrequncy;
	}
	
	public int getMinContigLength() {
		return minContigLength;
	}
	
	public void setMinContigLength(int minContigLength) {
		this.minContigLength = minContigLength;
	}
	
	public double getMinEdgeRatio() {
		return minEdgeRatio;
	}
	
	public void setMinEdgeRatio(double minEdgeRatio) {
		this.minEdgeRatio = minEdgeRatio;
	}
	
	public int getMaxPotentialContigs() {
		return maxPotentialContigs;
	}
	
	public void setMaxPotentialContigs(int maxPotentialContigs) {
		this.maxPotentialContigs = maxPotentialContigs;
	}
	
	public double getMinContigRatio() {
		return minContigRatio;
	}
	
	public void setMinContigRatio(double minContigRatio) {
		this.minContigRatio = minContigRatio;
	}
	
	public int getMinUniqueReads() {
		return minUniqueReads;
	}
	
	public void setMinUniqueReads(int minUniqueReads) {
		this.minUniqueReads = minUniqueReads;
	}
	
	public String getDescription() {
		StringBuffer str = new StringBuffer();
		
		appendSetting(str, "kmerSize", kmerSize);
		appendSetting(str, "minEdgeFrequency", minEdgeFrequency);
		appendSetting(str, "minNodeFrequncy", minNodeFrequncy);
		appendSetting(str, "minContigLength", minContigLength);
		appendSetting(str, "minEdgeRatio", minEdgeRatio);
		appendSetting(str, "maxPotentialContigs", maxPotentialContigs);
		appendSetting(str, "minContigRatio", minContigRatio);
		appendSetting(str, "minUniqueReads", minUniqueReads);
		
		return str.toString();
	}
	
	private void appendSetting(StringBuffer str, String setting, int value) {
		str.append(setting);
		str.append(": ");
		str.append(value);
		str.append('\n');
	}
	
	private void appendSetting(StringBuffer str, String setting, double value) {
		str.append(setting);
		str.append(": ");
		str.append(value);
		str.append('\n');
	}
}

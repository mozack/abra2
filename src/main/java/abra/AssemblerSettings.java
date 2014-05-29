/* Copyright 2013 University of North Carolina at Chapel Hill.  All rights reserved. */
package abra;

/**
 * Settings for Contig Assembly
 * @author Lisle E. Mose (lmose at unc dot edu)
 */
public class AssemblerSettings {

	private int[] kmerSize;
	private int minEdgeFrequency;
	private int minNodeFrequncy;
	private int minUnalignedNodeFrequency;
	private int minContigLength;
	private int maxPotentialContigs;
	private int minBaseQuality;
	private double minReadCandidateFraction;
	private int maxAverageDepth;
	private boolean searchForStructuralVariation;
	private int averageDepthCeiling;
		
	public int getAverageDepthCeiling() {
		return averageDepthCeiling;
	}

	public void setAverageDepthCeiling(int averageDepthCeiling) {
		this.averageDepthCeiling = averageDepthCeiling;
	}

	public boolean searchForStructuralVariation() {
		return searchForStructuralVariation;
	}

	public void setSearchForStructuralVariation(boolean searchForStructuralVariation) {
		this.searchForStructuralVariation = searchForStructuralVariation;
	}

	public int getMaxAverageDepth() {
		return maxAverageDepth;
	}

	public void setMaxAverageDepth(int maxAverageDepth) {
		this.maxAverageDepth = maxAverageDepth;
	}

	public double getMinReadCandidateFraction() {
		return minReadCandidateFraction;
	}

	public void setMinReadCandidateFraction(double minReadCandidateFraction) {
		this.minReadCandidateFraction = minReadCandidateFraction;
	}

	public int getMinBaseQuality() {
		return minBaseQuality;
	}

	public void setMinBaseQuality(int minBaseQuality) {
		this.minBaseQuality = minBaseQuality;
	}

	public int[] getKmerSize() {
		return kmerSize;
	}
	
	public void setKmerSize(int[] kmerSize) {
		this.kmerSize = kmerSize;
	}
	
	public int getMinEdgeFrequency() {
		return minEdgeFrequency;
	}
	
	public void setMinEdgeFrequency(int minEdgeFrequency) {
		this.minEdgeFrequency = minEdgeFrequency;
	}
	
	public void setMinUnalignedNodeFrequency(int minUnalignedNodeFrequency) {
		this.minUnalignedNodeFrequency = minUnalignedNodeFrequency;
	}
	
	public int getMinUnalignedNodeFrequency() {
		return minUnalignedNodeFrequency;
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
	
	public int getMaxPotentialContigs() {
		return maxPotentialContigs;
	}
	
	public void setMaxPotentialContigs(int maxPotentialContigs) {
		this.maxPotentialContigs = maxPotentialContigs;
	}
		
	public String getDescription() {
		StringBuffer str = new StringBuffer();
		
		for (int i=0; i<kmerSize.length; i++) {
			appendSetting(str, "kmer" + i, kmerSize[i]);
		}
		appendSetting(str, "minEdgeFrequency", minEdgeFrequency);
		appendSetting(str, "minNodeFrequncy", minNodeFrequncy);
		appendSetting(str, "minContigLength", minContigLength);
		appendSetting(str, "maxPotentialContigs", maxPotentialContigs);
		appendSetting(str, "minBaseQuality", minBaseQuality);
		appendSetting(str, "minReadCandidateFraction", minReadCandidateFraction);
		appendSetting(str, "maxAverageRegionDepth", maxAverageDepth);
		
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

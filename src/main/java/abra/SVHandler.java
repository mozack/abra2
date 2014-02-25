package abra;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMRecord;

/**
 * Handles processing of SV candidates.<br/>
 * This class is not threadsafe.
 * 
 * @author Lisle E. Mose (lmose at unc dot edu)
 */
public class SVHandler {
	
	private static final int MIN_MAPQ = 45;
	
	private int readLength;
	//TODO: Re-evaluate this range and consider not hardcoding.
	private static final int GROUP_RANGE = 1000;
	
	public SVHandler(int readLength) {
		this.readLength = readLength;
	}

	public boolean identifySVCandidates(String input, String output) throws IOException {
		boolean hasCandidates = false;
		SamMultiMappingReader reader = new SamMultiMappingReader(input);
		
		BufferedWriter writer = new BufferedWriter(new FileWriter(output, false));

		List<Breakpoint> breakpoints = new ArrayList<Breakpoint>();
		for (List<SAMRecord> readList : reader) {
			Breakpoint breakpoint = processRead(readList);
			
			if (breakpoint != null) {
//				System.out.println(breakpoint);
				breakpoints.add(breakpoint);
			}
		}
		
		List<BreakpointGroup> groups = getBreakpointGroups(breakpoints);
		
		for (BreakpointGroup group : groups) {
			for (Breakpoint breakpoint : group.getBreakpoints()) {
				String desc = ">BP_" + group.getGroupId() + "_" + breakpoint.getLabel();
				writer.write(desc + "\n");
				writer.write(breakpoint.getBases() + "\n");
//				System.out.println(desc);
//				System.out.println(breakpoint.getBases());
				hasCandidates = true;
			}
		}
		
		reader.close();
		writer.close();
		
		return hasCandidates;
	}
	
	private List<BreakpointGroup> getBreakpointGroups(List<Breakpoint> breakpoints) {
		Collections.sort(breakpoints);
		List<BreakpointGroup> cached = new ArrayList<BreakpointGroup>();
		List<BreakpointGroup> groups = new ArrayList<BreakpointGroup>();
		
		for (Breakpoint breakpoint : breakpoints) {
			List<BreakpointGroup> newCache = new ArrayList<BreakpointGroup>();
			
			boolean isMerged = false;
			for (BreakpointGroup group : cached) {
				if (group.leftChr.equals(breakpoint.leftChr) && group.leftStart > breakpoint.leftPos-GROUP_RANGE) {
					newCache.add(group);
				}

				if (!isMerged && group.shouldMerge(breakpoint)) {
					group.addBreakpoint(breakpoint);
					isMerged = true;
				}
			}
			
			if (!isMerged) {
				BreakpointGroup group = new BreakpointGroup(breakpoint);
				groups.add(group);
				newCache.add(group);
			}
			
			cached = newCache;
		}
		
		return groups;
	}
	
	protected Breakpoint processRead(List<SAMRecord> readList) {
		if (readList.size() != 2) {
			return null;
		}
		
		SAMRecord read1 = readList.get(0);
		SAMRecord read2 = readList.get(1);
		
		if (read1.getMappingQuality() < MIN_MAPQ || read2.getMappingQuality() < MIN_MAPQ) {
			return null;
		}
		
		SAMRecord primary = null;
		SAMRecord secondary = null;
		
		if (SAMRecordUtils.isPrimary(read1) && !SAMRecordUtils.isPrimary(read2)) {
			primary = read1;
			secondary = read2;
		} else if (!SAMRecordUtils.isPrimary(read1) && SAMRecordUtils.isPrimary(read2)) {
			primary = read2;
			secondary = read1;
		} else {
			return null;
		}
		
		if (!isInTargetRegions(primary, secondary)) {
			return null;
		}
		
//		System.out.println("-------------------------------------------");
//		System.out.println(primary.getSAMString());
//		System.out.println(secondary.getSAMString());
		
		return getBreakpoint(primary, secondary);		
	}
	
	private boolean isInTargetRegions(SAMRecord primary, SAMRecord secondary) {
		boolean isInTargetRegions = false;
		
		//TODO: Handle underscore in chromosome name.
		String[] fields = primary.getReadName().split("_");
		String chr1 = fields[0];
		int start1 = Integer.parseInt(fields[1]);
		int stop1 = Integer.parseInt(fields[2]);
		
		String chr2 = fields[4];
		int start2 = Integer.parseInt(fields[5]);
		int stop2 = Integer.parseInt(fields[6]);
		
		Feature region1 = new Feature(chr1, start1, stop1);
		Feature region2 = new Feature(chr2, start2, stop2);
		
		if (region1.overlapsRead(primary) || region1.overlapsRead(secondary)) {
			if (region2.overlapsRead(primary) || region2.overlapsRead(secondary)) {
				isInTargetRegions = true;
			}
		}
		
		return isInTargetRegions;
	}
	
	private Breakpoint getBreakpoint(SAMRecord primary, SAMRecord secondary) {
//		System.out.println("pl: " + primary.getReadLength());
//		System.out.println("sl: " + secondary.getReadLength());
		
		int readIdx;
		String left;
		String right;
		
		String leftChr;
		String rightChr;
		int leftPos;
		int rightPos;
		
		if (secondary.getCigar().getCigarElement(0).getOperator() == CigarOperator.HARD_CLIP) {
			readIdx = secondary.getCigar().getCigarElement(0).getLength();
			left = primary.getReferenceName() + "_" + (primary.getAlignmentStart() + readIdx) + "_" + (primary.getReadNegativeStrandFlag() ? "R" : "F");
			right = secondary.getReferenceName() + "_" + secondary.getAlignmentStart() + "_" + (secondary.getReadNegativeStrandFlag() ? "R" : "F");
			
			leftChr = primary.getReferenceName();
			leftPos = primary.getAlignmentStart() + readIdx;
			rightChr = secondary.getReferenceName();
			rightPos = secondary.getAlignmentStart();
		} else {
			readIdx = secondary.getReadLength();
			left = secondary.getReferenceName() + "_" + (secondary.getAlignmentStart() + readIdx) + "_" + (secondary.getReadNegativeStrandFlag() ? "R" : "F");
			right = primary.getReferenceName() + "_" + primary.getAlignmentStart() + "_" + (primary.getReadNegativeStrandFlag() ? "R" : "F");
			
			leftChr = secondary.getReferenceName();
			leftPos = secondary.getAlignmentStart() + readIdx;
			rightChr = primary.getReferenceName();
			rightPos = primary.getAlignmentStart();
		}
		
		int startIdx = Math.max(readIdx - readLength, 0);
		int stopIdx = Math.min(readIdx + readLength, primary.getReadLength()-1);
		String svBases = primary.getReadString().substring(startIdx, stopIdx);
		
//		String label = ">" + left + "_" + right + "_" + primary.getReadName();
//		return label + "\n" + svBases;
		
		return new Breakpoint(leftChr, leftPos, rightChr, rightPos, svBases, primary.getReadName());
	}
	
	static class BreakpointGroup {
		private List<Breakpoint> breakpoints = new ArrayList<Breakpoint>();
		private String leftChr;
		private int leftStart;
		private int leftStop;
		private String rightChr;
		private int rightStart;
		private int rightStop;
		private int groupId;
		private Set<String> sequences = new HashSet<String>();
		
		private static int groupCounter = 1;
		
		public BreakpointGroup(Breakpoint breakpoint) {
			this.leftChr = breakpoint.leftChr;
			this.leftStart = breakpoint.leftPos;
			this.leftStop = breakpoint.leftPos;
			this.rightChr = breakpoint.rightChr;
			this.rightStart = breakpoint.rightPos;
			this.rightStop = breakpoint.rightPos;
			
			this.breakpoints.add(breakpoint);
			sequences.add(breakpoint.getBases());
			groupId = groupCounter++;
		}
		
		void addBreakpoint(Breakpoint breakpoint) {
			if (!sequences.contains(breakpoint.getBases())) {
				
				if (breakpoint.leftPos < leftStart) {
					leftStart = breakpoint.leftPos;
				}
				
				if (breakpoint.leftPos > leftStop) {
					leftStop = breakpoint.leftPos;
				}
				
				if (breakpoint.rightPos < rightStart) {
					rightStart = breakpoint.rightPos;
				}
				
				if (breakpoint.rightPos > rightStop) {
					rightStop = breakpoint.rightPos;
				}
				
				breakpoints.add(breakpoint);
				sequences.add(breakpoint.getBases());
			}
		}
		
		boolean shouldMerge(Breakpoint breakpoint) {
			if (leftChr.equals(breakpoint.leftChr) && rightChr.equals(breakpoint.rightChr)) {
				if (breakpoint.leftPos >= leftStart-GROUP_RANGE && breakpoint.leftPos <= leftStart+GROUP_RANGE &&
					breakpoint.rightPos >= rightStart-GROUP_RANGE && breakpoint.rightPos <= rightStart+GROUP_RANGE) {
					
					return true;
				}
			}
			
			return false;
		}
		
		int getGroupId() {
			return groupId;
		}
		
		List<Breakpoint> getBreakpoints() {
			return breakpoints;
		}
	}
	
	static class Breakpoint implements Comparable<Breakpoint> {
		private String leftChr;
		private int leftPos;
		private String rightChr;
		private int rightPos;
		private String bases;
		private String readName;
		
		Breakpoint(String leftChr, int leftPos, String rightChr, int rightPos, String bases, String readName) {
			this.leftChr = leftChr;
			this.leftPos = leftPos;
			this.rightChr = rightChr;
			this.rightPos = rightPos;
			this.bases = bases;
			this.readName = readName;
		}

		@Override
		public int compareTo(Breakpoint that) {
			int compare = this.leftChr.compareTo(that.leftChr);
			
			if (compare == 0) {
				compare = this.leftPos - that.leftPos;
			}
			if (compare == 0) {
				compare = this.rightChr.compareTo(that.rightChr);
			}
			if (compare == 0) {
				compare = this.rightPos - that.rightPos;
			}
			
			return compare;
		}
		
		String getLabel() {
			return leftChr + "_" + leftPos + "_" + rightChr + "_" + rightPos + "_" + readName;
		}
		
		String getBases() {
			return bases;
		}
		
		public String toString() {
			return getLabel();
		}
	}
	
	public static void main(String[] args) throws Exception {
		SVHandler svh = new SVHandler(100);
		svh.identifySVCandidates("/home/lmose/dev/abra/sv/sv_contigs.sam", "/home/lmose/dev/abra/sv/sv_candidates.fa");
	}
}

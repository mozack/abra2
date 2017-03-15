package abra.cadabra;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import abra.SAMRecordUtils;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;

public class SpliceJunctionCounter {
	
	Map<SpliceJunction, Integer> uniqueReads   = new HashMap<SpliceJunction, Integer>();
	Map<SpliceJunction, Integer> multiMapReads = new HashMap<SpliceJunction, Integer>();

	public void countSplices(String input) {
		SamReader reader = SAMRecordUtils.getSamReader(input);
		
		for (SAMRecord read : reader) {
			if (read.getCigarString().contains("N")) {
				for (SpliceJunction junc : getJunctions(read)) {
					incrementCount(junc, read);
				}
			}
		}
		
		List<SpliceJunction> junctions = new ArrayList<SpliceJunction>(uniqueReads.keySet());
		Collections.sort(junctions, new SpliceJunctionComparator(reader.getFileHeader()));
		
		for (SpliceJunction junction : junctions) {
			String rec = String.format("%s\t%d\t%d\t.\t.\t.\t%d\t%d\t.", junction.chrom, junction.start, junction.stop,
					uniqueReads.get(junction), multiMapReads.get(junction));
			
			System.out.println(rec);
		}
	}
	
	private void incrementCount(SpliceJunction junction, SAMRecord read) {
		
		if (!uniqueReads.containsKey(junction)) {
			uniqueReads.put(junction, 0);
		}
		
		if (!multiMapReads.containsKey(junction)) {
			multiMapReads.put(junction, 0);
		}
		
		// TODO: Hardcoded to STAR values here.
		if (read.getMappingQuality() == 255) {
			uniqueReads.put(junction, uniqueReads.get(junction)+1);
		} else {
			multiMapReads.put(junction, multiMapReads.get(junction)+1);
		}
	}
	
	private List<SpliceJunction> getJunctions(SAMRecord read) {
		List<SpliceJunction> junctions = new ArrayList<SpliceJunction>();
		
		int pos = read.getAlignmentStart();
		
		for (CigarElement elem : read.getCigar()) {
			switch (elem.getOperator()) {
				case D:
				case M:
					pos += elem.getLength();
					break;
				case N:
					junctions.add(new SpliceJunction(read.getReferenceName(), pos, pos+elem.getLength()-1));
					pos += elem.getLength();
					break;
				case S:
				case I:
				case H:
					// NOOP
					break;
				default:
					throw new UnsupportedOperationException("Unsupported Cigar Operator: " + elem.getOperator());
			}
		}
		
		return junctions;
	}
	
	static class SpliceJunctionComparator implements Comparator<SpliceJunction> {
		
		private SAMFileHeader header;

		public SpliceJunctionComparator(SAMFileHeader header) {
			this.header = header;
		}
		
		@Override
		public int compare(SpliceJunction j1, SpliceJunction j2) {
			
			int idx1 = header.getSequenceIndex(j1.chrom);
			int idx2 = header.getSequenceIndex(j2.chrom);
			
			if (idx1 != idx2) {
				return idx1-idx2;
			}
			
			if (j1.start != j2.start) {
				return j1.start - j2.start;
			}
			
			return j1.stop - j2.stop;
		}
		
	}
	
	static class SpliceJunction {
		String chrom;
		int start;
		int stop;
		
		public SpliceJunction(String chrom, int start, int stop) {
			this.chrom = chrom;
			this.start = start;
			this.stop  = stop;
		}
		
		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result + ((chrom == null) ? 0 : chrom.hashCode());
			result = prime * result + start;
			result = prime * result + stop;
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
			SpliceJunction other = (SpliceJunction) obj;
			if (chrom == null) {
				if (other.chrom != null)
					return false;
			} else if (!chrom.equals(other.chrom))
				return false;
			if (start != other.start)
				return false;
			if (stop != other.stop)
				return false;
			return true;
		}
	}
	
	public static void main(String[] args) {
		String input = args[0];

		SpliceJunctionCounter counter = new SpliceJunctionCounter();
		counter.countSplices(input);
	}
}

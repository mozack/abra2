package abra;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import abra.SSWAligner.SSWAlignerResult;
import abra.SimpleMapper.SimpleMapperResult;

public class ReadEvaluator {

	// key = SimpleMapper with cached contig, value = contig SW alignment result
	private Map<SimpleMapper, SSWAlignerResult> mappedContigs;
	
	public ReadEvaluator(Map<SimpleMapper, SSWAlignerResult> mappedContigs) {
		this.mappedContigs = mappedContigs;
	}
	
	/**
	 * If an improved alignment exists for the input read, return it.
	 * Returns null if there is no improved alignment
	 * If multiple alignments exist for the read with the same number of mismatches,
	 * the alignments are ambiguous and null is returned.
	 * A read may align to multiple contigs, but result in the same alignment in
	 * the context of the reference.  In this case the alignment is considered distinct. 
	 */
	public Alignment getImprovedAlignment(int origEditDist, String read) {
		Alignment result = null;
		
		List<AlignmentHit> alignmentHits = new ArrayList<AlignmentHit>();
		
		int bestMismatches = origEditDist;
		
		// Map read to all contigs, caching the hits with the smallest number of mismatches
		for (SimpleMapper mapper : mappedContigs.keySet()) {
			SimpleMapperResult mapResult = mapper.map(read);
			
			// TODO: Ambiguous alignment within single contig handling??
			if (mapResult.getMismatches() < bestMismatches) {
				bestMismatches = mapResult.getMismatches();
				alignmentHits.clear();
				alignmentHits.add(new AlignmentHit(mapResult, mapper));
			} else if (mapResult.getMismatches() == bestMismatches && bestMismatches < origEditDist) {
				alignmentHits.add(new AlignmentHit(mapResult, mapper));
			}
		}
		
		// If multiple "best" hits, check to see if they agree.
		Set<Alignment> alignments = new HashSet<Alignment>();
		
		for (AlignmentHit alignmentHit : alignmentHits) {
			SSWAlignerResult contigAlignment = mappedContigs.get(alignmentHit.mapper); 
			
			// Read position in the local reference
			int readRefPos = contigAlignment.getRefPos() + alignmentHit.mapResult.getPos();
			String cigar = CigarUtils.subsetCigarString(alignmentHit.mapResult.getPos(), read.length(), contigAlignment.getCigar());
			
			Alignment readAlignment = new Alignment(readRefPos, cigar, contigAlignment.getRefPos(), contigAlignment.getCigar());
			alignments.add(readAlignment);
		}
		
		// If there is more than 1 distinct alignment, we have an ambiguous result which will not be used.
		if (alignments.size() == 1) {
			result = alignments.iterator().next();
		}
		
		return result;
	}
	
	
	static class AlignmentHit {
		SimpleMapperResult mapResult;
		SimpleMapper mapper;
		
		AlignmentHit(SimpleMapperResult mapResult, SimpleMapper mapper) {
			this.mapResult = mapResult;
			this.mapper = mapper;
		}
	}
	
	//TODO: Genericize this and share
	static class Alignment {
		int pos;
		String cigar;
		
		int contigPos;
		String contigCigar;
		
		Alignment(int pos, String cigar, int contigPos, String contigCigar) {
			this.pos = pos;
			this.cigar = cigar;
			
			this.contigPos = contigPos;
			this.contigCigar = contigCigar;
		}
		
		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result + ((cigar == null) ? 0 : cigar.hashCode());
			result = prime * result + pos;
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
			Alignment other = (Alignment) obj;
			if (cigar == null) {
				if (other.cigar != null)
					return false;
			} else if (!cigar.equals(other.cigar))
				return false;
			if (pos != other.pos)
				return false;
			return true;
		}
	}
}

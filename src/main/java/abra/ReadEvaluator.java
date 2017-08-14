package abra;

import htsjdk.samtools.SAMRecord;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import abra.ContigAligner.ContigAlignerResult;
import abra.SimpleMapper.Orientation;
import abra.SimpleMapper.SimpleMapperResult;

public class ReadEvaluator {

	// key = SimpleMapper with cached contig, value = contig SW alignment result
	private Map<Feature, Map<SimpleMapper, ContigAlignerResult>> mappedContigs;
	
	public ReadEvaluator(Map<Feature, Map<SimpleMapper, ContigAlignerResult>> mappedContigs) {
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
	public Alignment getImprovedAlignment(int origEditDist, SAMRecord samRecord, CompareToReference2 c2r) {
		Alignment result = getImprovedAlignment(origEditDist, samRecord.getReadString());
		
		if (result == null) {
			// If soft clipped, attempt to map only the unclipped portion of the read
			// Relying on the original mapper to handle clipping here
			if (samRecord.getCigarString().contains("S")) {
				
				// Use mapper's edit distance which should be for the unclipped portion of the read only
				int unclippedEditDist = SAMRecordUtils.getEditDistance(samRecord, c2r, false);
				
				// Don't bother remapping if alignment cannot be improved
				if (unclippedEditDist > 0) {
					String unclipped = SAMRecordUtils.getMappedReadPortion(samRecord);
					result = getImprovedAlignment(unclippedEditDist, unclipped);
					if (result != null) {
						result.cigar = SAMRecordUtils.getLeadingClips(samRecord) + result.cigar + SAMRecordUtils.getTrailingClips(samRecord);
					}
				}
			}
		}
		
		return result;
	}
	
	public Alignment getImprovedAlignment(int origEditDist, String read) {
		Alignment result = null;
		
		List<AlignmentHit> alignmentHits = new ArrayList<AlignmentHit>();
		
		int bestMismatches = read.length()+1;
		
		// Map read to all contigs, caching the hits with the smallest number of mismatches
		for (Feature region : mappedContigs.keySet()) {
			Map<SimpleMapper, ContigAlignerResult> regionContigs = mappedContigs.get(region);
			
			for (SimpleMapper mapper : regionContigs.keySet()) {
				SimpleMapperResult mapResult = mapper.map(read);
				
				if (mapResult.getMismatches() < bestMismatches) {
					bestMismatches = mapResult.getMismatches();
					alignmentHits.clear();
					alignmentHits.add(new AlignmentHit(mapResult, mapper, region));
				} else if (mapResult.getMismatches() == bestMismatches) {
					alignmentHits.add(new AlignmentHit(mapResult, mapper, region));
				}
			}
		}
		
		// If multiple "best" hits, check to see if they agree.
		Set<Alignment> alignments = new HashSet<Alignment>();
		
		for (AlignmentHit alignmentHit : alignmentHits) {
			ContigAlignerResult contigAlignment = mappedContigs.get(alignmentHit.region).get(alignmentHit.mapper); 
			
			int readRefPos = alignmentHit.mapResult.getPos();
			String cigar = "";
			
			// Read position in the local reference
			if (alignmentHit.mapResult.getPos() >= 0) {
				StringBuffer cigarBuf = new StringBuffer();
				int readPosInCigarRelativeToRef = CigarUtils.subsetCigarString(alignmentHit.mapResult.getPos(), read.length(), contigAlignment.getCigar(), cigarBuf);
				readRefPos = contigAlignment.getGenomicPos() + readPosInCigarRelativeToRef;
				cigar = cigarBuf.toString();
			}
			
			Alignment readAlignment = new Alignment(contigAlignment.getChromosome(), readRefPos, cigar, alignmentHit.mapResult.getOrientation(), bestMismatches, 
					contigAlignment.getGenomicPos(), contigAlignment.getCigar(), contigAlignment.isSecondary());
			
			if (alignments.size() == 1) {
				Alignment alignment2 = alignments.iterator().next();
				
				int choice = CigarUtils.testEquivalenceAndSelectIntronPreferred(readAlignment, alignment2);
				
				if (choice == 0) {
					// Non-equivalent contigs means we have ambiguous alignments.  Break and return null;
					alignments.clear();
					if (bestMismatches == origEditDist) {
						// Mark as ambiguous if this is same is original alignment.
						return Alignment.AMBIGUOUS;
					}
					break;
				}
				
				// Choose the first alignment
				if (choice == 1) {
					alignments.clear();
					alignments.add(readAlignment);
				}
				
				// Second alignment is already in the Set.
				
			} else if (alignments.size() == 0) {
				alignments.add(readAlignment);
			} else {
				alignments.clear();
				break;
			}
		}
		
		// If there is more than 1 distinct alignment, we have an ambiguous result which will not be used.
		if (alignments.size() == 1) {
			result = alignments.iterator().next();
			// If the result was ambiguously mapped within a single contig, return null.
			if (result.pos < 0) {
				result = null;
			}
		} else if (alignments.size() > 1) {
			// We have multiple possible alignments that are equal to the original alignment
			if (bestMismatches == origEditDist) {
				result = Alignment.AMBIGUOUS;
			}
		}
		
		return result;
	}
	
	
	static class AlignmentHit {
		SimpleMapperResult mapResult;
		SimpleMapper mapper;
		Feature region;
		
		AlignmentHit(SimpleMapperResult mapResult, SimpleMapper mapper, Feature region) {
			this.mapResult = mapResult;
			this.mapper = mapper;
			this.region = region;
		}
	}
	
	//TODO: Genericize this and share ?
	static class Alignment {
		String chromosome;
		int pos;
		String cigar;
		Orientation orientation;
		int numMismatches;
		
		int contigPos;
		String contigCigar;
		boolean isSecondary;
		
		static final Alignment AMBIGUOUS = new Alignment();
		
		Alignment() {
		}
		
		Alignment(String chromosome, int pos, String cigar, Orientation orientation, int numMismatches, int contigPos, String contigCigar, boolean isSecondary) {
			this.chromosome = chromosome;
			this.pos = pos;
			this.cigar = cigar;
			this.orientation = orientation;
			this.numMismatches = numMismatches;
			
			this.contigPos = contigPos;
			this.contigCigar = contigCigar;
			this.isSecondary = isSecondary;
		}

		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result + ((cigar == null) ? 0 : cigar.hashCode());
			result = prime * result
					+ ((orientation == null) ? 0 : orientation.hashCode());
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
			if (orientation != other.orientation)
				return false;
			if (pos != other.pos)
				return false;
			return true;
		}

	}
}

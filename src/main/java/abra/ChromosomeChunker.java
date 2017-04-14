package abra;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import abra.Logger.Level;

public class ChromosomeChunker {
	
	// Minimum chromosome chunk size
	private static final int MIN_CHUNK_SIZE = 50000000;
	
	private CompareToReference2 c2r;
	
	// Chromosome chunks
	private List<Feature> chunks;
	
	// Chromosome chunks grouped by chromosome
	private Map<String, List<Feature>> chunkGroups;
	
	public ChromosomeChunker(CompareToReference2 c2r) {
		this.c2r = c2r;
	}

	// Identify chunks for processing
	// Split chunks at N regions
	public void init() {
		
		chunks = new ArrayList<Feature>();
		chunkGroups = new HashMap<String, List<Feature>>();
		
		for (String chromosome : c2r.getChromosomes()) {
			long currStart = 1;
			int chromosomeLength = c2r.getChromosomeLength(chromosome);
			
			while (currStart < chromosomeLength) {
				Feature chunk = new Feature(chromosome, currStart, Math.min(currStart+MIN_CHUNK_SIZE-1, chromosomeLength));
				currStart = chunk.getEnd()+1;
				chunks.add(chunk);
				
				if (!chunkGroups.containsKey(chromosome)) {
					chunkGroups.put(chromosome, new ArrayList<Feature>());
				}
				
				chunkGroups.get(chromosome).add(chunk);
			}
						
/*			
			List<Feature> nRegions = c2r.getUndefinedRegions().get(chromosome);
			for (Feature nRegion : nRegions) {
				// Only consider N regions of reasonable size
				if (nRegion.getEnd() - nRegion.getStart() > 1000) {
					if (nRegion.getStart() - currStart > MIN_CHUNK_SIZE && nRegion.getStart()+MIN_CHUNK_SIZE < chromosomeLength) {
						Feature chunk = new Feature(chromosome, currStart, nRegion.getStart());
						chunks.add(chunk);
						currStart = nRegion.getStart()+1;
					}
				}
			}
			
			if (currStart < chromosomeLength) {
				chunks.add(new Feature(chromosome, currStart, chromosomeLength));
			}
*/
		}
		
		Logger.info("Chromosome chunks:");
		for (Feature chunk : chunks) {
			Logger.info(chunk.toString());
		}
	}
	
	public List<Feature> getChunks() {
		return chunks;
	}

	public Map<String, List<Feature>> getChunkGroups() {
		return chunkGroups;
	}

	public static void main(String[] args) throws Exception {
		Logger.LEVEL = Level.TRACE;
		//String ref = "/home/lmose/dev/reference/hg38/hg38.fa";
		String ref = "/home/lmose/dev/reference/hg38/chr1.fa";
		
		CompareToReference2 c2r = new CompareToReference2();
		c2r.init(ref);
		
		ChromosomeChunker cc = new ChromosomeChunker(c2r);
		
		cc.init();
		for (Feature chunk : cc.getChunks()) {
			System.out.println(chunk);
		}
	}
}

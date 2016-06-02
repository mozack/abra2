package abra.utils;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import abra.Feature;

/**
 * Orders the input bed file by start and stop position.  Chromosome ordering is preserved.
 * Can be used to correct ordering for use in KmerSizeEvaluator.
 * 
 * @author lmose
 */
public class OrderBed {
	
	public static void orderBed(String bedFile) throws IOException {
		BufferedReader reader = new BufferedReader(new FileReader(bedFile));
		
		try {
			
			String lastChr = "";
			
			String line = reader.readLine();
			
			List<Feature> features = new ArrayList<Feature>();
			
			while (line != null) {
				if(line.startsWith("#")) {
					line = reader.readLine();
					continue;
				}
				String[] fields = line.split("\t");
				
				String chromosome = fields[0];
				long startPos = Long.valueOf(fields[1]);
				long endPos = Long.valueOf(fields[2]);
				
				if (startPos > endPos) {
					throw new IllegalArgumentException("Region end must be greater than or equal to region start in target BED file: " + line);
				}
				
				Feature region = new Feature(chromosome, startPos, endPos);
				
				if (chromosome.equals(lastChr)) {
					features.add(region);
				} else {
					// Sort and output the current chromosome's features
					outputFeatures(features);
				}
				
				line = reader.readLine();
			}
			
			// Output last chromosome
			outputFeatures(features);
		} finally {
			reader.close();
		}
	}

	private static void outputFeatures(List<Feature> features) {
		Collections.sort(features, new FeatureComparator());
		for (Feature feature : features) {
			System.out.println(feature.getSeqname() + "\t" + feature.getStart() + "\t" + feature.getEnd());
		}
		features.clear();
	}
	
	private static class FeatureComparator implements Comparator<Feature> {

		@Override
		public int compare(Feature o1, Feature o2) {
			if (!o1.getSeqname().equals(o2.getSeqname())) {
				throw new IllegalArgumentException(
						"Can't compare regions from different chromosomes: " + o1 + "," + o2);
			}
			
			if (o1.getStart() < o2.getStart()) {
				return -1;
			} else if (o1.getStart() > o2.getStart()) {
				return 1;
			} else if (o1.getEnd() < o2.getEnd()) {
				return -1;
			} else if (o1.getEnd() > o2.getEnd()) {
				return 1;
			} else {
				return 0;
			}
		}

	}
	
	public static void main(String[] args) throws Exception {
		System.err.println("Starting bed ordering.");
		String bedFile = args[0];
		orderBed(bedFile);
		System.err.println("Bed ordering done.");
		
	}
}

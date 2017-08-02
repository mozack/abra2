package abra;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

public class Variant implements Comparable<Variant> {
	
	private String chr;
	private int position;
	private String ref;
	private String alt;
	
	Variant(String chr, int position, String ref, String alt) {
		this.chr = chr;
		this.position = position;
		this.ref = ref;
		this.alt = alt;
	}

	public String getChr() {
		return chr;
	}

	public int getPosition() {
		return position;
	}

	public String getRef() {
		return ref;
	}

	public String getAlt() {
		return alt;
	}
	
	@Override
	public String toString() {
		return "Variant [chr=" + chr + ", position=" + position + ", ref=" + ref + ", alt=" + alt + "]";
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((alt == null) ? 0 : alt.hashCode());
		result = prime * result + ((chr == null) ? 0 : chr.hashCode());
		result = prime * result + position;
		result = prime * result + ((ref == null) ? 0 : ref.hashCode());
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
		Variant other = (Variant) obj;
		if (alt == null) {
			if (other.alt != null)
				return false;
		} else if (!alt.equals(other.alt))
			return false;
		if (chr == null) {
			if (other.chr != null)
				return false;
		} else if (!chr.equals(other.chr))
			return false;
		if (position != other.position)
			return false;
		if (ref == null) {
			if (other.ref != null)
				return false;
		} else if (!ref.equals(other.ref))
			return false;
		return true;
	}

	public int getRefSpan() {
		return ref.length() - alt.length() + 1;
	}
	
	@Override
	public int compareTo(Variant that) {
		return this.position - that.position;
	}
	
	// Assumes all variants on same chromosome.
	static Map<Feature, List<Variant>> groupByRegion(List<Feature> regions, List<Variant> variants) {
		
		// Key = region, value = variants with start in region
		Map<Feature, List<Variant>> regionVariants = new HashMap<Feature, List<Variant>>();
		
		// Key = position, value = list of variants
		Map<Integer, List<Variant>> posVariants = new HashMap<Integer, List<Variant>>();
		
		// Track assigned variants, so that we process only once.
		Set<Variant> assignedVariants = new HashSet<Variant>();
		
		for (Variant variant : variants) {
			if (!posVariants.containsKey(variant.getPosition())) {
				posVariants.put(variant.getPosition(), new ArrayList<Variant>());
			}
			posVariants.get(variant.getPosition()).add(variant);
		}
		
		// For each region
		for (Feature region : regions) {
			// For each position in region
			for (int i= (int) region.getStart(); i < (int) region.getEnd(); i++) {
				List<Variant> currVariants = posVariants.get(i);
				if (currVariants != null) {
					// For each variant at current position
					for (Variant variant : currVariants) {
						
						// If variant not already assigned, assign to curr region.
						if (!assignedVariants.contains(variant)) {
							if (!regionVariants.containsKey(region)) {
								regionVariants.put(region,  new ArrayList<Variant>());
							}
							
							regionVariants.get(region).add(variant);
							assignedVariants.add(variant);
						}
					}
				}
			}
		}
		
		return regionVariants;
	}
	
	/**
	 *  Load variants from VCF return map with key = chromosome, value = variant list sorted by position 
	 */
	static Map<String, List<Variant>> loadFromFile(String vcfFile) throws FileNotFoundException, IOException {
		
		Map<String, List<Variant>> variants = new HashMap<String, List<Variant>>();
		
		BufferedReader reader = new BufferedReader(new FileReader(vcfFile));
		
		String line = reader.readLine();
		while (line != null) {
			
			if (!line.startsWith("#") && !line.trim().isEmpty()) {
				String[] fields = line.split("\\t");
				String chr = fields[0];
				int pos = Integer.parseInt(fields[1]);
				String ref = fields[3];
				String alt = fields[4];
				String filt = fields[6];
				
				if (filt.equalsIgnoreCase("PASS")) {
					Variant variant = new Variant(chr,pos,ref,alt);
					
					if (!variants.containsKey(chr)) {
						variants.put(chr, new ArrayList<Variant>());
					}
					
					variants.get(chr).add(variant);
				}
			}
			
			line = reader.readLine();
		}
		
		reader.close();
		
		for (String chr : variants.keySet()) {
			Collections.sort(variants.get(chr));
		}
		
		return variants;
	}	
}

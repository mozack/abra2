/* Copyright 2013 University of North Carolina at Chapel Hill.  All rights reserved. */
package abra;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * Loads GTF file content into memory.
 * 
 * @author Lisle Mose (lmose at unc dot edu)
 */
public class GtfLoader {
	
	private static final int SEQNAME_IDX = 0;
	private static final int START_IDX   = 3;
	private static final int END_IDX     = 4;

	public List<Feature> load(String gtfFile) throws FileNotFoundException, IOException {
		List<Feature> features = new ArrayList<Feature>();
		
		BufferedReader reader = new BufferedReader(new FileReader(gtfFile));
		
		String line = reader.readLine();
		
		while (line != null) {
			String[] fields = line.split("\t");			
			features.add(new Feature(fields[SEQNAME_IDX], Long.valueOf(fields[START_IDX]), Long.valueOf(fields[END_IDX])));
			
			line = reader.readLine();
		}
		
		reader.close();
		
		return features;
	}
	
	private List<Mut> getLocations(String file) throws Exception {
		
		List<Mut> locations = new ArrayList<Mut>();
		BufferedReader reader = new BufferedReader(new FileReader(file));
		String line = reader.readLine();
		
		while (line != null) {
			String[] fields = line.split("\t");
			if (fields.length >= 2) {
				locations.add(new Mut(fields[0], Integer.parseInt(fields[1])));
			}
			
			line = reader.readLine();
		}
		
		reader.close();
		
		return locations;
	}
	
	public void printSubset(String locationsFile, String gtfFile) throws Exception {
		/*
		List<Mut> mutations = new ArrayList<Mut>();
		Mut mut;
		
		mut = new Mut("chr1",6166686);
		mutations.add(mut);
		mut = new Mut("chr2",16086282);
		mutations.add(mut);
		mut = new Mut("chr3",12627213);
		mutations.add(mut);
		mut = new Mut("chr4",55095482);
		mutations.add(mut);
		mut = new Mut("chr5",1255521);
		mutations.add(mut);
		mut = new Mut("chr6",32164773);
		mutations.add(mut);
		mut = new Mut("chr7",33003202);
		mutations.add(mut);
		mut = new Mut("chr8",27170777);
		mutations.add(mut);
		mut = new Mut("chr9",4985979);
		mutations.add(mut);
		mut = new Mut("chr10",30728218);
		mutations.add(mut);
		mut = new Mut("chr11",4116181);
		mutations.add(mut);
		mut = new Mut("chr12",12868034);
		mutations.add(mut);
		mut = new Mut("chr13",26956993);
		mutations.add(mut);
		mut = new Mut("chr14",36988759);
		mutations.add(mut);
		mut = new Mut("chr15",66735694);
		mutations.add(mut);
		mut = new Mut("chr16",2589116);
		mutations.add(mut);
		mut = new Mut("chr17",7573991);
		mutations.add(mut);
		mut = new Mut("chr18",45364526);
		mutations.add(mut);
		mut = new Mut("chr19",932455);
		mutations.add(mut);
		mut = new Mut("chr20",9543644);
		mutations.add(mut);
		mut = new Mut("chr21",36162224);
		mutations.add(mut);
		mut = new Mut("chr22",21288485);
		mutations.add(mut);
		*/
		
		List<Mut> mutations = getLocations(locationsFile);
		
		BufferedReader reader = new BufferedReader(new FileReader(gtfFile));
		
		String line = reader.readLine();
		
		while (line != null) {
			String[] fields = line.split("\t");			
			Feature f = new Feature(fields[SEQNAME_IDX], Long.valueOf(fields[START_IDX]), Long.valueOf(fields[END_IDX]));
			
			for (Mut mut1 : mutations) {
				if ((f.getSeqname().equals(mut1.chr)) && (f.getStart() <= mut1.pos) && (f.getEnd() >= mut1.pos)) {
					System.out.println(line);
				}
			}
			
			line = reader.readLine();
		}
		
		reader.close();

	}
	
	static class Mut {
		String chr;
		int pos;
		
		Mut(String chr, int pos) {
			this.chr = chr;
			this.pos = pos;
		}
	}
	
	public static void main(String[] args) throws Exception {
		GtfLoader gl = new GtfLoader();
		gl.printSubset("/home/lmose/dev/ayc/sim/sim528/regions.txt", "/home/lmose/dev/ayc/regions/clinseq5/uncseq5.gtf");
	}
}

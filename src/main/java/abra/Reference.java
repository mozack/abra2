/* Copyright 2013 University of North Carolina at Chapel Hill.  All rights reserved. */
package abra;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

/**
 * Representation of reference
 * @author Lisle E. Mose (lmose at unc dot edu)
 */
public class Reference {

	//TODO: Delete this class in favor CompareToReference2
	private Map<String, StringBuffer> refMap = new HashMap<String, StringBuffer>();
	
	public Reference(String reference) throws IOException {
		
		BufferedReader refReader = new BufferedReader(new FileReader(reference));
		String line = refReader.readLine();
		
		String key = null;
		StringBuffer seq = null;
		
		while (line != null) {
			
			if (line.startsWith(">")) {
				if (key != null) {
					refMap.put(key, seq);
				}
				
				int endIdx = line.indexOf(' ');
				if (endIdx < 0) {
					endIdx = line.length();
				}
				key = line.substring(1, endIdx);
				seq = new StringBuffer();
			} else {
				seq.append(line);
			}
			
			line = refReader.readLine();
		}
		
		if (key != null) {
			refMap.put(key, seq);
		}
		
		refReader.close();
	}
	
	public String getSequence(String chromosome, int position, int length) {
		StringBuffer ref = refMap.get(chromosome);
		
		if (ref == null) {
			System.out.println("No ref for chromosome: " + chromosome);
		}
		
		position -= 1;
		
		return ref.substring(Math.max(position, 0), Math.min(position + length, ref.length()));
	}
}

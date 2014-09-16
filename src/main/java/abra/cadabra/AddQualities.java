package abra.cadabra;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

/**
 * Temporary code for adding quality scores to existing VCF files.
 * 
 * @author lmose
 */
public class AddQualities {

	public void run(String filename, boolean isMerged) throws FileNotFoundException, IOException {
		BufferedReader reader = new BufferedReader(new FileReader(filename));
		String line = reader.readLine();
		
		while (line != null) {
			if (!line.startsWith("#") && (line.contains("CADABRA") || !isMerged)) {
				line = addQual(line);
				System.out.println(line);
			}
			
			line = reader.readLine();
		}
		
		reader.close();
	}
	
	private String addQual(String line) {
		String[] fields = line.split("\t");
		String[] normalAd = fields[9].split(":")[1].split(",");
		String[] tumorAd = fields[10].split(":")[1].split(",");
		
		int normalRef = Integer.parseInt(normalAd[0]);
		int normalAlt = Integer.parseInt(normalAd[1]);
		int tumorRef = Integer.parseInt(tumorAd[0]);
		int tumorAlt = Integer.parseInt(tumorAd[1]);
		
		double qual = Cadabra.calcPhredScaledQuality(normalRef, normalAlt, tumorRef, tumorAlt);
		fields[5] = String.valueOf(qual);
		
		StringBuffer updatedLine = new StringBuffer();
		for (String field : fields) {
			updatedLine.append(field);
		}
		
		return updatedLine.toString();
	}
	
	public static void main(String[] args) throws Exception {
		String filename = args[0];
		boolean isMerged = args.length > 1 && args[1].equalsIgnoreCase("merged");
		
		AddQualities aq = new AddQualities();
		aq.run(filename,  isMerged);
		
		System.err.println("Add qualities done.");
	}
}

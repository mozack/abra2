/* Copyright 2013 University of North Carolina at Chapel Hill.  All rights reserved. */
package abra;

import java.io.File;
import java.io.IOException;

import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMRecord;

public class FixEditDistance {

	public void fix(String input, String output, String reference) throws IOException {
		SAMFileReader reader = new SAMFileReader(new File(input));
		reader.setValidationStringency(ValidationStringency.SILENT);
		
		SAMFileWriter writer = new SAMFileWriterFactory().makeSAMOrBAMWriter(
				reader.getFileHeader(), true, new File(output));
		
		CompareToReference c2r = new CompareToReference();
		c2r.init(reference);
		
		int i=0;
		
		for (SAMRecord read : reader) {
			int numMismatches = c2r.numMismatches(read);
			int numIndelBases = getNumIndelBases(read);
			int editDistance = numMismatches + numIndelBases;
			
			read.setAttribute("NM", editDistance);
			
			writer.addAlignment(read);
			
			if ((i % 1000000) == 0) {
				System.out.println("Processed: " + i + " reads.");
			}
			i++;
		}
		
		writer.close();
		reader.close();
		
		System.out.println("Done.");
	}
	
	private int getNumIndelBases(SAMRecord read) {
		int numIndelBases = 0;
		
		for (CigarElement element : read.getCigar().getCigarElements()) {
			if ((element.getOperator() == CigarOperator.D) || (element.getOperator() == CigarOperator.I)) {
				numIndelBases += element.getLength();
			}
		}
		
		return numIndelBases;
	}
	
	public static void main(String[] args) throws Exception {
		new FixEditDistance().fix(args[0], args[1], args[2]);
	}
}

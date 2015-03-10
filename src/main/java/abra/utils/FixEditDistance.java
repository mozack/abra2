/* Copyright 2013 University of North Carolina at Chapel Hill.  All rights reserved. */
package abra.utils;

import java.io.File;
import java.io.IOException;

import abra.CompareToReference;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.SAMRecord;

/**
 * Sets NM tag to num mismatches plus total indel length.
 * 
 * @author Lisle E. Mose (lmose at unc dot edu)
 */
public class FixEditDistance {

	//TODO: Move to utils package
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

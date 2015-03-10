/* Copyright 2013 University of North Carolina at Chapel Hill.  All rights reserved. */
package htsjdk.samtools;

import java.io.ByteArrayInputStream;
import java.io.InputStream;

import htsjdk.samtools.ValidationStringency;

public class MyReader {

	public static SAMRecord getRead(String str, SAMFileHeader header) {
		
		
		SAMFileReader.setDefaultValidationStringency(ValidationStringency.SILENT);
		SAMTextReader2 reader = new SAMTextReader2(str, ValidationStringency.SILENT, new DefaultSAMRecordFactory(), header);
//		SAMTextReader2 reader = new SAMTextReader2(is, ValidationStringency.SILENT, new DefaultSAMRecordFactory(), header);
		
		SAMRecord read = reader.getIterator().next();
		
		return read;
	}
}

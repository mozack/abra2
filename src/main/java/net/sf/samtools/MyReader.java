package net.sf.samtools;

import java.io.ByteArrayInputStream;
import java.io.InputStream;

import net.sf.samtools.SAMFileReader.ValidationStringency;

public class MyReader {

	public static SAMRecord getRead(String str, SAMFileHeader header) {
		
		
		SAMFileReader.setDefaultValidationStringency(ValidationStringency.SILENT);
		SAMTextReader2 reader = new SAMTextReader2(str, ValidationStringency.SILENT, new DefaultSAMRecordFactory(), header);
//		SAMTextReader2 reader = new SAMTextReader2(is, ValidationStringency.SILENT, new DefaultSAMRecordFactory(), header);
		
		SAMRecord read = reader.getIterator().next();
		
		return read;
	}
}

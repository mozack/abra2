package abra.utils;

import java.io.File;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMRecord;

public class CalcReadMovement {

	public static void getReadDistance(String file) {
		SAMFileReader reader = new SAMFileReader(new File(file));
		reader.setValidationStringency(ValidationStringency.SILENT);
		
		int numReads = 0;
		for (SAMRecord read : reader) {
			numReads++;
			String yo = (String) read.getAttribute("YO");
			//chr17:7579659:-:100M
			
			if (yo != null) {
			
				if (yo.equals("N/A")) {
					System.out.println("N/A\t" + read.getReferenceName());
				} else {
					String[] fields = yo.split(":");
					String chr = fields[0];
					int pos = Integer.parseInt(fields[1]);
					String strand = fields[2];
				
					StringBuffer s = new StringBuffer();
					if (!read.getReferenceName().equals(chr)) {
						s.append('C');
					}
					
					String readStrand = read.getReadNegativeStrandFlag() ? "-" : "+";
					
					if (!strand.equals(readStrand)) {
						s.append("S");
					}
					
					int diff = Math.abs(read.getAlignmentStart() - pos);
					s.append("\t");
					s.append(diff);
					s.append("\t");
					
					if (!read.getReferenceName().equals(chr)) {
						s.append(chr + "->" + read.getReferenceName());
					}
					
					System.out.println(s.toString());
				}
			}
		}
	}
	
	public static void main(String[] args) {
		String pid = args[0];
		String baseDir = args[1];
		String abraBam = args[2];

		String bam = baseDir + "/" + pid + "/" + abraBam;
		getReadDistance(bam);
	}
}

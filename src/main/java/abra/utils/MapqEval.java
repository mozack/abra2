package abra.utils;

import java.io.File;
import java.io.IOException;
import java.util.Iterator;
import java.util.List;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMFileReader.ValidationStringency;

import abra.Feature;
import abra.RegionLoader;

public class MapqEval {

	public void eval(String input, String bed) throws IOException {
		RegionLoader loader = new RegionLoader();
		List<Feature> regions = loader.load(bed);
		SAMFileReader reader = new SAMFileReader(new File(input));
		reader.setValidationStringency(ValidationStringency.SILENT);

		int[] thresholds = new int[] { 0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 255 };
		
		
		for (Feature region : regions) {
			int[] counts = new int[thresholds.length];
			Iterator<SAMRecord> iter = reader.queryOverlapping(region.getSeqname(), (int) region.getStart(), (int) region.getEnd());
			while (iter.hasNext()) {
				SAMRecord read = iter.next();
				for (int i=0; i<thresholds.length; i++) {
					if (read.getMappingQuality() <= thresholds[i]) {
						counts[i] += 1;
					}
				}
			}
			
			StringBuffer output = new StringBuffer();
			output.append(region.getDescriptor());
			
			for (int count : counts) {
				output.append('\t');
				output.append(count);
			}
			
			System.out.println(output.toString());
		}
		
		reader.close();
	}
	
	public static void main(String[] args) throws Exception {
		new MapqEval().eval(args[0], args[1]);
	}
}

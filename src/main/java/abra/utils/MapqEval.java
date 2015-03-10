package abra.utils;

import java.io.File;
import java.io.IOException;
import java.util.Iterator;
import java.util.List;

import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.ValidationStringency;

import abra.Feature;
import abra.RegionLoader;

public class MapqEval {

	public void eval(String input, String bed) throws IOException {
		RegionLoader loader = new RegionLoader();
		List<Feature> regions = loader.load(bed, false);

		int[] thresholds = new int[] { 0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 255 };
		
		
		for (Feature region : regions) {
			SAMFileReader reader = new SAMFileReader(new File(input));
			reader.setValidationStringency(ValidationStringency.SILENT);

			int[] counts = new int[thresholds.length];
			System.err.println(region.getDescriptor());
			Iterator<SAMRecord> iter = reader.queryOverlapping(region.getSeqname(), (int) region.getStart(), (int) region.getEnd());
			while (iter.hasNext()) {
				SAMRecord read = iter.next();
				if (!read.getReadUnmappedFlag()) {
					for (int i=0; i<thresholds.length; i++) {
						if (read.getMappingQuality() <= thresholds[i]) {
							counts[i] += 1;
						}
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
			reader.close();
		}
	}
	
	public static void main(String[] args) throws Exception {
		new MapqEval().eval(args[0], args[1]);
	}
}

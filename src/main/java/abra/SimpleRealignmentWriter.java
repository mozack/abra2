package abra;

import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMRecord;

public class SimpleRealignmentWriter implements RealignmentWriter {

	private SAMFileWriter writer;
	private ReAligner realigner;
	private int realignCount = 0;
	
	public SimpleRealignmentWriter(ReAligner realigner, SAMFileWriter writer) {
		this.writer = writer;
		this.realigner = realigner;
	}
	
	@Override
	public void addAlignment(SAMRecord contigAlignedRead,
			SAMRecord updatedRead, SAMRecord origRead) {
		
		if (updatedRead != null) {
			// Output realigned read
			writer.addAlignment(updatedRead);
			if (updatedRead.getAttribute("YO") != null) {
				realignCount += 1;
			}
		} else {
			// Output original read
			realigner.adjustForStrand(contigAlignedRead.getReadNegativeStrandFlag(), origRead);
			writer.addAlignment(origRead);
		}
	}

	@Override
	public int flush() {
		return realignCount;
	}

}

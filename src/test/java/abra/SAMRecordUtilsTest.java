package abra;

import net.sf.samtools.SAMRecord;

import org.testng.Assert;
import org.testng.annotations.Test;

public class SAMRecordUtilsTest {

	@Test (groups = "unit" )
	public void testRemoveSoftClips_withDeletionAndSoftClipAtEnd() {
		SAMRecord read = new SAMRecord(null);
		
		read.setReadName("TEST1");
		read.setCigarString("5M2D3M2S");
		read.setReadString("CCCCCCAGCC");
		
		SAMRecordUtils.removeSoftClips(read);
		
		Assert.assertEquals(read.getCigarString(), "5M2D3M");
		Assert.assertEquals(read.getReadString(), "CCCCCCAG");
	}
	
	@Test (groups = "unit" )
	public void testRemoveSoftClips_withDeletionAndSoftClipAtStart() {
		SAMRecord read = new SAMRecord(null);
		
		read.setReadName("TEST1");
		read.setCigarString("2S5M2D3M");
		read.setReadString("CCCCCCAGCC");
		
		SAMRecordUtils.removeSoftClips(read);
		
		Assert.assertEquals(read.getCigarString(), "5M2D3M");
		Assert.assertEquals(read.getReadString(), "CCCCAGCC");
	}
}

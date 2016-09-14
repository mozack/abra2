package abra;

import java.util.HashMap;
import java.util.Map;

import static org.testng.Assert.assertEquals;
import org.testng.annotations.Test;

import abra.ReadEvaluator.Alignment;
import abra.SSWAligner.SSWAlignerResult;

public class ReadEvaluatorTest {

	@Test (groups="unit")
	public void testSingleAlignmentSingleContig() {
		String contig1 = "ATCGAAAAAATTTTTTCCCCCCGGGGGGATCGGCTAATCG";
		String read    =     "ATAAAATTTTTTCCCCCCGGGGGGATCG";  // matches contig at 0 based position 4 with 1 mismatch
		
		SimpleMapper sm1 = new SimpleMapper(contig1);
		SSWAlignerResult swc1 = new SSWAlignerResult(10, "10M1D30M");
		
		Map<SimpleMapper, SSWAlignerResult> mappedContigs = new HashMap<SimpleMapper, SSWAlignerResult>();
		mappedContigs.put(sm1, swc1);
		
		ReadEvaluator re = new ReadEvaluator(mappedContigs);
		
		// 1 mismatch in alignment to contig versus edit distance 2 in original read
		// should result in an improved alignment
		Alignment alignment = re.getImprovedAlignment(2, read);
		assertEquals(alignment.pos, 14);  // Alignment pos = 10 + 4
		assertEquals(alignment.cigar, "6M1D22M");
		assertEquals(alignment.numMismatches, 1);
	}
	
	@Test (groups="unit")
	public void testSingleAlignmentSingleContig_noImprovement() {
		String contig1 = "ATCGAAAAAATTTTTTCCCCCCGGGGGGATCGGCTAATCG";
		String read    =     "ATAAAATTTTTTCCCCCCGGGGGGATCG";  // matches contig at 0 based position 4 with 1 mismatch
		
		SimpleMapper sm1 = new SimpleMapper(contig1);
		SSWAlignerResult swc1 = new SSWAlignerResult(10, "10M1D30M");
		
		Map<SimpleMapper, SSWAlignerResult> mappedContigs = new HashMap<SimpleMapper, SSWAlignerResult>();
		mappedContigs.put(sm1, swc1);
		
		ReadEvaluator re = new ReadEvaluator(mappedContigs);
		
		// 1 mismatch in alignment to contig versus edit distance 1 in original read
		// should result in no improved alignment
		Alignment alignment = re.getImprovedAlignment(1, read);
		assertEquals(alignment, null);
	}
	
	@Test (groups="unit")
	public void testSelectBestAlignment() {
		String contig1 = "ATCGAAAAAATTTTTTCCCCCCGGGGGGATCGGCTAATCG";
		String contig2 = "ATCGAAAAAATTTTTTCCCCCCGGGGGGATCGGCTTATCG";
		String contig3 = "ATCGAAAAAATTTTTTCCCCCCGGGGGGATCGGCTCATCG";
		String contig4 = "ATCGATAAAATTTTTTCCCCCCGGGGGGATCGGCTAATCG";
		String read    =     "ATAAAATTTTTTCCCCCCGGGGGGATCG";  // matches contig at 0 based position 4 with 0 mismatches
		
		SimpleMapper sm1 = new SimpleMapper(contig1);
		SSWAlignerResult swc1 = new SSWAlignerResult(10, "10M1D30M");

		SimpleMapper sm2 = new SimpleMapper(contig2);
		SSWAlignerResult swc2 = new SSWAlignerResult(20, "10M1D30M");

		SimpleMapper sm3 = new SimpleMapper(contig3);
		SSWAlignerResult swc3 = new SSWAlignerResult(30, "10M1D30M");

		SimpleMapper sm4 = new SimpleMapper(contig4);
		SSWAlignerResult swc4 = new SSWAlignerResult(40, "10M1D30M");

		Map<SimpleMapper, SSWAlignerResult> mappedContigs = new HashMap<SimpleMapper, SSWAlignerResult>();
		mappedContigs.put(sm1, swc1);
		mappedContigs.put(sm2, swc2);
		mappedContigs.put(sm3, swc3);
		mappedContigs.put(sm4, swc4);
		
		ReadEvaluator re = new ReadEvaluator(mappedContigs);
		
		// Exact match to contig 4
		Alignment alignment = re.getImprovedAlignment(2, read);
		assertEquals(alignment.pos, 44);  // Alignment pos = 40 + 4
		assertEquals(alignment.cigar, "6M1D22M");
	}
	
	@Test (groups="unit")
	public void testMapToMultipleContigsSynonymously() {
		String contig1 = "ATCGAAAAAATTTTTTCCCCCCGGGGGGATCGGCTAATCG";
		String contig2 = "AATCGAAAAAATTTTTTCCCCCCGGGGGGATCGGCTTATCG";
		String contig3 = "TCGAAAAAATTTTTTCCCCCCGGGGGGATCGGCTCATCG";
		String contig4 = "ATCGAAAAAATTTTTTCCCCCCGGGGGGATCGGCTAATCG";
		String read    =     "ATAAAATTTTTTCCCCCCGGGGGGATCG";  // matches contig at 0 based position 4 with 0 mismatches
		
		SimpleMapper sm1 = new SimpleMapper(contig1);
		SSWAlignerResult swc1 = new SSWAlignerResult(10, "10M1D30M");

		SimpleMapper sm2 = new SimpleMapper(contig2);
		SSWAlignerResult swc2 = new SSWAlignerResult(9, "11M1D31M");

		SimpleMapper sm3 = new SimpleMapper(contig3);
		SSWAlignerResult swc3 = new SSWAlignerResult(11, "9M1D29M");

		SimpleMapper sm4 = new SimpleMapper(contig4);
		SSWAlignerResult swc4 = new SSWAlignerResult(10, "10M1D30M");

		Map<SimpleMapper, SSWAlignerResult> mappedContigs = new HashMap<SimpleMapper, SSWAlignerResult>();
		mappedContigs.put(sm1, swc1);
		mappedContigs.put(sm2, swc2);
		mappedContigs.put(sm3, swc3);
		mappedContigs.put(sm4, swc4);
		
		ReadEvaluator re = new ReadEvaluator(mappedContigs);
		
		// Maps to multiple contigs with a single mismatch, with each contig's
		// alignment result identical in the context of the reference
		Alignment alignment = re.getImprovedAlignment(2, read);
		assertEquals(alignment.pos, 14);  // Alignment pos = 40 + 4
		assertEquals(alignment.cigar, "6M1D22M");
	}
	
	@Test (groups="unit")
	public void testMultimapWithinContig() {
		String contig1 = "ATCGATCGATCGATCGATCGATCGATCGATCGATCG";
		String read    = "ACCGATCGATCGATCGATCGATCGATCGATCG";  // matches 2 locations with single mismatch
		
		SimpleMapper sm1 = new SimpleMapper(contig1);
		SSWAlignerResult swc1 = new SSWAlignerResult(100, "36M");
		
		Map<SimpleMapper, SSWAlignerResult> mappedContigs = new HashMap<SimpleMapper, SSWAlignerResult>();
		mappedContigs.put(sm1, swc1);
		
		ReadEvaluator re = new ReadEvaluator(mappedContigs);
		
		// 1 mismatch in alignment to contig versus edit distance 2 in original read
		// should result in an improved alignment
		Alignment alignment = re.getImprovedAlignment(2, read);
		assertEquals(alignment, null);
	}
}

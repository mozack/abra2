package abra;

import static abra.Logger.log;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMFileReader.ValidationStringency;

/**
 * Responsible for adjusting read alignments.
 * This class is used by multiple threads simultaneously.  Do not store thread
 * specific variables as members.
 * 
 * @author Lisle E. Mose (lmose at unc dot edu)
 */
public class ReadAdjuster {
	
	private int maxMapq;
	private ReverseComplementor reverseComplementor = new ReverseComplementor();
	private boolean isPairedEnd;
	private CompareToReference2 c2r;
	private int minInsertLen;
	private int maxInsertLen;
	
	private static final String ORIGINAL_ALIGNMENT_TAG = "YO";
	private static final String MISMATCHES_TO_CONTIG_TAG = "YM";
	private static final String CONTIG_QUALITY_TAG = "YQ";
	private static final String CONTIG_ALIGNMENT_TAG = "YA";
	
	public ReadAdjuster(boolean isPairedEnd, int maxMapq, CompareToReference2 c2r, int minInsertLen, int maxInsertLen) {
		this.isPairedEnd = isPairedEnd;
		this.maxMapq = maxMapq;
		this.minInsertLen = minInsertLen;
		this.maxInsertLen = maxInsertLen;
		this.c2r = c2r;
	}
	
	public void adjustReads(String alignedToContigSam, SAMFileWriter outputSam, boolean isTightAlignment,
			String tempDir, SAMFileHeader samHeader) throws IOException {
		
		log("Adjusting reads.");
		
		RealignmentWriter writer = getRealignmentWriter(outputSam, isTightAlignment, tempDir);
		
		SAMFileReader contigReader = new SAMFileReader(new File(alignedToContigSam));
		contigReader.setValidationStringency(ValidationStringency.SILENT);
		
		SamStringReader samStringReader = new SamStringReader(samHeader);
		
		for (SAMRecord read : contigReader) {
			
			String origSamStr = read.getReadName();
			origSamStr = origSamStr.replace(Sam2Fastq.FIELD_DELIMITER, "\t");
			SAMRecord orig;
			try {
				orig = samStringReader.getRead(origSamStr);
			} catch (RuntimeException e) {
				System.out.println("Error processing: [" + origSamStr + "]");
				System.out.println("Contig read: [" + read.getSAMString() + "]");
				e.printStackTrace();
				throw e;
			}
			orig.setHeader(samHeader);
			
			orig.setReadString(read.getReadString());
			orig.setBaseQualityString(read.getBaseQualityString());

			SAMRecord readToOutput = null;
			
			// Only adjust reads that align to contig with no indel and shorter edit distance than the original alignment
			String matchingString = read.getReadLength() + "M";
			if ((read.getCigarString().equals(matchingString)) &&
				(read.getReadUnmappedFlag() == false)  &&
				(!orig.getCigarString().contains("N")) &&  // Don't remap introns
				(SAMRecordUtils.getEditDistance(read, null) < SAMRecordUtils.getOrigEditDistance(orig)) &&
				(!isFiltered(orig))) {
				
				SAMRecord origRead = orig;
				String contigReadStr = read.getReferenceName();
				
				int numBestHits = SAMRecordUtils.getIntAttribute(read, "X0");
				int numSubOptimalHits = SAMRecordUtils.getIntAttribute(read, "X1");
				
				int totalHits = numBestHits + numSubOptimalHits;
				
				List<HitInfo> bestHits = getBestHits(contigReadStr, samStringReader, read, matchingString);
				
				Map<String, SAMRecord> outputReadAlignmentInfo = convertBestHitsToAlignmentInfo(bestHits, origRead);
				
				readToOutput = getUpdatedReadInfo(outputReadAlignmentInfo, read, 
						orig, origRead, c2r, totalHits, isTightAlignment);
				
			}
			
			adjustForStrand(read.getReadNegativeStrandFlag(), orig);
			writer.addAlignment(readToOutput, orig);
		}

		int realignedCount = writer.flush();
		contigReader.close();
		
		log("Done adjusting reads.  Number of reads realigned: " + realignedCount);
	}
	
	private SAMRecord getUpdatedReadInfo(Map<String, SAMRecord> outputReadAlignmentInfo, SAMRecord read, 
			SAMRecord orig, SAMRecord origRead, CompareToReference2 c2r, int totalHits, boolean isTightAlignment) {
		SAMRecord readToOutput = null;
		
		if (outputReadAlignmentInfo.size() == 1) {
			readToOutput = outputReadAlignmentInfo.values().iterator().next();
			
			// Check to see if the original read location was in a non-target region.
			// If so, compare updated NM to ref versus original NM to ref
			if ((c2r != null) && (!isImprovedAlignment(readToOutput, orig, c2r))) {
				readToOutput = null;
			} else {
				int origBestHits = SAMRecordUtils.getIntAttribute(readToOutput, "X0");
				int origSuboptimalHits = SAMRecordUtils.getIntAttribute(readToOutput, "X1");
				
				// If the read mapped to multiple locations, set mapping quality to zero.
				if ((outputReadAlignmentInfo.size() > 1) || (totalHits > 1000)) {
					readToOutput.setMappingQuality(0);
				}
				
				// This must happen prior to updateMismatchAndEditDistance
				adjustForStrand(read.getReadNegativeStrandFlag(), readToOutput);
				
				if (readToOutput.getAttribute(ORIGINAL_ALIGNMENT_TAG) != null) {
					// HACK: Only add X0 for final alignment.  Assembler skips X0 > 1
					if (isTightAlignment) {
						readToOutput.setAttribute("X0", outputReadAlignmentInfo.size());
					} else {
						readToOutput.setAttribute("X0", null);
					}
					readToOutput.setAttribute("X1", origBestHits + origSuboptimalHits);
					
					// Clear various tags
					readToOutput.setAttribute("XO", null);
					readToOutput.setAttribute("XG", null);
					readToOutput.setAttribute("MD", null);
					readToOutput.setAttribute("XA", null);
					readToOutput.setAttribute("XT", null);
					
					if (c2r != null) {
						updateMismatchAndEditDistance(readToOutput, c2r, origRead);
					}
				}
			}
		}
		
		return readToOutput;
	}
	
	private Map<String, SAMRecord> convertBestHitsToAlignmentInfo(List<HitInfo> bestHits, SAMRecord origRead) {
		
		Map<String, SAMRecord> outputReadAlignmentInfo = new HashMap<String, SAMRecord>();
		
		for (HitInfo hitInfo : bestHits) {
			
			SAMRecord contigRead = hitInfo.getRecord();
			int position = hitInfo.getPosition() - 1;

			List<ReadBlock> contigReadBlocks = ReadBlock.getReadBlocks(contigRead);
			
			ReadPosition readPosition = new ReadPosition(origRead, position, -1);
			SAMRecord updatedRead = updateReadAlignment(contigRead,
					contigReadBlocks, readPosition);
			
			if (updatedRead != null) {						
				if (updatedRead.getReadUnmappedFlag()) {
					updatedRead.setReadUnmappedFlag(false);
				}
				
				updatedRead.setReadNegativeStrandFlag(hitInfo.isOnNegativeStrand());
				
				// If the read's alignment info has been modified, record the original alignment.
				if (origRead.getReadUnmappedFlag() ||
					!origRead.getReferenceName().equals(updatedRead.getReferenceName()) ||
					origRead.getAlignmentStart() != updatedRead.getAlignmentStart() ||
					origRead.getReadNegativeStrandFlag() != updatedRead.getReadNegativeStrandFlag() ||
					!origRead.getCigarString().equals(updatedRead.getCigarString())) {
				
					if (SAMRecordUtils.isSoftClipEquivalent(origRead, updatedRead)) {
						// Restore Cigar and position
//							System.out.println("Re-setting [" + updatedRead.getSAMString() + "] --- [" + origRead.getSAMString() + "]");
						updatedRead.setAlignmentStart(origRead.getAlignmentStart());
						updatedRead.setCigar(origRead.getCigar());
						
					} else {
						String originalAlignment;
						if (origRead.getReadUnmappedFlag()) {
							originalAlignment = "N/A";
						} else {
							originalAlignment = origRead.getReferenceName() + ":" + origRead.getAlignmentStart() + ":" +
									(origRead.getReadNegativeStrandFlag() ? "-" : "+") + ":" + origRead.getCigarString();
						}
						
						// Read's original alignment position
						updatedRead.setAttribute(ORIGINAL_ALIGNMENT_TAG, originalAlignment);
					}
				}
				
				// Mismatches to the contig
				updatedRead.setAttribute(MISMATCHES_TO_CONTIG_TAG, hitInfo.getNumMismatches());
				
				// Contig's mapping quality
				updatedRead.setAttribute(CONTIG_QUALITY_TAG, hitInfo.getRecord().getMappingQuality());
				
				// Contig's length
//				updatedRead.setAttribute("YL", hitInfo.getRecord().getCigar().getReadLength());
				
				// Contig Position + CIGAR
				updatedRead.setAttribute(CONTIG_ALIGNMENT_TAG, contigRead.getReferenceName() + ":" + contigRead.getAlignmentStart() +
						":" + contigRead.getCigarString());
				
				//TODO: Check strand!!!
				String readAlignmentInfo = updatedRead.getReferenceName() + "_" + updatedRead.getAlignmentStart() + "_" +
						updatedRead.getCigarString();
				
				if (!outputReadAlignmentInfo.containsKey(readAlignmentInfo)) {
					outputReadAlignmentInfo.put(readAlignmentInfo, updatedRead);
				}
			}
		}
		
		return outputReadAlignmentInfo;
	}
	
	private List<HitInfo> getBestHits(String contigReadStr, SamStringReader samStringReader,
			SAMRecord read, String matchingString) {
		List<HitInfo> bestHits = new ArrayList<HitInfo>();

		contigReadStr = contigReadStr.substring(contigReadStr.indexOf('~')+1);
		contigReadStr = contigReadStr.replace('~', '\t');
		SAMRecord contigRead = samStringReader.getRead(contigReadStr);
		
		int bestMismatches = SAMRecordUtils.getIntAttribute(read, "XM");
		
		// Filter this hit if it aligns past the end of the contig
		// Must use cigar length instead of read length, because the
		// the contig read bases are not loaded.
		if (read.getAlignmentEnd() <= contigRead.getCigar().getReadLength()) {
			HitInfo hit = new HitInfo(contigRead, read.getAlignmentStart(),
					read.getReadNegativeStrandFlag() ? '-' : '+', bestMismatches);
			
			bestHits.add(hit);
		}
		
		int numBestHits = SAMRecordUtils.getIntAttribute(read, "X0");
		int numSubOptimalHits = SAMRecordUtils.getIntAttribute(read, "X1");
		
		int totalHits = numBestHits + numSubOptimalHits;
		
		if ((totalHits > 1) && (totalHits < 1000)) {
			
			// Look in XA tag.
			String alternateHitsStr = (String) read.getAttribute("XA");
			if (alternateHitsStr == null) {
				String msg = "best hits = " + numBestHits + ", but no XA entry for: " + read.getSAMString();
				System.out.println(msg);							
			} else {
				
				String[] alternates = alternateHitsStr.split(";");
				
				for (int i=0; i<alternates.length-1; i++) {
					String[] altInfo = alternates[i].split(",");
					String altContigReadStr = altInfo[0];
					char strand = altInfo[1].charAt(0);
					int position = Integer.parseInt(altInfo[1].substring(1));
					String cigar = altInfo[2];
					int mismatches = Integer.parseInt(altInfo[3]);
					
					if ((cigar.equals(matchingString)) && (mismatches < bestMismatches)) {
						System.out.println("MISMATCH_ISSUE: " + read.getSAMString());
					}
					
					if ((cigar.equals(matchingString)) && (mismatches == bestMismatches)) {
						altContigReadStr = altContigReadStr.substring(altContigReadStr.indexOf('~')+1);
						altContigReadStr = altContigReadStr.replace('~', '\t');
						contigRead = samStringReader.getRead(altContigReadStr);
						
						// Filter this hit if it aligns past the end of the contig
						// Must use cigar length instead of read length, because the
						// the contig read bases are not loaded.
						if ((position + read.getReadLength()) <= contigRead.getCigar().getReadLength()) {
							HitInfo hit = new HitInfo(contigRead, position, strand, mismatches);
							bestHits.add(hit);
						}
					}
				}
			}
		}

		return bestHits;
	}

	private void updateMismatchAndEditDistance(SAMRecord read, CompareToReference2 c2r, SAMRecord origRead) {
		if (read.getAttribute(ORIGINAL_ALIGNMENT_TAG) != null) {
			int numMismatches = c2r.numMismatches(read);				
			int numIndelBases = SAMRecordUtils.getNumIndelBases(read);
			read.setAttribute("XM", numMismatches);
			read.setAttribute("NM", numMismatches + numIndelBases);
			read.setMappingQuality(calcMappingQuality(read, origRead));			
		}
	}

	private int calcMappingQuality(SAMRecord read, SAMRecord origRead) {
		int mapq = 0;
		
		// Need original read here because updated read has already had 0x04 flag unset.
		
		if ((origRead.getReadUnmappedFlag()) || (read.getMappingQuality() > 0)) {
			int contigQuality = (Integer) read.getAttribute("CONTIG_QUALITY_TAG");
			int quality = Math.min(contigQuality, this.maxMapq);
			int mismatchesToContig = (Integer) read.getAttribute(MISMATCHES_TO_CONTIG_TAG);
			quality -= mismatchesToContig * 5;
			mapq = Math.max(quality, 1);
		}
		
		return mapq;
	}
	
	private boolean isImprovedAlignment(SAMRecord read, SAMRecord orig, CompareToReference2 c2r) {
		
		boolean isImproved = false;
			
		Integer yr = orig.getIntegerAttribute("YR");
		if ((yr != null) && (yr == 1)) {
			// Original alignment was outside of target region list.
			// Calc updated edit distance to reference and compare to original
			
//			int origEditDistance = getOrigEditDistance(orig);
//			int updatedEditDistance = c2r.numMismatches(read) + getNumIndelBases(read);
			
			double origEditDistance = SAMRecordUtils.getOrigEditDistance(orig);
			double updatedEditDistance = c2r.numMismatches(read) + (1.5 * SAMRecordUtils.getNumIndels(read));
			
			isImproved = updatedEditDistance < origEditDistance;
		} else {
			isImproved = true;
		}
		
		return isImproved;
	}
	
	private void adjustForStrand(boolean readAlreadyReversed, SAMRecord read) {
		if ( ((!readAlreadyReversed) && (read.getReadNegativeStrandFlag())) ||
			 ((readAlreadyReversed) && (!read.getReadNegativeStrandFlag())) ){
			read.setReadString(reverseComplementor.reverseComplement(read.getReadString()));
			read.setBaseQualityString(reverseComplementor.reverse(read.getBaseQualityString()));
		}
	}
	
	private SAMRecord updateReadAlignment(SAMRecord contigRead,
			List<ReadBlock> contigReadBlocks, ReadPosition orig) {
		List<ReadBlock> blocks = new ArrayList<ReadBlock>();
		SAMRecord read = SAMRecordUtils.cloneRead(orig.getRead());

		read.setReferenceName(contigRead.getReferenceName());

		int contigPosition = orig.getPosition();
		int accumulatedLength = 0;

		// read block positions are one based
		// ReadPosition is zero based
		
		int totalInsertLength = 0;

		for (ReadBlock contigBlock : contigReadBlocks) {
			if ((contigBlock.getReadStart() + contigBlock.getReferenceLength()) >= orig
					.getPosition() + 1) {
				ReadBlock block = contigBlock.getSubBlock(accumulatedLength,
						contigPosition, read.getReadLength()
								- accumulatedLength);
				
				block.setReferenceStart(block.getReferenceStart() - totalInsertLength);
				
				// If this is an insert, we need to adjust the alignment start
				if ((block.getType() == CigarOperator.I) && (block.getLength() != 0)) {
					contigPosition = contigPosition - (contigBlock.getLength() - block.getLength());
					block.setReferenceStart(block.getReferenceStart() - (contigBlock.getLength() - block.getLength()));
//					block = contigBlock.getSubBlock(accumulatedLength,
//								contigPosition, read.getReadLength()
//								- accumulatedLength);
					
					totalInsertLength += block.getLength();
				}
				
				//TODO: Drop leading and trailing delete blocks

				// TODO: Investigate how this could happen
				if (block.getLength() != 0) {
					blocks.add(block);

					if (block.getType() != CigarOperator.D) {
						accumulatedLength += block.getLength();
					}

					if (accumulatedLength > read.getReadLength()) {
						throw new IllegalStateException("Accumulated Length: "
								+ accumulatedLength
								+ " is greater than read length: "
								+ read.getReadLength());
					}

					if (accumulatedLength == read.getReadLength()) {
						break;
					}
				}
			}
		}

		if (blocks.size() > 0) {
			
			// If we've aligned past the end of the contig resulting in a short Cigar
			// length, append additonal M to the Cigar			
			ReadBlock.fillToLength(blocks, read.getReadLength());
			int newAlignmentStart = blocks.get(0).getReferenceStart();
			String newCigar = ReadBlock.toCigarString(blocks);

			read.setCigarString(newCigar);
			read.setAlignmentStart(newAlignmentStart);
		} else {
			// TODO: Investigate how this could happen.
			read = null;
		}

		return read;
	}
	
	private boolean isFiltered(SAMRecord read) {
		return SAMRecordUtils.isFiltered(isPairedEnd, read);
	}
	
	private RealignmentWriter getRealignmentWriter(SAMFileWriter outputReadsBam, boolean isTightAlignment, String tempDir) {
		RealignmentWriter writer;
		
		if (isTightAlignment && isPairedEnd) {
			writer = new BetaPairValidatingRealignmentWriter(c2r, outputReadsBam, tempDir, minInsertLen, maxInsertLen);
		} else {
			writer = new SimpleRealignmentWriter(c2r, outputReadsBam, isTightAlignment);
		}
		
		return writer;
	}
	
	static class HitInfo {
		private SAMRecord record;
		private int position;
		private char strand;
		private int mismatches;
		
		public HitInfo(SAMRecord record, int position, char strand, int mismatches) {
			this.record = record;
			this.position = position;
			this.strand = strand;
			this.mismatches = mismatches;
		}

		public SAMRecord getRecord() {
			return record;
		}

		public int getPosition() {
			return position;
		}

		public char getStrand() {
			return strand;
		}
		
		public boolean isOnNegativeStrand() {
			return strand == '-';
		}
		
		public int getNumMismatches() {
			return mismatches;
		}
	}

}

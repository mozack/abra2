package abra;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Set;
import java.util.UUID;

import abra.SequenceUtil.MatchResult;


import net.sf.samtools.SAMRecord;

public class Contig {

	private UUID id;
	private String descriptor;
	private String sequence = "";
	private List<Node> nodes = new ArrayList<Node>();
	
	public Contig() {
		this.id = UUID.randomUUID();
	}
	
	/*
	public Contig(String descriptor, String sequence) {
		this();
		this.descriptor = descriptor;
		this.sequence = sequence;
	}
	*/
	
	public Contig(Contig orig) {
		this();
		this.sequence = orig.sequence;
		this.nodes.addAll(orig.nodes);
	}
	
	public int hashCode() {
		return id.hashCode();
	}
	
	public boolean equals(Object object) {
		Contig that = (Contig) object;
		return this.id.equals(that.id);
	}
	
	public String getDescriptor() {
		return descriptor;
	}
	
	public String getSequence() {
		return sequence;
	}
	
	public void append(Node node, String sequence) {
		nodes.add(node);
		this.sequence += sequence;
	}

	public void prependSequence(String prependDescriptor, String prefix) {
		descriptor += "_p_" + prependDescriptor;
		sequence = prefix + sequence;
	}
	
	public void setDescriptor(String descriptor) {
		this.descriptor = descriptor;
	}
	
	public List<Node> getNodes() {
		return Collections.unmodifiableList(nodes);
	}
	
	/**
	 * Returns information about each read that contributed to this Contig and
	 * the location of the read within the contig.  Only those reads that exactly
	 * match the contig are included. 
	 */
	public List<ReadPosition> getFilteredReadPositions(int allowedMismatchesFromContig, List<SAMRecord> allReads) {
		List<ReadPosition> readPositions = new ArrayList<ReadPosition>();
		int index = 0;
		
		// Loop through each node's starting reads
		for (Node node : nodes) {
			//for (SAMRecord read : node.getStartingReads()) {
			for (SAMRecord read : allReads) {
				// Does this read match the contig?
				String readSequence = read.getReadString();
				
				if (sequence.length() >= index + readSequence.length()) {
					String contigSubstring = sequence.substring(index, index+readSequence.length());
					
					//if (contigSubstring.equals(readSequence)) {
					MatchResult matchResult = SequenceUtil.isMatch(contigSubstring, readSequence, allowedMismatchesFromContig);
					if (matchResult.isMatch()) {
						// We have a match, record this read's position in the contig
						readPositions.add(new ReadPosition(read, index, matchResult.getNumMismatches()));
					}
				}
			}
			
			index += 1;
		}
		
		return readPositions;
	}
}

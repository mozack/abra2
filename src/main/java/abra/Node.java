/* Copyright 2013 University of North Carolina at Chapel Hill.  All rights reserved. */
package abra;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class Node {

	private Sequence sequence;
	private int count = 1;
//	private Edge[] toEdges = new Edge[0];
//	private Edge[] fromEdges = new Edge[0];
	
	private Node[] toNodes = new Node[0];
	private Node[] fromNodes = new Node[0];
	
//	private List<Edge> toEdges = new ArrayList<Edge>(1);
//	private List<Edge> fromEdges = new ArrayList<Edge>(1);
//	private Map<Node, Edge> toEdges = new HashMap<Node, Edge>(2);
//	private Map<Node, Edge> fromEdges = new HashMap<Node, Edge>(2);
	private String contributingRead = null;
	private boolean hasMultipleUniqueReads = false;
	
	public Node(String sequence) {
		this(new Sequence(sequence));
	}
	
	public Node(Sequence sequence) {
		this.sequence = sequence;
	}
	
	public String toString() {
		return count + "_" + sequence;
	}
	
	public int hashCode() {
		return sequence.hashCode();
	}
	
	public boolean equals(Object object) {
		Node that = (Node) object;
		return this.sequence.equals(that.sequence);
	}
	
	public void incrementCount() {
		count++;
	}
	
	public int getCount() {
		return count;
	}
	
	public Collection<Node> getToNodes() {
		return Arrays.asList(toNodes);
	}
	
	public Collection<Node> getFromNodes() {
		return Arrays.asList(fromNodes);
	}
	
/*
	public Collection<Edge> getToEdges() {
		return Arrays.asList(toEdges);
	}
	
	public Collection<Edge> getFromEdges() {
		return Arrays.asList(fromEdges);
	}
*/
	
	public Sequence getSequence() {
		return sequence;
	}
	
	public void addReadSequence(String sequence) {
		if (!hasMultipleUniqueReads) {
			if (contributingRead == null) {
				contributingRead = sequence;
			} else {
				if (!contributingRead.equals(sequence)) {
					hasMultipleUniqueReads = true;
					contributingRead = null;
				}
			}
		}
	}
				
	public boolean hasMultipleUniqueReads() {
		return hasMultipleUniqueReads;
	}
	
		
	private Node[] addNode(Node[] nodes, Node node) {
//		int minCapacity = edges.length + 1;
		//Edge[] newEdges = new Edge[edges.length + 1];
		Node[] newNodes = Arrays.copyOf(nodes, nodes.length + 1);
		newNodes[nodes.length] = node;
		
		return newNodes;
		
//        int oldCapacity = edges.length;
//        if (minCapacity > oldCapacity) {
//            Object oldData[] = edges;
//            int newCapacity = (oldCapacity * 3)/2 + 1;
//            if (newCapacity < minCapacity)
//                newCapacity = minCapacity;
//            // minCapacity is usually close to size, so this is a win:
//            elementData = Arrays.copyOf(elementData, newCapacity);
//        }
	}
	
	public void addToNode(Node to) {
		Node node = findNode(to);
		
		if (node == null) {
//			node = new Edge(this, to);
			toNodes = addNode(toNodes, to);
			to.updateFromEdges(this);
		} else {
//			node.incrementCount();
//			edge.incrementCount();
		}
	}
	
	public boolean isRootNode() {
		return fromNodes.length == 0;
	}
	
	private void updateFromEdges(Node node) {
		fromNodes = this.addNode(fromNodes, node);
//		fromEdges = addEdge(fromEdges, edge);
	}
	
	public boolean isSingleton() {
		return fromNodes.length == 0 && toNodes.length == 0;
	}
	
	//TODO: Smarter array allocation
	public void removeToNode(Node to) {
//		this.toEdges.remove(edge.getTo());
		int index = this.findToNodeIndex(to);
		if (index > -1) {
			Node[] newNodes = new Node[toNodes.length-1];
			System.arraycopy(toNodes, 0, newNodes, 0, index);
			if ((index) < newNodes.length) { 
				System.arraycopy(toNodes, index+1, newNodes, index, newNodes.length-index);
			}
			
			toNodes = newNodes;
		}
	}
	
	public void removeFromNode(Node node) {
//		this.fromEdges.remove(edge.getFrom());
		int index = this.findFromNodeIndex(node);
		if (index > -1) {
			Node[] newNodes = new Node[fromNodes.length-1];
			System.arraycopy(fromNodes, 0, newNodes, 0, index);
			if ((index) < newNodes.length) {
				System.arraycopy(fromNodes, index+1, newNodes, index, newNodes.length-index);
			}
			
			fromNodes = newNodes;
		}
	}
	
	private int findToNodeIndex(Node to) {
		for (int i=0; i<toNodes.length; i++) {
			if (toNodes[i].equals(to)) {
				return i;
			}
		}
		
		return -1;
	}
	
	private int findFromNodeIndex(Node from) {
		for (int i=0; i<fromNodes.length; i++) {
			if (fromNodes[i].equals(from)) {
				return i;
			}
		}
		
		return -1;
	}

	
	private Node findNode(Node to) {
		for (Node node : toNodes) {
			if (node.equals(to)) {
				return node;
			}
		}
		
		return null;
	}
	
	public void remove() {
		for (Node toNode : toNodes) {
			toNode.removeFromNode(this);
		}
		
		for (Node fromNode : fromNodes) {
			fromNode.removeToNode(this);
		}
	}
}

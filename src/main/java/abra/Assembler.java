package abra;

import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;

public class Assembler {

    private static final int TOO_MANY_PATHS_FROM_ROOT = -1;
    private static final int TOO_MANY_CONTIGS = -2;
    private static final int STOPPED_ON_REPEAT = -3;
    private static final int TOO_MANY_NODES = -4;

    private static int MAX_CONTIG_SIZE = 5000;
    //private static int MAX_READ_LENGTH = 1000;

    private static int MAX_FREQUENCY = 32766;
    private static int MAX_QUAL_SUM = 255;
    //private static int MAX_KMER_LEN = 201;
    private static int MAX_SAMPLES = 16;
    private static int MIN_BASE_QUALITY = 8;

    static class Node {
        static int nextNodeId = 1;
        int nodeId;
        String kmer;
        String seq = null;
        String contributingRead;
        int contributingStrand;
        int sampleFrequency [] = new int[MAX_SAMPLES];
        int frequency = 1;
        boolean hasMultipleUniqueReads;
        int [] qualSums;
        LinkedList<Node> toNodes = new LinkedList<>();
        LinkedList<Node> fromNodes = new LinkedList<>();
        boolean isCondensed = false;
        boolean isFiltered = false;
        boolean isRoot = false;

        Node() {
            // empty
        }

        Node(char sampleId, String kmer, String contributingRead, int strand, String quals, int khmerSize) {
            this.nodeId = nextNodeId++;
            this.kmer = kmer;
            this.contributingRead = contributingRead;
            this.sampleFrequency[sampleId-1] = 1;
            this.hasMultipleUniqueReads = false;
            this.contributingStrand = strand;
            this.qualSums = new int[khmerSize];
            for (int i = 0; i < khmerSize; i++) {
                qualSums[i] = phred33(quals.charAt(i));
            }
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;
            Node node = (Node) o;
            return Objects.equals(kmer, node.kmer);
        }

        @Override
        public int hashCode() {
            return Objects.hash(kmer);
        }

        @Override
        public String toString() {
            return "Node{" +
                    "kmer='" + kmer + '\'' +
                    '}';
        }
    }

    private static void pruneGraph(HashMap<String, Node> nodes, int minNodeFreq, double minEdgeRatio, int kmerSize,
                                   int minBaseQuality) {
        //System.out.printf("  Initial nodes before pruning step 0: %d\n", nodes.size());
        // First prune kmers that do not reach base quality sum threshold

        Iterator<String> it = nodes.keySet().iterator();
        while (it.hasNext()) {
            String key = it.next();
            Node node = nodes.get(key);
            if (!isBaseQualityGood(node, kmerSize, minBaseQuality)) {
                removeNodeAndCleanup(node, it);
            }
        }

        // TODO Logger.debug
        //System.out.printf("  Remaining nodes after pruning step 1: %d\n", nodes.size());

        // Now go back through and ensure that each node reaches minimum frequency threshold.
        if (minNodeFreq > 1) {
            it = nodes.keySet().iterator();
            while (it.hasNext()) {
                String key = it.next();
                Node node = nodes.get(key);
                if (node.frequency < minNodeFreq || (!node.hasMultipleUniqueReads)) {
                    it.remove();
                }
            }
        }

        //System.out.printf("  Remaining nodes after pruning step 2: %d\n", nodes.size());

        pruneLowFrequencyEdges(nodes, minEdgeRatio);

        // Final pass through cleaning up nodes that are unreachable
        it = nodes.keySet().iterator();
        while (it.hasNext()) {
            String key = it.next();
            Node node = nodes.get(key);
            if (node.toNodes.size() == 0 && node.fromNodes.size() == 0) {
                it.remove();
            }
        }

        //System.out.printf("  Remaining nodes after edge pruning: %d\n", nodes.size());
    }

    private static void removeNodeAndCleanup(Node node, Iterator<String> it) {
        // Remove node from "from" lists
        for (Node toNode:node.toNodes) {
            toNode.fromNodes = removeNodeFromList(node, toNode.fromNodes);
        }

        // Remove node from "to" lists
        for (Node fromNode:node.fromNodes) {
            fromNode.toNodes = removeNodeFromList(node, fromNode.toNodes);
        }

        // Remove node from map
        it.remove();
        node.toNodes = new LinkedList<>();
        node.fromNodes = new LinkedList<>();
    }

    private static void pruneLowFrequencyEdges(HashMap<String, Node> nodes, double minEdgeRatio) {

        long removedEdgeCount = 0;

        for (Map.Entry<String, Node> entry : nodes.entrySet()) {
            Node currNode = entry.getValue();
            ////////////////////////////////////////////////
            // Check to node list for low frequency edges

            // Calculate total outgoing "edge" frequency
            int toNodeTotalFreq = 0;

            int[] perSampleTotalFreq = new int[MAX_SAMPLES]; // array of zero values

            for (Node toNode: currNode.toNodes) {
                // Using node frequency as proxy for edge frequency here...
                toNodeTotalFreq = toNodeTotalFreq + toNode.frequency;
                for (int i=0; i < MAX_SAMPLES; i++) {
                    perSampleTotalFreq[i] += toNode.sampleFrequency[i];
                }
            }

            // Identify edges to prune
            LinkedList<Node> toNodesToRemove = new LinkedList<>();
            for (Node toNode : currNode.toNodes) {
                boolean exceedsMinRatio = isMinEdgeRatioReached(perSampleTotalFreq, toNode, minEdgeRatio);

                if (!exceedsMinRatio) {
                    toNodesToRemove.addLast(toNode);
                }
            }

            // Remove edges from lists (node remains in HashMap)
            for (Node nodeToRemove : toNodesToRemove) {
                // Remove edges in each direction
                nodeToRemove.fromNodes = removeNodeFromList(currNode, nodeToRemove.fromNodes);
                currNode.toNodes = removeNodeFromList(nodeToRemove, currNode.toNodes);
                removedEdgeCount += 1;
                // TODO logger.debug
                //System.out.println("Removing kmer:" + nodeToRemove);
            }

            ////////////////////////////////////////////////
            // Check from node list for low frequency edges

            // Calculate total outgoing "edge" frequency
            int fromNodeTotalFreq = 0;
            perSampleTotalFreq = new int[MAX_SAMPLES];

            for (Node fromNode: currNode.fromNodes) {
                fromNodeTotalFreq =  fromNodeTotalFreq + fromNode.frequency;
                for (int i=0; i < MAX_SAMPLES; i++) {
                    perSampleTotalFreq[i] += fromNode.sampleFrequency[i];
                }
            }

            // Identify edges to prune
            LinkedList<Node> fromNodesToRemove = new LinkedList<>();
            for (Node fromNode: currNode.fromNodes) {
                boolean exceedsMinRatio = isMinEdgeRatioReached(perSampleTotalFreq, fromNode, minEdgeRatio);

                if (!exceedsMinRatio) {
                    fromNodesToRemove.addLast(fromNode);
                }
            }

            // Remove edges from lists (node remains in HashMap)
            for (Node nodeToRemove : fromNodesToRemove) {
                // Remove edges in each direction
                nodeToRemove.toNodes = removeNodeFromList(currNode, nodeToRemove.toNodes);
                currNode.fromNodes = removeNodeFromList(nodeToRemove, currNode.fromNodes);
                removedEdgeCount +=1;
                //System.out.println("  Removed from node: " + nodeToRemove);
            }
        }

        //System.err.println("  Pruned " + removedEdgeCount + " edges");
    }

    private static LinkedList<Node> removeNodeFromList(Node node, LinkedList<Node> nodeList) {
        nodeList.remove(node); // not always there
        return nodeList;
    }

    private static boolean isMinEdgeRatioReached(int[] perSampleTotalFreq, Node node, double minEdgeRatio) {
        for (int i=0; i < MAX_SAMPLES; i++) {


            if (perSampleTotalFreq[i] > 0 &&
                    ((double) node.sampleFrequency[i] / (double) perSampleTotalFreq[i] >= minEdgeRatio)) {
                return true;
            }

//            if (0 != perSampleTotalFreq[i] && node.sampleFrequency[i] != 0)
//                System.out.printf("  sample: %d, freq: %d, total_freq: %d, ratio: %f\n", i, node.sampleFrequency[i], perSampleTotalFreq[i], (double) node.sampleFrequency[i] / (double) perSampleTotalFreq[i]);

        }
        return false;
    }

    private static boolean isBaseQualityGood(Node node, int kmerSize, int minBaseQuality) {
        for (int i=0; i < kmerSize; i++) {
            if (node.qualSums[i] < minBaseQuality) {
                return false;
            }
        }
        return true;
    }

    private
    static void buildGraph2(String input, HashMap<String, Node> nodes, int readLength, int maxNodes, int kmerSize) {
        int recordLength = readLength * 2 + 2;
        int numRecords = input.length() / recordLength;
        int record = 0;

        while (record < numRecords && nodes.size() < maxNodes) {
            int startIndex = record * recordLength;

            char sampleId = input.charAt(startIndex);
            char strandChar = input.charAt(startIndex + 1);
            int strand = 0;
            if (strandChar == '0') {
                // no-op strand = 0;
            } else if (strandChar == '1')
                strand = 1;
            else {
                System.err.println("Initial char in input invalid. Exiting.");
                System.exit(-1);
            }

            String readStr = input.substring(startIndex + 2, startIndex + 2 + readLength);
            int qualStrStart = startIndex + 2 + readLength;
            String qualStr = input.substring(qualStrStart, qualStrStart+readLength);

            addToGraph(nodes, sampleId, strand, readStr, qualStr, kmerSize, readLength);

            record++;
        }
    }

    private static void addToGraph(HashMap<String, Node> nodes, char sampleId, int strand, String sequence, String qual,
                                     int kmerSize, int readLength) {
        Node prev = null;

        for (int i=0; i <= readLength-kmerSize; i++) {

            if (includeKmer(sequence, qual, i, kmerSize)) {
                String kmer = getKmer(i, sequence, kmerSize);
                String kmerQual = getKmer(i, qual, kmerSize);

                Node curr;
                if (nodes.containsKey(kmer)) {
                    curr = nodes.get(kmer);
                    incrementNodeFreq(sampleId, curr, sequence, strand, kmerQual, kmerSize);
                } else {
                    curr = new Node(sampleId, kmer, sequence, strand, kmerQual, kmerSize);
                    nodes.put(kmer, curr);
                }

                if (prev != null) {
                    linkNodes(prev, curr);
                }

                prev = curr;
            } else {
                prev = null;
            }
        }
    }

    private static void linkNodes(Node fromNode, Node toNode) {
        if (!fromNode.toNodes.contains(toNode)) {
            fromNode.toNodes.addFirst(toNode);
        }

        if (!toNode.fromNodes.contains(fromNode)) {
            toNode.fromNodes.addFirst(fromNode);
        }
    }

    private
    static void incrementNodeFreq(char sampleId, Node node, String sequence, int strand, String kmerQual, int kmerSize) {
        if (node.frequency < MAX_FREQUENCY - 1)
            node.frequency++;

        if (node.sampleFrequency[sampleId-1] < MAX_FREQUENCY - 1)
            node.sampleFrequency[sampleId-1]++;

        if (!node.hasMultipleUniqueReads &&
                (!node.contributingRead.equals(sequence) || node.contributingStrand != strand)) {
            node.hasMultipleUniqueReads = true;
        }

        for (int i=0; i < kmerSize; i++) {
            int phred33Qual = phred33(kmerQual.charAt(i));
            if (node.qualSums[i] + phred33Qual < MAX_QUAL_SUM) {
                node.qualSums[i] += phred33Qual;
            } else {
                node.qualSums[i] = MAX_QUAL_SUM;
            }
        }
    }

    private static String getKmer(int i, String sequence, int kmerSize) {
        return sequence.substring(i, i+kmerSize);
    }

    private static boolean includeKmer(String sequence, String qual, int idx, int kmerSize) {
        for (int i = idx; i < idx+kmerSize; i++) {
            // Discard kmers with ambiguous bases
            if (sequence.charAt(i) == 'N')
                return false;

            // Discard kmers with low base qualities
            if (phred33(qual.charAt(i)) < MIN_BASE_QUALITY)
                return false;
        }
        return true;
    }

    private static int phred33(char ch) {
        return (int) ch - '!';
    }

    // NOTE: From nodes are invalid after this step!!!
    private static void condenseGraph(HashMap<String, Node> nodes) {
        int nodesCondensed = 1;

        for (Map.Entry<String, Node> entry : nodes.entrySet()) {
            Node node = entry.getValue();

            // Starting point 0 or >1 incoming edges or previous node with multiple outgoing edges and curr node has 1 outgoing edge
            if ((!hasOneIncomingEdge(node) || prevHasMultipleOutgoingEdges(node)) && hasOneOutgoingEdge(node)) {
                Node next = node.toNodes.get(0); // because hasOneOutgoingEdge is true

                if (hasOneIncomingEdge(next)) {
                    List<Node> last = next.toNodes;

                    StringBuilder seq = new StringBuilder(node.kmer.length());

                    seq.append(node.kmer.charAt(0));

                    while (next != null && hasOneIncomingEdge(next) && nodesCondensed < MAX_CONTIG_SIZE) {
                        last = next.toNodes;

                        if (next.toNodes != null && next.toNodes.size() > 0) {

                            seq.append(next.kmer.charAt(0));
                        } else {
                            // end of path, copy entire kmer
                            seq.append(next.kmer);
                        }

                        Node temp = null;
                        if (hasOneOutgoingEdge(next)) {
                            temp = next.toNodes.get(0);
                        }

                        next.isFiltered = true;

                        next = temp;

                        nodesCondensed += 1;
                    }

                    // Update node
                    node.seq = seq.toString();
                    node.isCondensed = true;
                    //System.out.println("IDX: " + idx + " SEQ:" + seq + " kmer:" + node.kmer + " FREQ:" + node.frequency);

                    // Free original toNodes list
                    node.toNodes = null;

                    // copy last toNodes to this node
                    if (last == null || last.size() == 0) {
                        // no-op node.toNodes = null;
                    } else {
                        node.toNodes = new LinkedList<>(last);
                        //last = null;
                    }
                }
            }
        }
    }

    private static boolean hasOneIncomingEdge(Node node) {
        return (node.fromNodes != null && node.fromNodes.size() == 1);
    }

    private static boolean hasOneOutgoingEdge(Node node) {
        return (node.toNodes != null && node.toNodes.size() == 1);
    }

    private static boolean prevHasMultipleOutgoingEdges(Node node) {
        boolean prevBifurcates = false;
        if (hasOneIncomingEdge(node)) {
            Node prev = node.fromNodes.get(0);
            if (prev != null && prev.toNodes != null && prev.toNodes.size() >= 2)
                prevBifurcates = true;
        }
        return prevBifurcates;
    }

    private static LinkedList<Node> identifyRootNodes(HashMap<String, Node> nodes, LinkedList<Node> list) {
//        if (list.size() != 12) {
//            while(list.size() != 12)
//                list.add(new Node());
//        }

        nodes.forEach((key, node) -> {
            // TODO
            // hashmap iterator order changes the output, make order in which root nodes are identified unimportant
            // TODO
            if (isRoot(node)) {
                node.isRoot = true;

//                if (node.kmer.equals("ACTAGGCTGAAAA")) {
//                    list.set(0, node);
//                } else if (node.kmer.equals("CATCTCTTAATCC")){
//                    list.set(1, node);
//                } else if (node.kmer.equals("AGAATTCTAATAA")){
//                    list.set(2, node);
//                } else if (node.kmer.equals("GTTTGTAAGGCCT")){
//                    list.set(3, node);
//                } else if (node.kmer.equals("GTGGCAAACATCC")){
//                    list.set(4, node);
//                } else if (node.kmer.equals("CAGAATTCCAATA")){
//                    list.set(5, node);
//                } else if (node.kmer.equals("GGCCAGCCTCATA")){
//                    list.set(6, node);
//                } else if (node.kmer.equals("CTCTTTGGTTATT")){
//                    list.set(7, node);
//                } else if (node.kmer.equals("CTGGCTATTCCTT")){
//                    list.set(8, node);
//                } else if (node.kmer.equals("CTTCTTAGCACTT")){
//                    list.set(9, node);
//                } else if (node.kmer.equals("AAGAATAACCAAA")){
//                    list.set(10, node);
//                } else if (node.kmer.equals("GAATTCTGCGTTT")){
//                    list.set(11, node);
//                }

                list.addFirst(node);
            }
        });
        //System.out.println("Root Nodes:" + list.size());
        return list;
    }

    private static boolean isRoot(Node node) {
        if (node.fromNodes.size() == 0) {
            return true;
        } else {
            // Identify nodes that point to themselves with no other incoming edges.
            // This will be cleaned up during contig building.
            if (node.fromNodes.size() == 1 && node.fromNodes.get(0).kmer.equals(node.kmer)) {
                return true;
            }
        }
        return false;
    }

    static class Contig {
        LinkedList<String> fragments = new LinkedList<>();
        Node currNode;
        HashSet<Integer> visitedNodes = new HashSet<>();
        double score = 0;
        int realSize = 0;
        boolean isRepeat = false;

        public Contig(Node node) {
            this.currNode = node;
        }

        public Contig(Contig o) {
            this.fragments = new LinkedList<>(o.fragments);
            this.currNode = o.currNode;
            this.visitedNodes = new HashSet<>(o.visitedNodes);
            this.score = o.score;
            this.realSize = o.realSize;
            this.isRepeat = o.isRepeat;
        }

        @Override
        public String toString() {
            return "Contig{" +
                    "currNode=" + currNode +
                    ", score=" + score +
                    ", realSize=" + realSize +
                    ", isRepeat=" + isRepeat +
                    ", visitedNodes=" + visitedNodes +
                    ", fragments=" + fragments +
                    '}';
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            Contig contig = (Contig) o;

            if (Double.compare(contig.score, score) != 0) return false;
            if (realSize != contig.realSize) return false;
            if (isRepeat != contig.isRepeat) return false;
            if (fragments != null ? !fragments.equals(contig.fragments) : contig.fragments != null) return false;
            if (currNode != null ? !currNode.equals(contig.currNode) : contig.currNode != null) return false;
            return visitedNodes != null ? visitedNodes.equals(contig.visitedNodes) : contig.visitedNodes == null;
        }

        @Override
        public int hashCode() {
            int result;
            long temp;
            result = fragments != null ? fragments.hashCode() : 0;
            result = 31 * result + (currNode != null ? currNode.hashCode() : 0);
            result = 31 * result + (visitedNodes != null ? visitedNodes.hashCode() : 0);
            temp = Double.doubleToLongBits(score);
            result = 31 * result + (int) (temp ^ (temp >>> 32));
            result = 31 * result + realSize;
            result = 31 * result + (isRepeat ? 1 : 0);
            return result;
        }
    }

    static int buildContigs(
            Node root,
            AtomicInteger contigCount,
            String prefix,
            int maxPathsFromRoot,
            int maxContigs,
            boolean stopOnRepeat, // aka truncateOnRepeat
            StringBuilder contigStr,
            PriorityQueue<Double> contigScores,
            int minContigLength,
            int kmerSize) {
        int status = 0;

        Stack<Contig> contigsToOutput = new Stack<>();
        Stack<Contig> poppedContigs = new Stack<>();
        Stack<Contig> contigs = new Stack<>();
        Contig rootContig = new Contig(root);
        contigs.push(rootContig);

        LinkedList<String> allContigFragments = new LinkedList<>();

        int pathsFromRoot = 1;

        Logger.debug("Building contig for root node: " + root);

        while (contigs.size() > 0 && status == 0) {
            // Get contig from stack
            Contig contig = contigs.peek();

            if (isNodeVisited(contig, contig.currNode)) {
                // We've encountered a repeat
                contig.isRepeat = true;
                poppedContigs.push(contig);
                contigs.pop();

                if (stopOnRepeat) {
                    status = STOPPED_ON_REPEAT;
                }
            }
            else if (contig.currNode.toNodes == null || contig.currNode.toNodes.size() == 0 ||
                    contig.realSize >= (MAX_CONTIG_SIZE - 1)) {
                // We've reached the end of the contig.
                // Append entire current node.

                appendToContig(contig, allContigFragments, true, kmerSize);

                // Now, write the contig
                contigsToOutput.push(contig);
                updateContigScores(contigScores, contig.score);

                contigs.pop();
            }
            else {
                // Append first base from current node
                appendToContig(contig, allContigFragments, false, kmerSize);

                if (contig.realSize >= MAX_CONTIG_SIZE) {
                    System.err.println("Max contig size exceeded at node:"+contig.currNode.kmer+" contig size:"+contig.realSize);

                    // TODO: Provide different status
                    status = TOO_MANY_CONTIGS;
                    break;
                }

                contig.visitedNodes.add(contig.currNode.nodeId);

                // Count total edges
                int totalEdgeCount = 0;
                List<Node> toNodeList = contig.currNode.toNodes;

                for (Node toNode: toNodeList) {
                    totalEdgeCount = totalEdgeCount + toNode.frequency;
                }

                // Move current contig to next "to" node.
                toNodeList = contig.currNode.toNodes;
                contig.currNode = toNodeList.get(0);
                pathsFromRoot++;

                double prevContigScore = contig.score;

                double log10TotalEdgeCount = 0;

                // Update contig score if there is fork here.
                // Otherwise, we would just multiply by 1.
                if (toNodeList.size() > 1) {
                    log10TotalEdgeCount = Math.log10(totalEdgeCount);
                    contig.score = contig.score + Math.log10(contig.currNode.frequency) - log10TotalEdgeCount;
                    //System.out.println("CONTIG SCORE:" + contig.score + " FREQUENCY:" + contig.currNode.frequency + " EDGE COUNT:" + log10TotalEdgeCount + " TO NODES:" + toNodeList.size());
                }

                if (!isContigScoreOk(contigScores, contig.score)) {
                    poppedContigs.push(contig);
                    contigs.pop();
                }

                // If there are multiple "to" nodes, branch the contig and push on stack
                Iterator<Node> toLinkedNodeIter = toNodeList.iterator();
                Node toLinkedNode;
                toLinkedNodeIter.next();
                if (toLinkedNodeIter.hasNext())
                    toLinkedNode = toLinkedNodeIter.next();
                else
                    toLinkedNode = null;
                while (toLinkedNode != null) {

                    double contigBranchScore = prevContigScore + Math.log10(toLinkedNode.frequency) - log10TotalEdgeCount;

                    if (isContigScoreOk(contigScores, contigBranchScore)) {
                        Contig contigBranch = new Contig(contig);
                        contigBranch.currNode = toLinkedNode;
                        contigBranch.score = contigBranchScore;
                        contigs.push(contigBranch);
                    }

                    if (toLinkedNodeIter.hasNext())
                        toLinkedNode = toLinkedNodeIter.next();
                    else
                        toLinkedNode = null;
                    pathsFromRoot++;
                }
            }

            if (contigCount.get() >= maxContigs) {
                status = TOO_MANY_CONTIGS;
            }

            if (pathsFromRoot >= maxPathsFromRoot) {
                status = TOO_MANY_PATHS_FROM_ROOT;
            }
        }

        if (status == 0) {

            while (contigsToOutput.size() > 0) {
                Contig contig = contigsToOutput.peek();

                if (isContigScoreOk(contigScores, contig.score)) {
                    outputContig(contig, contigCount, prefix, contigStr, minContigLength);
                }

                contigsToOutput.pop();
                freeContig(contig);
            }
        }

        // Cleanup stranded contigs in case processing stopped.
        contigsToOutput.clear();

        contigs.clear();

        poppedContigs.clear();

        allContigFragments.clear();

        return status;
    }

    private static void freeContig(Contig contig) {
        contig.visitedNodes.clear();
        contig.fragments.clear();
    }

    private static void updateContigScores(PriorityQueue<Double> contigScores, double score) {
        if (contigScores.size() == 128 && score >= contigScores.peek()) {
            contigScores.remove();
            contigScores.add(score);
        }
        else if (contigScores.size() < 128) {
            //System.out.println("Adding contig score:" + score);
            contigScores.add(score);
        }
    }

    // Append curr node to contig
    private
    static void appendToContig(Contig contig, LinkedList<String> allContigFragments, boolean entireKmer, int kmerSize) {

        if (contig.currNode.isCondensed) {
            // Add condensed node sequence to fragment vector
            contig.fragments.addLast(contig.currNode.seq);
            contig.realSize += contig.currNode.seq.length();
        } else {

            if (!entireKmer) {
                contig.fragments.addLast(contig.currNode.kmer.substring(0,1));
                contig.realSize += 1;
            } else {
                String fragment = contig.currNode.kmer;
                //System.out.println("Appending fragment:" + fragment);
                contig.fragments.addLast(fragment);
                contig.realSize += kmerSize;
                allContigFragments.addLast(fragment);
            }
        }
    }

    private static boolean isNodeVisited(Contig contig, Node node) {
        return contig.visitedNodes.contains(node.nodeId);
    }

    private static void outputContig(Contig contig, AtomicInteger contigCount, String prefix, StringBuilder contigs,
                                     int minContigLength) {
        if (contig.realSize >= minContigLength) {
            if (!contig.isRepeat) {

                String idLine = String.format(">%s_%d_%f\n", prefix, contigCount.getAndIncrement(), contig.score);
                contigs.append(idLine);

                for(String s: contig.fragments) {
                    contigs.append(s);
                }

                contigs.append("\n");
                //System.out.println("Output contig:" + contig + " CONTIG COUNT::::" + contigCount);
            }
        }
    }

    // Return true if minimum contig score exceeded (where min = lowest of 128 top scoring contigs)
    // Update priority queue as needed
    private static boolean isContigScoreOk(PriorityQueue<Double> contigScores, double score) {

        //System.out.printf("contig_scores size: %d\n", contigScores.size());

        // 128 = max contigs
        // TODO: parameterize
        if (contigScores.size() == 128) {
            double minScore = contigScores.peek();
            return score >= minScore;
        } else if (contigScores.size() < 128) {
            return true;
        } else {
            System.err.printf("ERROR: Invalid contig score size: %d\n", contigScores.size());
        }
        return false;
    }

    // for http://www.webgraphviz.com/ or similar
    public static void dumpGraph(Map<String, Node> nodes, String filename) {
        try {
            System.out.printf("Filename: %s\n", filename);
            FileWriter fileWriter = new FileWriter(filename);
            PrintWriter fp = new PrintWriter(fileWriter);

            // Output edges
            fp.printf("digraph vdjer {\n//\tEdges\n");
            nodes.forEach((key, currNode) -> {
                if (!currNode.isFiltered) {
                    if (currNode.toNodes != null) {
                        currNode.toNodes.forEach((tonode) -> {
                            fp.printf("\tv_%d -> v_%d\n", currNode.nodeId, tonode.nodeId);
                        });
                    }
                }
            });

            int numVertices = 0;
            int numCondensed = 0;

            // Output vertices
            fp.printf("//\tVertices\n");
            for (Map.Entry<String, Node> entry : nodes.entrySet()) {
                Node currNode = entry.getValue();

                // Skip orphans
                if (!currNode.isFiltered) {
                    if (currNode.isCondensed) {
                        if (currNode.isRoot) {
                            fp.printf("\tv_%d [label=\"%s\",shape=box,color=green] %d\n", currNode.nodeId, currNode.seq, currNode.frequency);
                        } else {
                            fp.printf("\tv_%d [label=\"%s\",shape=box,color=blue] %d\n", currNode.nodeId, currNode.seq, currNode.frequency);
                        }
                        numCondensed += 1;
                    } else {
                        if (currNode.isRoot) {
                            fp.printf("\tv_%d [label=\"%s\",shape=box,color=red]\n", currNode.nodeId, currNode.kmer);
                        } else {
                            fp.printf("\tv_%d [label=\"%s\",shape=box] %d\n", currNode.nodeId, currNode.kmer, currNode.frequency);
                        }
                    }

                    numVertices += 1;
                }
            }

            fp.printf("}\n");

            fp.flush();
            fp.close();

            System.out.printf("Num traversable vertices: %d\n", numVertices);
            System.out.printf("Num condensed vertices: %d\n", numCondensed);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    /**
     *
     * @return
     */
    public static String assemble(String input,
                                  String output,
                                  String prefix,
                                  boolean truncateOnRepeat,
                                  int maxContigs,
                                  int maxPathsFromRoot,
                                  int inputReadLength,
                                  int inputKmerSize,
                                  int minNodeFreq,
                                  int minBaseQuality,
                                  double minEdgeRatio,
                                  boolean debug,
                                  int maxNodes) {

        int readLength = inputReadLength;
        int kmerSize = inputKmerSize;
        int minContigLength = readLength + 1;
        //TODO: Parameterize mcl - shorter for unaligned region?

        if (debug) {
            System.err.printf("Abra non-JNI entry point, prefix: %s, read_length: %d, kmer_size: %d, min_node_freq: %d, min_base_qual: %d, min_edge_ratio %f, max_nodes: %d\n",
                    prefix, readLength, kmerSize, minNodeFreq, minBaseQuality, minEdgeRatio, maxNodes);
        }
        long startTime = System.currentTimeMillis();

        HashMap<String, Node> nodes = new HashMap<>();

        buildGraph2(input, nodes, readLength, maxNodes, kmerSize);

        int status = 0;
        if (nodes.size() >= maxNodes) {
            status = TOO_MANY_NODES;
            Logger.debug("Graph too complex for region:");
        }

        boolean isUnalignedRegion = !truncateOnRepeat;
        pruneGraph(nodes, minNodeFreq, minEdgeRatio, kmerSize, minBaseQuality);

        LinkedList<Node> rootNodes = new LinkedList<>();
        if (status != TOO_MANY_NODES) {
            rootNodes = identifyRootNodes(nodes, rootNodes);
        }

        condenseGraph(nodes);

        AtomicInteger contigCount = new AtomicInteger(0);
        //boolean truncateOutput = false; // unused cpp variable

        StringBuilder contigStr = new StringBuilder();
        PriorityQueue<Double> contigScores = new PriorityQueue<>();

        for (Node rootNode: rootNodes) {
            status = buildContigs(rootNode, contigCount, prefix, maxPathsFromRoot, maxContigs, truncateOnRepeat,
                    contigStr, contigScores, minContigLength, kmerSize);

            switch (status) {
                case TOO_MANY_CONTIGS:
                    System.err.printf("TOO_MANY_CONTIGS: %s\n", prefix);
                    break;
                case STOPPED_ON_REPEAT:
                    if (debug)
                        System.err.printf("STOPPED_ON_REPEAT: %s\n", prefix);
                    break;
                case TOO_MANY_PATHS_FROM_ROOT:
                    System.err.printf("TOO_MANY_PATHS_FROM_ROOT: %s - %s\n", prefix, rootNodes.get(0).kmer);
                    break;
            }

            // If too many contigs or abort due to repeat, break out of loop and truncate output.
            if ((status == TOO_MANY_CONTIGS) || (status == STOPPED_ON_REPEAT)) {
                //truncateOutput = true;
                break;
            }
        }

        rootNodes.clear();
        nodes.clear();

        long stopTime = System.currentTimeMillis();

        if (debug)
          System.err.printf("Done assembling(%d ms): %s contigs=%d \n %s \n", (stopTime-startTime), output, contigCount.get(), contigStr.toString());

        if (status == 0 || status == TOO_MANY_PATHS_FROM_ROOT)
            return contigStr.toString();
        else if (status == STOPPED_ON_REPEAT)
            return "<REPEAT>";
        else
            return "<ERROR>";
    }
}
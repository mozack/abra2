package abra;

import org.junit.Test;

import java.util.ArrayList;
import java.util.Random;

import static org.junit.Assert.assertEquals;

public class ReverseComplementorTest {
//C->G, G->C, T->A, A->T

    @Test
    public void reverseComplement() {
        ReverseComplementor r = new ReverseComplementor();
        String x = r.reverseComplement("AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT");
        assertEquals("AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT", x);
    }

    @Test
    public void reverseComplementNotAllDNA() {
        ReverseComplementor r = new ReverseComplementor();
        String x = r.reverseComplement("NAATGANNN");
        assertEquals("NNNTCATTN", x);
        System.out.println("REVERSE:" + x);
    }

    @Test
    public void reverseComplementTimeTest() {
        ReverseComplementor r = new ReverseComplementor();
        ArrayList<String> dnaList = new ArrayList<String>();
        int nStringToReverse = 1000000;

        for (int i=0; i < nStringToReverse; i++) {
            String s = makeRandomDNA();
            dnaList.add(s);
        }

        long start = System.currentTimeMillis();

        for (int i=0; i < nStringToReverse; i++)
            r.reverseComplement(dnaList.get(i));

        long stop = System.currentTimeMillis();

        System.out.println("Elapsed time to reverse (ms):" + (stop-start));
    }

    public static String makeRandomDNA() {
        Random generator = new Random(0);
        int numSymbols = (int)(generator.nextDouble() * 41) + 40;  // random integer between 40 and 80 inclusive

        // Make a random DNA sequence with numSymbols symbols
        String DNAletters = new String("ACTG");
        StringBuilder result = new StringBuilder();

        for (int i = 1; i <= numSymbols; i++) {
            // Pick a random letter from "ACTG" string
            char symbol = DNAletters.charAt((int)(generator.nextDouble()*4));
            // concatenate next random symbol on to result
            result.append(symbol);
        }

        return result.toString();
    }
}
package abra.utils;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.zip.GZIPInputStream;

public class DbSnpSearcher {

	private SortedSet<Variant> variants = new TreeSet<Variant>();
	
	public void run(String dbSnpFile, String variantFile) throws Exception {
		load(dbSnpFile);
		
		BufferedReader reader = new BufferedReader(new FileReader(variantFile));
		
		String line = reader.readLine();
		
		int matchCount = 0;
		
		while (line != null) {
			String[] fields = line.split("\t");
			String chr = "chr" + fields[1];
			int pos = Integer.parseInt(fields[2]);
			int length = Math.max(Integer.parseInt(fields[4]), Integer.parseInt(fields[5]));
			
//			int fuzz = 100 + length;
			int fuzz = 100;
			
			Variant floor = new Variant("", chr, pos-fuzz, length);
			Variant ceil = new Variant("", chr, pos+fuzz, length);
			
			Set<Variant> matches = variants.subSet(floor, ceil);
			
			StringBuffer matchStr = new StringBuffer();
			
			boolean isMatch = false;
			
			for (Variant match : matches) {
//				if (match.length == length) {
				if (Math.abs(match.length - length) <= (length/3)+1) {
					matchStr.append(match.toString());
					matchStr.append(',');
					isMatch = true;
				}
			}
			
			if (isMatch) {
				matchCount++;
			}
			
			System.out.println(line.trim() + "\t" + matchStr);
			line = reader.readLine();
		}
		
		reader.close();
		
		System.out.println("match count: " + matchCount);
	}
	
	private void load(String dbSnpFile) throws FileNotFoundException, IOException {
		
		
		GZIPInputStream gzip = new GZIPInputStream(new FileInputStream(dbSnpFile));
		BufferedReader br = new BufferedReader(new InputStreamReader(gzip));
		
		String line = br.readLine();
		
		int count = 0;
		
		while (line != null) {
			String[] fields = line.split("\t");
			String chr = fields[1];
			int pos = Integer.parseInt(fields[2]);
			String name = fields[4];
			int length = fields[8].length();
			
			Variant var = new Variant(name, chr, pos, length);
			
			variants.add(var);
			
			line = br.readLine();
			count++;
		}
		
		br.close();
		gzip.close();
		
//		System.err.println("" + count + " variants loaded");
	}
	
	static class Variant implements Comparable<Variant> {

		String name;
		String chr;
		int pos;
		int length;
		
		public Variant(String name, String chr, int pos, int length) {
			this.name = name;
			this.chr = chr;
			this.pos = pos;
			this.length = length;
		}
				
		@Override
		public int compareTo(Variant that) {
			// TODO Auto-generated method stub
			int compare = this.chr.compareTo(that.chr);
			
			if (compare == 0) {
				compare = this.pos - that.pos;
			}
			
			if (compare == 0) {
				compare = this.length - that.length;
			}
			
			if (compare == 0) {
				compare = this.name.compareTo(that.name);
			}
			
			return compare;
		}
		
		@Override
		public int hashCode() {
			return name.hashCode();
		}
		
		@Override
		public boolean equals(Object obj) {
			Variant that = (Variant) obj;
			
			return ((this.name.equals(that.name)) &&
				(this.chr.equals(that.chr)) &&
				(this.pos == that.pos) &&
				(this.length == that.length));
		}
		
		public String toString() {
			return name + ":" + chr + ":" + pos + ":" + length;
		}
	}
	
	public static void main(String[] args) throws Exception {
		DbSnpSearcher s = new DbSnpSearcher();
//		s.run("/home/lmose/dev/ayc/germline_analysis/dbsnp_indels.txt.gz", "/home/lmose/dev/ayc/germline_analysis/abra_only.txt");
//		s.run("/home/lmose/dev/ayc/germline_analysis/dbsnp_indels.txt.gz", "/home/lmose/dev/ayc/germline_analysis/all_abra_only.txt");
//		s.run("/home/lmose/dev/ayc/germline_analysis/dbsnp_indels.txt.gz", "/home/lmose/dev/ayc/germline_analysis/abra_only_qual.txt");
//		s.run("/home/lmose/dev/ayc/germline_analysis/dbsnp_indels.txt.gz", "/home/lmose/dev/ayc/germline_analysis/abra_only_high_freq.txt");
//		s.run("/home/lmose/dev/ayc/germline_analysis/dbsnp_indels.txt.gz", "/home/lmose/Documents/abra/abra46_results/abra_germline_indels2.txt");
//		s.run("/home/lmose/dev/ayc/germline_analysis/dbsnp_indels.txt.gz", "/home/lmose/Documents/abra/abra46_results/gl3.tsv");
//		s.run("/home/lmose/dev/ayc/germline_analysis/dbsnp_indels.txt.gz", "/home/lmose/dev/ayc/germline_analysis/round2/abra_only.txt");
		s.run("/home/lmose/dev/ayc/germline_analysis/dbsnp_indels.txt.gz", "/home/lmose/dev/abra/calls2/germline/data/all_germline_sorted.txt");
	}
}

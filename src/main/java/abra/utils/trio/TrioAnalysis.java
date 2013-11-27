package abra.utils.trio;

import java.io.IOException;

public class TrioAnalysis {
	
	public void run(String chromosomes, String father, String mother, String child, TrioVcfReader.Caller caller) throws IOException {
		TrioVcfReader rdr = new TrioVcfReader(chromosomes, father, mother, child, caller);
    	
    	int l = 0;
    	for (TrioGenotype gt : rdr) {
    		if (gt.hasVariant()) {
    			
    			System.out.println(gt.summary());
    		}
    		
//    		if (l++ > 500) break;
    	}
	}
	
	public static void usage() {
		System.out.println("TrioAnalysis <caller: [fb|gatk]> <chromosome_file> <child.vcf> <mother.vcf> <father.vcf>");
		System.exit(-1);
	}

    public static void main(String[] args) throws Exception {
    	
    	if (args.length != 5) {
    		usage();
    	}
    	String caller = args[0];
    	String chromosomes = args[1];
    	String child = args[2];
    	String mother = args[3];
    	String father = args[4];
    	
    	TrioVcfReader.Caller callerType = null;
    	if (caller.equals("fb")) {
    		callerType = TrioVcfReader.Caller.FREEBAYES;
    	} else if (caller.equals("gatk")) {
    		callerType = TrioVcfReader.Caller.GATK;
    	} else {
    		usage();
    	}
    	
    	TrioAnalysis ta = new TrioAnalysis();
    	ta.run(chromosomes, father, mother, child, callerType);
    	
//    	String chromosomes = "/home/lmose/dev/abra/trio/chromosomes.txt";
    	
    	// Isaac

//    	String father = "/home/lmose/dev/abra/1000g/isaac/sift/father.sift.vcf";
//    	String mother = "/home/lmose/dev/abra/1000g/isaac/sift/mother.sift.vcf";
//    	String child = "/home/lmose/dev/abra/1000g/isaac/sift/child.sift.vcf";
    	
//    	String father = "/home/lmose/dev/abra/1000g/isaac/sift/father.abra.sift.vcf";
//    	String mother = "/home/lmose/dev/abra/1000g/isaac/sift/mother.abra.sift.vcf";
//    	String child = "/home/lmose/dev/abra/1000g/isaac/sift/child.abra.sift.vcf";


//    	String father = "/home/lmose/dev/abra/1000g/isaac/father.abra.pass.vcf";
//    	String mother = "/home/lmose/dev/abra/1000g/isaac/mother.abra.pass.vcf";
//    	String child = "/home/lmose/dev/abra/1000g/isaac/child.abra.pass.vcf";

    	
//    	String father = "/home/lmose/dev/abra/1000g/isaac/father.vcf";
//    	String mother = "/home/lmose/dev/abra/1000g/isaac/mother.vcf";
//    	String child = "/home/lmose/dev/abra/1000g/isaac/child.vcf";
    	
//    	String father = "/home/lmose/dev/abra/1000g/isaac/father.abra.vcf";
//    	String mother = "/home/lmose/dev/abra/1000g/isaac/mother.abra.vcf";
//    	String child = "/home/lmose/dev/abra/1000g/isaac/child.abra.vcf";


    	
//    	String father = "/home/lmose/dev/abra/trio_wxs/isaac/no_sift/father.wxs.vcf";
//    	String mother = "/home/lmose/dev/abra/trio_wxs/isaac/no_sift/mother.wxs.vcf";
//    	String child = "/home/lmose/dev/abra/trio_wxs/isaac/no_sift/child.wxs.vcf";
    	
//    	String father = "/home/lmose/dev/abra/trio_wxs/isaac/no_sift/father.abra.wxs.vcf";
//    	String mother = "/home/lmose/dev/abra/trio_wxs/isaac/no_sift/mother.abra.wxs.vcf";
//    	String child = "/home/lmose/dev/abra/trio_wxs/isaac/no_sift/child.abra.wxs.vcf";

//    	String father = "/home/lmose/dev/abra/trio_wxs/isaac/no_sift/father.abra.trio.wxs.vcf";
//    	String mother = "/home/lmose/dev/abra/trio_wxs/isaac/no_sift/mother.abra.trio.wxs.vcf";
//    	String child = "/home/lmose/dev/abra/trio_wxs/isaac/no_sift/child.abra.trio.wxs.vcf";

//    	String father = "/home/lmose/dev/abra/trio_wxs/isaac/father.wxs.sift.vcf";
//    	String mother = "/home/lmose/dev/abra/trio_wxs/isaac/mother.wxs.sift.vcf";
//    	String child = "/home/lmose/dev/abra/trio_wxs/isaac/child.wxs.sift.vcf";

//    	String father = "/home/lmose/dev/abra/trio_wxs/isaac/father.abra.wxs.sift.vcf";
//    	String mother = "/home/lmose/dev/abra/trio_wxs/isaac/mother.abra.wxs.sift.vcf";
//    	String child = "/home/lmose/dev/abra/trio_wxs/isaac/child.abra.wxs.sift.vcf";
    	
//    	String father = "/home/lmose/dev/abra/trio_wxs/isaac/father.abra.pass.vcf";
//    	String mother = "/home/lmose/dev/abra/trio_wxs/isaac/mother.abra.pass.vcf";
//    	String child = "/home/lmose/dev/abra/trio_wxs/isaac/child.abra.pass.vcf";
    	
//    	String father = "/home/lmose/dev/abra/trio_wxs/isaac/father.pass.vcf";
//    	String mother = "/home/lmose/dev/abra/trio_wxs/isaac/mother.pass.vcf";
//    	String child = "/home/lmose/dev/abra/trio_wxs/isaac/child.pass.vcf";


    	
//    	String father = "/home/lmose/dev/abra/trio_wxs/isaac/father.abra.trio.wxs.sift.vcf";
//    	String mother = "/home/lmose/dev/abra/trio_wxs/isaac/mother.abra.trio.wxs.sift.vcf";
//    	String child = "/home/lmose/dev/abra/trio_wxs/isaac/child.abra.trio.wxs.sift.vcf";

    	
    	// UnifiedGenotyper

//    	String father = "/home/lmose/dev/abra/hapmap/round3.4/ug/father.abra3.trio.sort.ug.coding.vcf";
//    	String mother = "/home/lmose/dev/abra/hapmap/round3.4/ug/mother.abra3.trio.sort.ug.coding.vcf";
//    	String child = "/home/lmose/dev/abra/hapmap/round3.4/ug/child.abra3.trio.sort.ug.coding.vcf";
    	
//    	String father = "/home/lmose/dev/abra/hapmap/round3.4/ug/father.wxs.abra3.sort.ug.coding.vcf";
//    	String mother = "/home/lmose/dev/abra/hapmap/round3.4/ug/mother.wxs.abra3.sort.ug.coding.vcf";
//    	String child = "/home/lmose/dev/abra/hapmap/round3.4/ug/child.wxs.abra3.sort.ug.coding.vcf";

    	
//    	String father = "/home/lmose/dev/abra/hapmap/round3.4/ug/father.wxs.ug.coding.vcf";
//    	String mother = "/home/lmose/dev/abra/hapmap/round3.4/ug/mother.wxs.ug.coding.vcf";
//    	String child = "/home/lmose/dev/abra/hapmap/round3.4/ug/child.wxs.ug.coding.vcf";

    	
//    	String father = "/home/lmose/dev/abra/hapmap/round3.2/ug/father.wxs.abra3.sort.ug.coding.vcf";
//    	String mother = "/home/lmose/dev/abra/hapmap/round3.2/ug/mother.wxs.abra3.sort.ug.coding.vcf";
//    	String child = "/home/lmose/dev/abra/hapmap/round3.2/ug/child.wxs.abra3.sort.ug.coding.vcf";
    	
//    	String father = "/home/lmose/dev/abra/hapmap/round3.2/ug/father.wxs.ug.coding.vcf";
//    	String mother = "/home/lmose/dev/abra/hapmap/round3.2/ug/mother.wxs.ug.coding.vcf";
//    	String child = "/home/lmose/dev/abra/hapmap/round3.2/ug/child.wxs.ug.coding.vcf";

    	
//    	String father = "/home/lmose/dev/abra/hapmap/round3/ug/father.abra2.trio.coding.vcf";
//    	String mother = "/home/lmose/dev/abra/hapmap/round3/ug/mother.abra2.trio.coding.vcf";
//    	String child = "/home/lmose/dev/abra/hapmap/round3/ug/child.abra2.trio.coding.vcf";

    	
//    	String father = "/home/lmose/dev/abra/hapmap/round3/ug/father.abra2.coding.vcf";
//    	String mother = "/home/lmose/dev/abra/hapmap/round3/ug/mother.abra2.coding.vcf";
//    	String child = "/home/lmose/dev/abra/hapmap/round3/ug/child.abra2.coding.vcf";

    	
//    	String father = "/home/lmose/dev/abra/hapmap/round3/ug/father.coding.vcf";
//    	String mother = "/home/lmose/dev/abra/hapmap/round3/ug/mother.coding.vcf";
//    	String child = "/home/lmose/dev/abra/hapmap/round3/ug/child.coding.vcf";

    	
//    	String father = "/home/lmose/dev/abra/hapmap/round3/ug/father.abra2.coding.vcf";
//    	String mother = "/home/lmose/dev/abra/hapmap/round3/ug/mother.abra2.coding.vcf";
//    	String child = "/home/lmose/dev/abra/hapmap/round3/ug/child.abra2.coding.vcf";

    	
//    	String father = "/home/lmose/dev/abra/trio_wxs/gatk_ug/father.wxs.ug.sift.vcf";
//    	String mother = "/home/lmose/dev/abra/trio_wxs/gatk_ug/mother.wxs.ug.sift.vcf";
//    	String child = "/home/lmose/dev/abra/trio_wxs/gatk_ug/child.wxs.ug.sift.vcf";
    	
//    	String father = "/home/lmose/dev/abra/trio_wxs/gatk_ug/father.abra.wxs.ug.sift.vcf";
//    	String mother = "/home/lmose/dev/abra/trio_wxs/gatk_ug/mother.abra.wxs.ug.sift.vcf";
//    	String child = "/home/lmose/dev/abra/trio_wxs/gatk_ug/child.abra.wxs.ug.sift.vcf";

    	// HaplotypeCaller

    	
//    	String father = "/home/lmose/dev/abra/hapmap/round5/hc/father.wxs.hc.coding.vcf";
//    	String mother = "/home/lmose/dev/abra/hapmap/round5/hc/mother.wxs.hc.coding.vcf";
//    	String child = "/home/lmose/dev/abra/hapmap/round5/hc/child.wxs.hc.coding.vcf";

    	//    	String father = "/home/lmose/dev/abra/hapmap/round3.4/hc/father.wxs.hc.coding.vcf";
//    	String mother = "/home/lmose/dev/abra/hapmap/round3.4/hc/mother.wxs.hc.coding.vcf";
//    	String child = "/home/lmose/dev/abra/hapmap/round3.4/hc/child.wxs.hc.coding.vcf";

    	
//    	String father = "/home/lmose/dev/abra/hapmap/round3.3/hc/father.wxs.hc.coding.vcf";
//    	String mother = "/home/lmose/dev/abra/hapmap/round3.3/hc/mother.wxs.hc.coding.vcf";
//    	String child = "/home/lmose/dev/abra/hapmap/round3.3/hc/child.wxs.hc.coding.vcf";
    	
//    	String father = "/home/lmose/dev/abra/hapmap/round3/hc/father.coding.vcf";
//    	String mother = "/home/lmose/dev/abra/hapmap/round3/hc/mother.coding.vcf";
//    	String child = "/home/lmose/dev/abra/hapmap/round3/hc/child.coding.vcf";

    	
//    	String father = "/home/lmose/dev/abra/hapmap/round2/hc/father.wxs.hc.filt.ann.coding.vcf";
//    	String mother = "/home/lmose/dev/abra/hapmap/round2/hc/mother.wxs.hc.filt.ann.coding.vcf";
//    	String child = "/home/lmose/dev/abra/hapmap/round2/hc/child.wxs.hc.filt.ann.coding.vcf";

    	
//    	String father = "/home/lmose/dev/abra/hapmap/round1/hc/father.wxs.hc.filt.vcf";
//    	String mother = "/home/lmose/dev/abra/hapmap/round1/hc/mother.wxs.hc.filt.vcf";
//    	String child = "/home/lmose/dev/abra/hapmap/round1/hc/child.wxs.hc.filt.vcf";
    	
//    	String father = "/home/lmose/dev/abra/1000g/round5/hc/father.hc.coding.vcf";
//    	String mother = "/home/lmose/dev/abra/1000g/round5/hc/mother.hc.coding.vcf";
//    	String child = "/home/lmose/dev/abra/1000g/round5/hc/child.hc.coding.vcf";
    	
//    	String father = "/home/lmose/dev/abra/1000g/round4/hc/father.hc.ann.vcf";
//    	String mother = "/home/lmose/dev/abra/1000g/round4/hc/mother.hc.ann.vcf";
//    	String child = "/home/lmose/dev/abra/1000g/round4/hc/child.hc.ann.vcf";

    	
//    	String father = "/home/lmose/dev/abra/1000g/round3/hc/father.hc.nofp.vcf";
//    	String mother = "/home/lmose/dev/abra/1000g/round3/hc/mother.hc.nofp.vcf";
//    	String child = "/home/lmose/dev/abra/1000g/round3/hc/child.hc.nofp.vcf";
    	
//    	String father = "/home/lmose/dev/abra/1000g/gatk.hc/father.hc.filt.vcf";
//    	String mother = "/home/lmose/dev/abra/1000g/gatk.hc/mother.hc.filt.vcf";
//    	String child = "/home/lmose/dev/abra/1000g/gatk.hc/child.hc.filt.vcf";

    	
//    	String father = "/home/lmose/dev/abra/trio_wxs/gatk_hc/father.wxs.hc.sift.vcf";
//    	String mother = "/home/lmose/dev/abra/trio_wxs/gatk_hc/mother.wxs.hc.sift.vcf";
//    	String child = "/home/lmose/dev/abra/trio_wxs/gatk_hc/child.wxs.hc.sift.vcf";

    	
    	// HaplotypeCaller 2000 ...
//    	String father = "/home/lmose/dev/abra/trio_wxs/gatk_hc2000/father.wxs.hc.sift.vcf";
//    	String mother = "/home/lmose/dev/abra/trio_wxs/gatk_hc2000/mother.wxs.hc.sift.vcf";
//    	String child = "/home/lmose/dev/abra/trio_wxs/gatk_hc2000/child.wxs.hc.sift.vcf";
    	
//    	String father = "/home/lmose/dev/abra/trio_wxs/gatk_hc2000/father.abra.trio.wxs.hc.sift.vcf";
//    	String mother = "/home/lmose/dev/abra/trio_wxs/gatk_hc2000/mother.abra.trio.wxs.hc.sift.vcf";
//    	String child = "/home/lmose/dev/abra/trio_wxs/gatk_hc2000/child.abra.trio.wxs.hc.sift.vcf";

    	// Freebayes F .001
//    	String father = "/home/lmose/dev/abra/trio_wxs/freebayes.f0/father.abra.wxs.prim.sift.vcf";
//    	String mother = "/home/lmose/dev/abra/trio_wxs/freebayes.f0/mother.abra.wxs.prim.sift.vcf";
//    	String child = "/home/lmose/dev/abra/trio_wxs/freebayes.f0/child.abra.wxs.prim.sift.vcf";
    	
//    	String father = "/home/lmose/dev/abra/trio_wxs/freebayes.f0/father.wxs.prim.sift.vcf";
//    	String mother = "/home/lmose/dev/abra/trio_wxs/freebayes.f0/mother.wxs.prim.sift.vcf";
//    	String child = "/home/lmose/dev/abra/trio_wxs/freebayes.f0/child.wxs.prim.sift.vcf";


    	
    	// Freebayes
    	
//    	String father = "/home/lmose/dev/abra/hapmap/round5/fb/father.53.73.abra.sort.F05.coding.vcf";
//    	String mother = "/home/lmose/dev/abra/hapmap/round5/fb/mother.53.73.abra.sort.F05.coding.vcf";
//    	String child = "/home/lmose/dev/abra/hapmap/round5/fb/child.53.73.abra.sort.F05.coding.vcf";
    	
//    	String father = "/home/lmose/dev/abra/hapmap/round5/fb/father.43_83.abra.sort.F05.coding.vcf";
//    	String mother = "/home/lmose/dev/abra/hapmap/round5/fb/mother.43_83.abra.sort.F05.coding.vcf";
//    	String child = "/home/lmose/dev/abra/hapmap/round5/fb/child.43_83.abra.sort.F05.coding.vcf";

//    	String father = "/home/lmose/dev/abra/hapmap/round5/hc/father.wxs.hc.coding.a5.vcf";
//    	String mother = "/home/lmose/dev/abra/hapmap/round5/hc/mother.wxs.hc.coding.a5.vcf";
//    	String child = "/home/lmose/dev/abra/hapmap/round5/hc/child.wxs.hc.coding.a5.vcf";
    	
    	
    	    	
//    	String father = "/home/lmose/dev/abra/hapmap/round5/fb/father.43_83.abra.sort.F05.coding.a5.vcf";
//    	String mother = "/home/lmose/dev/abra/hapmap/round5/fb/mother.43_83.abra.sort.F05.coding.a5.vcf";
//    	String child = "/home/lmose/dev/abra/hapmap/round5/fb/child.43_83.abra.sort.F05.coding.a5.vcf";

    	
//    	String father = "/home/lmose/dev/abra/hapmap/round4/fb/father.abra.sort.coding.vcf";
//    	String mother = "/home/lmose/dev/abra/hapmap/round4/fb/mother.abra.sort.coding.vcf";
//    	String child = "/home/lmose/dev/abra/hapmap/round4/fb/child.abra.sort.coding.vcf";

    	
//    	String father = "/home/lmose/dev/abra/hapmap/round3.4/fb/father.wxs.abra3.sort.coding.vcf";
//    	String mother = "/home/lmose/dev/abra/hapmap/round3.4/fb/mother.wxs.abra3.sort.coding.vcf";
//    	String child = "/home/lmose/dev/abra/hapmap/round3.4/fb/child.wxs.abra3.sort.coding.vcf";

    	
//    	String father = "/home/lmose/dev/abra/hapmap/round3.4/fb/father.wxs.coding.vcf";
//    	String mother = "/home/lmose/dev/abra/hapmap/round3.4/fb/mother.wxs.coding.vcf";
//    	String child = "/home/lmose/dev/abra/hapmap/round3.4/fb/child.wxs.coding.vcf";


//    	String father = "/home/lmose/dev/abra/hapmap/round3.3/fb001/father.wxs.abra3.sort.coding.vcf";
//    	String mother = "/home/lmose/dev/abra/hapmap/round3.3/fb001/mother.wxs.abra3.sort.coding.vcf";
//    	String child = "/home/lmose/dev/abra/hapmap/round3.3/fb001/child.wxs.abra3.sort.coding.vcf";

    	
//    	String father = "/home/lmose/dev/abra/hapmap/round3.3/fb001/father.wxs.coding.vcf";
//    	String mother = "/home/lmose/dev/abra/hapmap/round3.3/fb001/mother.wxs.coding.vcf";
//    	String child = "/home/lmose/dev/abra/hapmap/round3.3/fb001/child.wxs.coding.vcf";

    	
//    	String father = "/home/lmose/dev/abra/hapmap/round3.3/fb/father.wxs.coding.vcf";
//    	String mother = "/home/lmose/dev/abra/hapmap/round3.3/fb/mother.wxs.coding.vcf";
//    	String child = "/home/lmose/dev/abra/hapmap/round3.3/fb/child.wxs.coding.vcf";

    	
//    	String father = "/home/lmose/dev/abra/hapmap/round3.3/fb/father.wxs.abra3.sort.coding.vcf";
//    	String mother = "/home/lmose/dev/abra/hapmap/round3.3/fb/mother.wxs.abra3.sort.coding.vcf";
//    	String child = "/home/lmose/dev/abra/hapmap/round3.3/fb/child.wxs.abra3.sort.coding.vcf";


//    	String father = "/home/lmose/dev/abra/hapmap/round3.2/fb/father.abra3.trio.sort.coding.vcf";
//    	String mother = "/home/lmose/dev/abra/hapmap/round3.2/fb/mother.abra3.trio.sort.coding.vcf";
//    	String child = "/home/lmose/dev/abra/hapmap/round3.2/fb/child.abra3.trio.sort.coding.vcf";

    	
//    	String father = "/home/lmose/dev/abra/hapmap/round3.2/fb/father.wxs.coding.vcf";
//    	String mother = "/home/lmose/dev/abra/hapmap/round3.2/fb/mother.wxs.coding.vcf";
//    	String child = "/home/lmose/dev/abra/hapmap/round3.2/fb/child.wxs.coding.vcf";
    	
//    	String father = "/home/lmose/dev/abra/hapmap/round3.2/fb/father.wxs.abra3.sort.coding.vcf";
//    	String mother = "/home/lmose/dev/abra/hapmap/round3.2/fb/mother.wxs.abra3.sort.coding.vcf";
//    	String child = "/home/lmose/dev/abra/hapmap/round3.2/fb/child.wxs.abra3.sort.coding.vcf";

    	
//    	String father = "/home/lmose/dev/abra/hapmap/round3/fb/father.coding.vcf";
//    	String mother = "/home/lmose/dev/abra/hapmap/round3/fb/mother.coding.vcf";
//    	String child = "/home/lmose/dev/abra/hapmap/round3/fb/child.coding.vcf";

    	
//    	String father = "/home/lmose/dev/abra/hapmap/round3/fb/father.abra2.coding.vcf";
//    	String mother = "/home/lmose/dev/abra/hapmap/round3/fb/mother.abra2.coding.vcf";
//    	String child = "/home/lmose/dev/abra/hapmap/round3/fb/child.abra2.coding.vcf";

    	
//    	String father = "/home/lmose/dev/abra/hapmap/round2/fb/father.wxs.filt.prim.ann.coding.vcf";
//    	String mother = "/home/lmose/dev/abra/hapmap/round2/fb/mother.wxs.filt.prim.ann.coding.vcf";
//    	String child = "/home/lmose/dev/abra/hapmap/round2/fb/child.wxs.filt.prim.ann.coding.vcf";

    	
//    	String father = "/home/lmose/dev/abra/hapmap/round2/fb/father.abra.wxs.filt.prim.ann.coding.vcf";
//    	String mother = "/home/lmose/dev/abra/hapmap/round2/fb/mother.abra.wxs.filt.prim.ann.coding.vcf";
//    	String child = "/home/lmose/dev/abra/hapmap/round2/fb/child.abra.wxs.filt.prim.ann.coding.vcf";

    	
//    	String father = "/home/lmose/dev/abra/hapmap/round1/fb/father.abra.wxs.filt.prim.vcf";
//    	String mother = "/home/lmose/dev/abra/hapmap/round1/fb/mother.abra.wxs.filt.prim.vcf";
//    	String child = "/home/lmose/dev/abra/hapmap/round1/fb/child.abra.wxs.filt.prim.vcf";

    	
//    	String father = "/home/lmose/dev/abra/hapmap/round1/fb/father.wxs.filt.prim.vcf";
//    	String mother = "/home/lmose/dev/abra/hapmap/round1/fb/mother.wxs.filt.prim.vcf";
//    	String child = "/home/lmose/dev/abra/hapmap/round1/fb/child.wxs.filt.prim.vcf";
    	
//    	String father = "/home/lmose/dev/abra/1000g/round5/freebayes/father.abra.fb.coding.vcf";
//    	String mother = "/home/lmose/dev/abra/1000g/round5/freebayes/mother.abra.fb.coding.vcf";
//    	String child = "/home/lmose/dev/abra/1000g/round5/freebayes/child.abra.fb.coding.vcf";
    	
//    	String father = "/home/lmose/dev/abra/1000g/round5/freebayes/father.fb.coding.vcf";
//    	String mother = "/home/lmose/dev/abra/1000g/round5/freebayes/mother.fb.coding.vcf";
//    	String child = "/home/lmose/dev/abra/1000g/round5/freebayes/child.fb.coding.vcf";
    	
//    	String father = "/home/lmose/dev/abra/1000g/round4/fb/freebayes_father.abra.sort.prim.ann.vcf";
//    	String mother = "/home/lmose/dev/abra/1000g/round4/fb/freebayes_mother.abra.sort.prim.ann.vcf";
//    	String child = "/home/lmose/dev/abra/1000g/round4/fb/freebayes_child.abra.sort.prim.ann.vcf";

    	
//    	String father = "/home/lmose/dev/abra/1000g/round4/fb/freebayes_father.prim.ann.vcf";
//    	String mother = "/home/lmose/dev/abra/1000g/round4/fb/freebayes_mother.prim.ann.vcf";
//    	String child = "/home/lmose/dev/abra/1000g/round4/fb/freebayes_child.prim.ann.vcf";

    	
//    	String father = "/home/lmose/dev/abra/1000g/round3/fb/father.fb.abra.prim.nofp.vcf";
//    	String mother = "/home/lmose/dev/abra/1000g/round3/fb/mother.fb.abra.prim.nofp.vcf";
//    	String child = "/home/lmose/dev/abra/1000g/round3/fb/child.fb.abra.prim.nofp.vcf";
    	
//    	String father = "/home/lmose/dev/abra/1000g/freebayes2/freebayes_father.prim.vcf";
//    	String mother = "/home/lmose/dev/abra/1000g/freebayes2/freebayes_mother.prim.vcf";
//    	String child = "/home/lmose/dev/abra/1000g/freebayes2/freebayes_child.prim.vcf";
    	
    	
//    	String father = "/home/lmose/dev/abra/1000g/freebayes2/freebayes_father.abra.sort.prim.vcf";
//    	String mother = "/home/lmose/dev/abra/1000g/freebayes2/freebayes_mother.abra.sort.prim.vcf";
//    	String child = "/home/lmose/dev/abra/1000g/freebayes2/freebayes_child.abra.sort.prim.vcf";
    	
//    	String father = "/home/lmose/dev/abra/1000g/freebayes/freebayes_father.filt.prim.vcf";
//    	String mother = "/home/lmose/dev/abra/1000g/freebayes/freebayes_mother.filt.prim.vcf";
//    	String child = "/home/lmose/dev/abra/1000g/freebayes/freebayes_child.filt.prim.vcf";
    	
//    	String father = "/home/lmose/dev/abra/1000g/freebayes/freebayes_father.abra.sort.filt.prim.vcf";
//    	String mother = "/home/lmose/dev/abra/1000g/freebayes/freebayes_mother.abra.sort.filt.prim.vcf";
//    	String child = "/home/lmose/dev/abra/1000g/freebayes/freebayes_child.abra.sort.filt.prim.vcf";

    	
    	
//    	String father = "/home/lmose/dev/abra/trio_wxs/freebayes/father.abra.wxs.prim.sift.vcf";
//    	String mother = "/home/lmose/dev/abra/trio_wxs/freebayes/mother.abra.wxs.prim.sift.vcf";
//    	String child = "/home/lmose/dev/abra/trio_wxs/freebayes/child.abra.wxs.prim.sift.vcf";

//    	String father = "/home/lmose/dev/abra/trio_wxs/freebayes/father.wxs.prim.sift.vcf";
//    	String mother = "/home/lmose/dev/abra/trio_wxs/freebayes/mother.wxs.prim.sift.vcf";
//    	String child = "/home/lmose/dev/abra/trio_wxs/freebayes/child.wxs.prim.sift.vcf";

    	
    	// Freebayes...
//    	String father = "/home/lmose/dev/abra/trio_wxs/prim/father.abra3.prim.vcf";
//    	String mother = "/home/lmose/dev/abra/trio_wxs/prim/mother.abra3.prim.vcf";
//    	String child = "/home/lmose/dev/abra/trio_wxs/prim/child.abra3.prim.vcf";
    	
//    	String father = "/home/lmose/dev/abra/trio_wxs/prim/father.abra.prim.vcf";
//    	String mother = "/home/lmose/dev/abra/trio_wxs/prim/mother.abra.prim.vcf";
//    	String child = "/home/lmose/dev/abra/trio_wxs/prim/child.abra.prim.vcf";
    	
//    	String father = "/home/lmose/dev/abra/trio_wxs/prim/father.prim.vcf";
//    	String mother = "/home/lmose/dev/abra/trio_wxs/prim/mother.prim.vcf";
//    	String child = "/home/lmose/dev/abra/trio_wxs/prim/child.prim.vcf";
    	
    	//////////////////////////////////////////////////////////////
    	
//    	String father = "/home/lmose/dev/abra/trio_wxs/father.vcf";
//    	String mother = "/home/lmose/dev/abra/trio_wxs/mother.vcf";
//    	String child = "/home/lmose/dev/abra/trio_wxs/child.vcf";
    	
//    	String father = "/home/lmose/dev/abra/trio_wxs/father.abra.sort.vcf";
//    	String mother = "/home/lmose/dev/abra/trio_wxs/mother.abra.sort.vcf";
//    	String child = "/home/lmose/dev/abra/trio_wxs/child.abra.sort.vcf";

//    	String father = "/home/lmose/dev/abra/trio_wxs/trio/father.abra.trio2.sort.vcf";
//    	String mother = "/home/lmose/dev/abra/trio_wxs/trio/mother.abra.trio2.sort.vcf";
//    	String child = "/home/lmose/dev/abra/trio_wxs/trio/child.abra.trio2.sort.vcf";


//    	String father = "/home/lmose/dev/abra/triofb/freebayes.father_orig.vcf";
//    	String mother = "/home/lmose/dev/abra/triofb/freebayes.mother_orig.vcf";
//    	String child = "/home/lmose/dev/abra/triofb/freebayes.child_orig.vcf";

    	
//    	String father = "/home/lmose/dev/abra/trio/father.vcf";
//    	String mother = "/home/lmose/dev/abra/trio/mother.vcf";
//    	String child = "/home/lmose/dev/abra/trio/child.vcf";

////    	String chromosomes = "/home/lmose/dev/abra/trio_wxs/chromosomes.txt";
//    	String chromosomes = "/home/lmose/dev/abra/1000g/chromosomes.txt";
    	
//    	String father = "/home/lmose/dev/abra/trio/pre/father.vcf";
//    	String mother = "/home/lmose/dev/abra/trio/pre/mother.vcf";
//    	String child = "/home/lmose/dev/abra/trio/pre/child.vcf";
    	
//    	String father = "/home/lmose/dev/abra/trio/multi/father.abra.sort.genome.vcf";
//    	String mother = "/home/lmose/dev/abra/trio/multi/mother.abra.sort.genome.vcf";
//    	String child = "/home/lmose/dev/abra/trio/multi/child.abra.sort.genome.vcf";

//    	String father = "/home/lmose/dev/abra/trio/trio41/father.genome.vcf";
//    	String mother = "/home/lmose/dev/abra/trio/trio41/mother.genome.vcf";
//    	String child = "/home/lmose/dev/abra/trio/trio41/child.genome.vcf";

//    	String father = "/home/lmose/dev/abra/trio/sep_wxs/father.vcf";
//    	String mother = "/home/lmose/dev/abra/trio/sep_wxs/mother.vcf";
//    	String child = "/home/lmose/dev/abra/trio/sep_wxs/child.vcf";

//    	String father = "/home/lmose/dev/abra/trio/pre/sep_wxs/father.vcf";
//    	String mother = "/home/lmose/dev/abra/trio/pre/sep_wxs/mother.vcf";
//    	String child = "/home/lmose/dev/abra/trio/pre/sep_wxs/child.vcf";
    	
//    	String father = "/home/lmose/dev/abra/trio/trio3/father.vcf";
//    	String mother = "/home/lmose/dev/abra/trio/trio3/mother.vcf";
//    	String child = "/home/lmose/dev/abra/trio/trio3/child.vcf";
    	
//    	String father = "/home/lmose/dev/abra/trio/trio4/father.genome.vcf";
//    	String mother = "/home/lmose/dev/abra/trio/trio4/mother.genome.vcf";
//    	String child = "/home/lmose/dev/abra/trio/trio4/child.genome.vcf";


    	
//    	TrioAnalysis ta = new TrioAnalysis();
//    	ta.run(chromosomes, father, mother, child, TrioVcfReader.Caller.FREEBAYES);
    }
}


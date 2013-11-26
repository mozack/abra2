/* Copyright 2013 University of North Carolina at Chapel Hill.  All rights reserved. */
package abra.utils;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.HashMap;
import java.util.Map;

public class ProcessFreebayesVcf {
	
	public void processIndels(String participantId, String source, String vcf) throws Exception {
		BufferedReader reader = new BufferedReader(new FileReader(vcf));
		
		String line = reader.readLine();
		while (line != null) {
			
			if (!line.startsWith("#") &&
				((line.contains("HIGH") || line.contains("MODERATE")))) {
				
				String[] fields = line.split("\t");
				String chr = fields[0];
				Integer pos = Integer.parseInt(fields[1]);
				String ref = fields[3];
				String alt = fields[4];
				if (!alt.contains(",")) {
					int qual = (int) Float.parseFloat(fields[5]);
//					String filter = fields[6];
					String infoStr = fields[7];
					String normalStr = fields[9];
					
					String indelType = "UNK";
					if (ref.length() > 1) {
						indelType = quote("DEL");
					} else if (alt.length() > 1) {
						indelType = quote("INS");
					}
					
					String info = parseInfo(infoStr);
					
					String countsAndGt = parseFormat(normalStr);
					
					String varId = chr + ":" + pos + ":" + ref + ":" + alt + ":" + participantId;
					
					StringBuffer str = new StringBuffer();
					str.append(quote(varId));
					str.append('\t');
					str.append(quote(participantId));
					str.append('\t');				
					str.append(quote(chr));
					str.append('\t');
					str.append(pos);
					str.append('\t');
					str.append(quote(ref));
					str.append('\t');
					str.append(quote(alt));
					str.append('\t');
					str.append(indelType);
					str.append('\t');
					str.append(info);
					str.append('\t');
					str.append(qual);
					str.append('\t');
					str.append(countsAndGt);
					
					System.out.println(sql(str.toString()));
				}
			}
			
			line = reader.readLine();
		}
		
		reader.close();
	}
	
	private String sql(String rec) {
		rec = rec.replaceAll("\t", ",");
		
//		String sql = "INSERT INTO somatic_indel (" +
//				"var_id,participant_id,chromosome,pos,gref,alt,source,indel_type,dp_filter,repeat_filter,ihpol_filter,bcnoise_filter,qsi_ref_filter,qsi,tqsi,nt,qsi_nt,tqsi_nt,sgt,ru,rc,ic,ihp,effect,impact,genes,gene,normal_dp1,normal_dp2,normal_tar1,normal_tar2,normal_tir1,normal_tir2,tumor_dp1,tumor_dp2,tumor_tar1,tumor_tar2,tumor_tir1,tumor_tir2" +
//				") VALUES ("+
//				rec + 
//				");";
		
		String sql = "INSERT INTO foo_brca_germline_indel (" +
				"var_id,participant_id,chromosome,pos,gref,alt,indel_type,effect,impact,genes,gene,qual,depth,ref_cnt,alt_cnt,gt" +
				") VALUES ("+
				rec + 
				");";

		
		return sql;
	}
	
	private String parseFormat(String format) {
		StringBuffer str = new StringBuffer();
		//DP:DP2:TAR:TIR:TOR:DP50:FDP50:SUBDP50
		//394:394:205,205:171,172:21,21:384.36:0.51:0.00
		
		// GT:DP:RO:QR:AO:QA:GL	
		// 0/1:23:13:443:10:338:-30.758,-0.865239,-40.2108
		
		// primitives
		// GT
		// 1|1
		
		format = format.replace("|", "/");
		
		String gt = format.split(":")[0];
		
		return "0\t0\t0\t" + quote(gt);
		
		/*
		String[] fields = format.split(":");
		
		int depth = Integer.parseInt(fields[1]);
		int refCount = Integer.parseInt(fields[2]);
		int altCount = Integer.parseInt(fields[4]);
		
		str.append(depth);
		str.append('\t');
		str.append(refCount);
		str.append('\t');
		str.append(altCount);
		
		return str.toString();
		*/
	}
	
	private String parseFilter(String filter) {
		// DP, Repeat, iHpol, BCNoise, QSI_ref
		StringBuffer output = new StringBuffer();
		output.append(filter.contains("DP") ? 1 : 0);
		output.append('\t');
		output.append(filter.contains("Repeat") ? 1 : 0);
		output.append('\t');
		output.append(filter.contains("iHpol") ? 1 : 0);
		output.append('\t');
		output.append(filter.contains("BCNoise") ? 1 : 0);
		output.append('\t');
		output.append(filter.contains("QSI_ref") ? 1 : 0);
		
		return output.toString();
	}
	
	private Map<String, String> getInfoMap(String info) {
		
		Map<String, String> map = new HashMap<String, String>();
		
		String[] fields = info.split(";");
		
		for (int i=0; i<fields.length; i++) {
			String[] field = fields[i].split("=");
			if (field.length == 2) {
				map.put(field[0], field[1]);
			}
		}
		
		return map;
	}
	
	private String appendString(String orig, String addition) {
		if (orig == null) {
			orig = addition;
		} else {
			orig = orig + "," + addition;
		}
		
		return orig;
	}
	
	private String parseEffect(String effectStr) {
		
		if (effectStr == null) {
			return "''\t''\t''\t''";
		}
		
		String[] effects = effectStr.split(",");
		
		String varEffects = null;
		String impacts = null;
		String classes = null;
		String genes = null;
		String collapsedGene = null;
		
		for (String effect : effects) {
			String[] split1 = effect.split("\\(");
			String varEffect = split1[0];
			if ((!varEffect.equals("UPSTREAM")) && (!varEffect.equals("DOWNSTREAM"))) {
				String fieldStr = split1[1];
				String[] fields = fieldStr.split("\\|");
				String impact = fields[0];
				String functionalClass = fields[1];
				String gene = fields[4];
				
				varEffects = appendString(varEffects, varEffect);
				impacts = appendString(impacts, impact);
				classes = appendString(classes, functionalClass);
				genes = appendString(genes, gene);
				
				if (collapsedGene == null) {
					collapsedGene = gene;
				} else if (!collapsedGene.equals(gene)) {
					collapsedGene = appendString(collapsedGene, gene);
				}
			}
		}
		
//		return quote(varEffects) + "\t" + quote(impacts) + "\t" + quote(classes) + "\t" + quote(genes) + "\t" + quote(collapsedGene);
		return quote(varEffects) + "\t" + quote(impacts) + "\t" + quote(genes) + "\t" + quote(collapsedGene);
	}
	
	private String quote(String str) {
		return "'" + str + "'";
	}
	
	private String parseInfo(String info) {
		// EFF
		
		Map<String, String> map = getInfoMap(info);
		
		StringBuffer str = new StringBuffer();
		
		str.append(parseEffect(map.get("EFF")));
		
		return str.toString();
	}

	public static void main(String[] args) throws Exception {
		ProcessFreebayesVcf p = new ProcessFreebayesVcf();
//		String vcf = "/home/lmose/dev/vcf/all.somatic.indels.ann.vcf";
		String participantId = args[0];
		String vcf = args[1];
		
//		String participantId = "8631969d-4ac6-4fe1-9db3-d47db604494a";
//		String participantId = "125b0a9a-6cbe-4207-9178-e7ff03487e87";
		String source = "";
//		String vcf = "/home/lmose/dev/vcf/all.germline.indels.abra.ann.vcf";
//		String vcf = "/home/lmose/dev/vcf/freebayes/fb.test.vcf";
		
		p.processIndels(participantId, source, vcf);
	}
}

package abra.utils;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.HashMap;
import java.util.Map;

public class ProcessStrelkaVcf {
	
	public void processIndels(String participantId, String source, String vcf) throws Exception {
		BufferedReader reader = new BufferedReader(new FileReader(vcf));
		
		String line = reader.readLine();
		while (line != null) {
			
			if (!line.startsWith("#")) {
//				((line.contains("HIGH") || line.contains("MODERATE")))) {
				
				String[] fields = line.split("\t");
				String chr = fields[0];
				Integer pos = Integer.parseInt(fields[1]);
				String ref = fields[3];
				String alt = fields[4];
				String filterStr = fields[6];
				String infoStr = fields[7];
				String normalStr = fields[9];
				String tumorStr = fields[10];
				
				String indelType = "UNK";
				if (ref.length() > 1) {
					indelType = quote("DEL");
				} else if (alt.length() > 1) {
					indelType = quote("INS");
				}
				
				String filter = parseFilter(filterStr);
				String info = parseInfo(infoStr);
				
				String normalCounts = parseFormat(normalStr);
				String tumorCounts = parseFormat(tumorStr);
				
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
				str.append(quote(source));
				str.append('\t');
				str.append(indelType);
				str.append('\t');
				str.append(filter);
				str.append('\t');
				str.append(info);
				str.append('\t');
				str.append(normalCounts);
				str.append('\t');
				str.append(tumorCounts);
				
				System.out.println(sql(str.toString()));
			}
			
			line = reader.readLine();
		}
		
		reader.close();
	}
	
	private String sql(String rec) {
		rec = rec.replaceAll("\t", ",");
		
		String sql = "INSERT INTO somatic_indel (" +
				"var_id,participant_id,chromosome,pos,gref,alt,source,indel_type,dp_filter,repeat_filter,ihpol_filter,bcnoise_filter,qsi_ref_filter,qsi,tqsi,nt,qsi_nt,tqsi_nt,sgt,ru,rc,ic,ihp,effect,impact,genes,gene,normal_dp1,normal_dp2,normal_tar1,normal_tar2,normal_tir1,normal_tir2,tumor_dp1,tumor_dp2,tumor_tar1,tumor_tar2,tumor_tir1,tumor_tir2" +
				") VALUES ("+
				rec + 
				");";
		
		return sql;
	}
	
	private String parseFormat(String format) {
		StringBuffer str = new StringBuffer();
		//DP:DP2:TAR:TIR:TOR:DP50:FDP50:SUBDP50
		//394:394:205,205:171,172:21,21:384.36:0.51:0.00
		
		String[] fields = format.split(":");
		
		String dp = fields[0];
		String dp2 = fields[1];
		String[] tar = fields[2].split(",");
		String[] tir = fields[3].split(",");
		String[] tor = fields[4].split(",");
		String dp50 = fields[5];
		String fdp50 = fields[6];
		String subdp50 = fields[7];

		str.append(dp);
		str.append('\t');
		str.append(dp2);
		str.append('\t');
		str.append(tar[0]);
		str.append('\t');
		str.append(tar[1]);
		str.append('\t');
		str.append(tir[0]);
		str.append('\t');
		str.append(tir[1]);
//		str.append('\t');
//		str.append(tor[0]);
//		str.append('\t');
//		str.append(tor[1]);
//		str.append('\t');
//		str.append(dp50);
//		str.append('\t');
//		str.append(fdp50);
//		str.append('\t');
//		str.append(subdp50);
		
		return str.toString();
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
		// QSI,TQSI,NT,QSI_NT,TQSI_NT,SGT,RU,RC,IC,IHP,EFF
		
		Map<String, String> map = getInfoMap(info);
		
		StringBuffer str = new StringBuffer();
		
		str.append(map.get("QSI"));
		str.append('\t');
		str.append(map.get("TQSI"));
		str.append('\t');
		str.append(quote(map.get("NT")));
		str.append('\t');
		str.append(map.get("QSI_NT"));
		str.append('\t');
		str.append(map.get("TQSI_NT"));
		str.append('\t');
		str.append(quote(map.get("SGT")));
		str.append('\t');
		str.append(quote(map.get("RU")));
		str.append('\t');
		str.append(map.get("RC"));
		str.append('\t');
		str.append(map.get("IC"));
		str.append('\t');
		str.append(map.get("IHP"));
		str.append('\t');
		str.append(parseEffect(map.get("EFF")));
		
		return str.toString();
	}

	public static void main(String[] args) throws Exception {
		ProcessStrelkaVcf p = new ProcessStrelkaVcf();
//		String vcf = "/home/lmose/dev/vcf/all.somatic.indels.ann.vcf";
		String participantId = args[0];
		String source = args[1];
		String vcf = args[2];
		
		p.processIndels(participantId, source, vcf);
	}
}

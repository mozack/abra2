package abra;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class ChromosomeRegex {
	
	public static final String DEFAULT_SKIP_REGEX = "GL.*|hs37d5|chr.*random|chrUn.*|chrEBV|CMV|HBV|HCV.*|HIV.*|KSHV|HTLV.*|MCV|SV40|HPV.*";
	
	private static Pattern p = null;
	
	public ChromosomeRegex(String regex) {
		if (!regex.equals("none")) {
			p = Pattern.compile(regex);
		}
	}
	
	public boolean matches(String chrom) {
		boolean ret = false;
		
		if (p != null) {
			Matcher m = p.matcher(chrom);
			ret = m.matches();
		}
		
		return ret;
	}
}

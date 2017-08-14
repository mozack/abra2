package abra.cadabra;

public class RepeatUtils {

	/**
	 * Return smallest repeat unit in input string.
	 * i.e. ATATATAT returns AT.  ATCTAGC returns ATCTAGC
	 */
	public static String getRepeatUnit(String seq) {
		String unit = seq;
		boolean isRepeat = false;
		
		for (int length = 1; length <= seq.length()/2; length++) {
			
			unit = seq.substring(0, length);
			isRepeat = true;
			
			for (int i=length; i<seq.length(); i+=length) {
				if ((seq.length()-i < length) || !seq.substring(i, i+length).equals(unit)) {
					isRepeat = false;
					break;
				}
			}
			
			if (isRepeat) {
				break;
			}
		}
		
		if (!isRepeat) {
			unit = seq;
		}
		
		return unit;
	}
	
	/**
	 * Return number of times input bases occurs without gap from beginning of ref sequence.
	 * i.e. AT in ATATATGAT is a period of 3 
	 */
	public static int getRepeatPeriod(String bases, String ref) {
		int maxLen = Math.min(ref.length(), bases.length()*100);
		
		int period = 0;
		
		if (bases.length() > 0) {
			int index = 0;
			while ((index < maxLen-bases.length()+1) && (bases.equals(ref.substring(index, index+bases.length())))) {
				period += 1;
				index += bases.length();
			}
		}
		
		return period;
	}
}

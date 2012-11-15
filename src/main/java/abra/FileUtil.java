package abra;

import java.io.FileOutputStream;
import java.io.IOException;

public class FileUtil {

	public static void truncateFile(String filename) throws IOException {
		FileOutputStream fos = new FileOutputStream(filename, false);
//		FileChannel outChan = fos.getChannel();
//		outChan.truncate(0);
//		outChan.close();
		fos.close();
	}
	
	public static void main(String[] args) throws Exception {
		truncateFile("/home/lmose/dev/ayc/trunc/1.txt");
	}
}

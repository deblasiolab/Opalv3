package opal.IO;

import com.traviswheeler.libs.LogWriter;

public class FastaWriter extends AlignmentWriter {

	public FastaWriter(String[] namesA, String[] namesB, char[][] alignment, int K, int L, boolean toUpper) {
		super(namesA, namesB, alignment, K, L, toUpper);
	}

	public FastaWriter(String[] names, char[][] alignment, int K, boolean toUpper) {
		super(names, alignment, K, toUpper);
	}
	
	final protected void setDefaultOutputWidth() {
		width = 80;
	}
		
	final public void write() {
		for (int x=0; x<K+L; x++) {
			String name = x<K ? namesA[x] : namesB[x-K];
			String seq = new String(alignment[x]);
			if (toUpper)
				seq = seq.toUpperCase();
			LogWriter.stdOutLogln(">" + name);
			int start = 0;
			int end = width;
			while (end < seq.length()) {
				LogWriter.stdOutLogln( seq.substring(start, end));
				start += width;
				end += width;
			}
			LogWriter.stdOutLogln( seq.substring(start)+"\n");
		}
		LogWriter.stdOutLogln("");		
	}

}

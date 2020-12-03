package opal.IO;


import java.io.PrintStream;

import com.traviswheeler.libs.LogWriter;

public class FastaWriter extends AlignmentWriter {

	public FastaWriter(String[] namesA, String[] namesB, char[][] alignment, int K, int L, boolean toUpper, PrintStream outStream) {
		super(namesA, namesB, alignment, K, L, toUpper, outStream);
	}

	public FastaWriter(String[] names, char[][] alignment, int K, boolean toUpper, PrintStream outStream) {
		super(names, alignment, K, toUpper, outStream);
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
			//LogWriter.stdOutLogln(">" + name);
			out.println(">" + name);
			int start = 0;
			int end = width;
			while (end < seq.length()) {
				//LogWriter.stdOutLogln( seq.substring(start, end));
				out.println(seq.substring(start, end));
				start += width;
				end += width;
			}
			//LogWriter.stdOutLogln( seq.substring(start)+"\n");
			out.println(seq.substring(start)+"\n");
		}
		//LogWriter.stdOutLogln("");
		out.println("");		
	}

}

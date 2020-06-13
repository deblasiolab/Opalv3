package opal.IO;

import opal.align.Aligner;
import com.traviswheeler.libs.LogWriter;

public class ClustalWriter extends AlignmentWriter {

	int nameWidth = 16;
	

	
	public ClustalWriter(String[] namesA, String[] namesB, char[][] alignment, int K, int L, boolean toUpper) {
		super(namesA, namesB, alignment, K, L, toUpper);
		
		int maxNameLen = 0;
		for(int i=0; i<namesA.length; i++)
			if (namesA[i].length() > maxNameLen) maxNameLen = namesA[i].length() ; 
		for(int i=0; i<namesB.length; i++)
			if (namesB[i].length() > maxNameLen) maxNameLen = namesB[i].length() ;
		
		if (maxNameLen < nameWidth) {
			width += nameWidth - maxNameLen;
			nameWidth = maxNameLen;
		}
	}

	public ClustalWriter(String[] names, char[][] alignment, int K, boolean toUpper) {
		super(names, alignment, K, toUpper);
	}	

	final protected void setDefaultOutputWidth() {
		width = 55;
	}
		
	final public void write() {
		int N = alignment[0].length;
		String[] seqs = new String[K+L];
		int totalLength = N;
		int[] starts = new int[K+L];
		int[] ends = new int[K+L];
		int[] first = new int[K+L];
		int[] last = new int[K+L];
		
		
		for (int i=0; i<K+L ; i++) {
			first[i] = -1;
			for (int j=0; j<N; j++){
				if (alignment[i][j] != SequenceConverter.GAP_CHAR ) {
					last[i] = j;
					if (first[i] == -1)
						first[i] = j;
				}
			}
			starts[i] = 0;
			seqs[i] = new String(alignment[i]);
			if (toUpper)
				seqs[i] = seqs[i].toUpperCase();
		}
		
		int maxPosWidth = (totalLength+"").length();

		for (int start=0; start<totalLength; start+=width) {
			int end = start+width > totalLength ? totalLength : start+width ;
			
			for (int i=0; i<K; i++) {
				String nm = pad(namesA[i], nameWidth, true);
				String seq = seqs[i].substring(start,end);
				ends[i] = starts[i] + seq.replaceAll("-| ","").length()-1;
			
				LogWriter.stdOutLogln(nm + " " + pad(starts[i]+(start>last[i] ? 0 : 1),maxPosWidth, false) + 
						             " " + seq + " " + pad(ends[i]+1,maxPosWidth, false));
				
				starts[i] = ends[i] + 1;
			}			
			
			
			if (null != path && L > 0) {
				//print path line
				LogWriter.stdOutLog( pad("",nameWidth+2+maxPosWidth, true));
				for (int i=start; i<end; i++) {
					Aligner.Direction direc = path.get(i);
					LogWriter.stdOutLog(direc == Aligner.Direction.diag ? "|" : " ");
				}
				LogWriter.stdOutLog("\n");
			}

			for (int i=0; i<L; i++) {
				String nm = pad(namesB[i], nameWidth, true);
				String seq = seqs[K+i].substring(start,end);			
				ends[K+i] = starts[K+i] + seq.replaceAll("-| ","").length()-1;
				
				LogWriter.stdOutLogln(nm + " " + pad(starts[K+i]+(start>last[K+i] ? 0 : 1),maxPosWidth, false) + 
			             " " + seq + " " + pad(ends[K+i]+1,maxPosWidth, false));
				starts[K+i] = ends[K+i] + 1;

			}
			
			LogWriter.stdOutLogln("\n");
		}
		LogWriter.stdOutLogln("");
		

	}

	private String pad(String nm, int width, boolean left) {
		StringBuffer sb = new StringBuffer(width>nm.length() ? nm : nm.substring(0,width));
		int cols = width-sb.length();
		for (int i=0; i<cols; i++)
				sb.insert(left?sb.length():0 , " ");
		
		return sb.toString();
	}

	private String pad(int n, int width, boolean left) {
		return pad(n+"", width, left);
	}
	
}

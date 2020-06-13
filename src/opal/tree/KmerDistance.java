package opal.tree;

import java.util.Enumeration;
import java.util.Hashtable;
import opal.align.Aligner;
import opal.align.Alignment;

public class KmerDistance extends Distance {
 	
	public KmerDistance (Aligner al) {
		super(al);
	}
	
	public float calcDistance(Alignment A, Alignment B) {
		Hashtable<String, Integer> kA = A.kmers;
		Hashtable<String, Integer> kB = B.kmers;
		
		Enumeration<String> e = kA.keys();
		int cnt = 0;
		int i=0;
		while( e.hasMoreElements() ){
			String key = e.nextElement();
			i++;
			Integer a = kA.get(key);
			Integer b = kB.get(key);
			cnt += Math.min(a==null?0:a.intValue(), b==null?0:b.intValue());
		}
		int denom = Math.min(A.M, B.M) - A.kmerK + 1;
		float dist = 1 - ((float)cnt / denom);
				
		return dist;
	}

}

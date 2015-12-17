package opal.tree;

import opal.IO.SequenceConverter;
import opal.align.Aligner;
import opal.align.Alignment;

public class PctIDDistance extends Distance {
 
	int maxLen = 1;
	
	public PctIDDistance (Aligner al) {
		super(al);
	}

	public float calcDistance(Alignment A, Alignment B) {

		aligner.setAlignments(A, B);
		aligner.align();
		int[][] result = aligner.config.sc.convertPathToIntAlignment(aligner.getPath(), A, B);

		int subs = 0;
		int ids = 0;
		
		for (int i=0; i<result[0].length; i++ ){
			if (result[0][i] != SequenceConverter.GAP_VAL && result[1][i] != SequenceConverter.GAP_VAL ) {
				subs++;
				if (result[0][i] == result[1][i])
					ids++;
			}
		}
		
		float pctID = (float) ids / subs;
						
		return pctID;
	}

}

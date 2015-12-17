package opal.tree;

import com.traviswheeler.libs.LogWriter;

import opal.IO.SequenceConverter;
import opal.align.Aligner;
import opal.align.Alignment;
import opal.exceptions.GenericOpalException;

public class NormCostDistance extends Distance {
	
	public NormCostDistance (Aligner aligner) {
		super(aligner);
	}

	public float calcDistance(Alignment A) {
		return calcDistance(A, null);
	}
	
	public float calcDistance(Alignment A, Alignment B) {
		long totalCost;
		float normCost;
		if (null == B ) {
			if (A.K != 2) {
				LogWriter.stdErrLogln("Internal error. Calculating normcost distance on an invalid number of sequences");
				throw new GenericOpalException("Internal error. Calculating normcost distance on an invalid number of sequences");
			}	
			totalCost = Aligner.calcCost(A.seqs, A.seqIds, aligner.config);
			int cnt = 0;
			for (int i=0; i<A.M; i++){
				if (A.seqs[0][i] != SequenceConverter.GAP_VAL) cnt++; 
				if (A.seqs[1][i] != SequenceConverter.GAP_VAL) cnt++; 
			}			
			normCost = (float) totalCost / cnt;
		} else {
			if (A.K > 1 || B.K > 1) {
				LogWriter.stdErrLogln("Internal error. Calculating normcost distance on an invalid number of sequences");
				throw new GenericOpalException("Internal error. Calculating normcost distance on an invalid number of sequences");
			}	
			aligner.setAlignments(A, B);
			aligner.align();
			totalCost = aligner.getTrueCost();
			normCost = (float) totalCost / (A.M + B.M);
		}
		
		return normCost;
	}

}

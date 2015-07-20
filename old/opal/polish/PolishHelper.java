package opal.polish;

import com.traviswheeler.libs.LogWriter;

import opal.IO.SequenceConverter;
import opal.align.Aligner;
import opal.align.Alignment;

public class PolishHelper {

	long cost;
	int[][] result;
	int[] order;
	Alignment alignment;
	int[] idsAll;
	
	public PolishHelper(Alignment X, Alignment Y, Alignment Z, int[] order, Aligner al) {
		this.order = order;
		al.setAlignments(X, Y);
		al.align();
		int[] idsXY = new int[X.K + Y.K];
		idsAll = new int[X.K + Y.K + Z.K];

		for (int i=0; i<X.K; i++) idsXY[i] = idsAll[i] = X.seqIds[i];
		for (int i=0; i<Y.K; i++) idsXY[X.K + i] = idsAll[X.K + i] = Y.seqIds[i];
		for (int i=0; i<Z.K; i++) idsAll[X.K + Y.K + i] = Z.seqIds[i];
		
		Alignment XY = Alignment.buildNewAlignment(SequenceConverter.convertPathToIntAlignment(al.getPath(), X, Y), idsXY);

		al.setAlignments(XY, Z);
		al.align();
		result = SequenceConverter.convertPathToIntAlignment(al.getPath(), XY, Z);

		//cost = al.getTrueCost();
		
		cost = Aligner.calcCost(XY.seqs, X.K, Y.K, idsXY)
				+ Aligner.calcCost(result, X.K+Y.K, Z.K, idsAll);
		
		setAlignment ();
	}

	public PolishHelper(Alignment X, Alignment Y, int[] order, Aligner al) {
		this.order = order;
		al.setAlignments(X, Y);
		al.align();

		result = SequenceConverter.convertPathToIntAlignment(al.getPath(), X, Y);

		idsAll = new int[X.K + Y.K];
		for (int i=0; i<X.K; i++) idsAll[i] = X.seqIds[i];
		for (int i=0; i<Y.K; i++) idsAll[X.K + i] = Y.seqIds[i];

		cost = Aligner.calcCost(result, X.K, Y.K, idsAll);

		setAlignment ();
	}
	
	public long getCost () {
		return cost;
	}

	public Alignment getAlignment () {
		return alignment;
	}
	
	public void setAlignment () {
		int K = result.length;
	
		int[] ids = new int[result.length];
		int[][] al = new int[result.length][];
		for (int j=0; j<K; j++) {
			al[order[j]] = result[j];
			ids[order[j]] = idsAll[j];
		}		
		alignment = Alignment.buildNewAlignment(al, ids);
	}

}


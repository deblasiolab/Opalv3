package opal.align.shapes;

import opal.IO.SequenceConverter;
import opal.align.Aligner;
import opal.align.ConsistencyModifiers_AllPairs;
import opal.align.ConsistencyModifiers_Pair;

public class ConsistencyShape extends Shape {

	public static ConsistencyModifiers_AllPairs mods;
	
	public ConsistencyShape() {
		// TODO Auto-generated constructor stub
	}

	public ConsistencyShape(Shape s) {
		super(s);
		// TODO Auto-generated constructor stub
	}

	
	 /* procedure: gapBoundaryCost
	 *   Input:
	 *     a, b: column indices into A and B respectively
	 *        s: integer array representation of a shape
	 *   Output: gc, a int representing the sum of pairwise
	 *           gap open/close penalties between A and B (not within A, nor B)
	 *           resulting from appending column a of A and column b of B 
	 *           to shape s (a < 0 indicates deletion, b < 0 indicates insertion).
	 *  This is a O(KL) version.
	 */
	final protected long gapBoundaryCost (int a, int b, Aligner.Direction dir) {
		long cost = 0;
		for (int p = 0; p < K; p++) { 
			for (int q = 0; q < L; q++) {
				ConsistencyModifiers_Pair modpair =  mods.modifiers[p][q];						
				int ii = A.posInUgappedString[p][a==-1?aPos:a];
				int jj = B.posInUgappedString[q][b==-1?bPos:b];				
				
				if (Aligner.Direction.horiz == dir)  {
					if (seqBlocks[p] >= seqBlocks[q+K]  &&  SequenceConverter.GAP_VAL != B.seqs[q][b-1]) { 
						cost += modpair.hGammaOpens[ii][jj];
						if (seqBlocks[p] > seqBlocks[q+K]) //strictly overhanging
							cost += modpair.vGammaCloses[ii][jj-1];
					}
				} else if (Aligner.Direction.vert == dir)  {
					if (seqBlocks[q+K] >= seqBlocks[p]  &&  SequenceConverter.GAP_VAL != A.seqs[p][a-1]) {
						cost += modpair.vGammaOpens[ii][jj];
						if (seqBlocks[q+K] > seqBlocks[p]) // strictly underhanging
							cost += modpair.hGammaCloses[ii-1][jj];
					}
				} else { //if (Aligner.DIAG == dir)
					if ( seqBlocks[p] >= seqBlocks[q+K]  &&  SequenceConverter.GAP_VAL != B.seqs[q][b-1] && SequenceConverter.GAP_VAL == A.seqs[p][a-1]) {
						cost += modpair.hGammaOpens[ii][jj];						
					} else if (seqBlocks[q+K] >= seqBlocks[p]  &&  SequenceConverter.GAP_VAL != A.seqs[p][a-1] && SequenceConverter.GAP_VAL == B.seqs[q][b-1])  {
						cost += modpair.vGammaOpens[ii][jj];						
					}
					if ( seqBlocks[p] > seqBlocks[q+K]  &&  SequenceConverter.GAP_VAL != B.seqs[q][b-1] ) {
						int iii = ii-1;
						if (SequenceConverter.GAP_VAL == A.seqs[p][a-1])
							iii = ii; 
						cost += modpair.vGammaCloses[iii][jj-1];	
					} else if (seqBlocks[q+K] > seqBlocks[p]  &&  SequenceConverter.GAP_VAL != A.seqs[p][a-1] )  {
						int jjj = jj-1;
						if (SequenceConverter.GAP_VAL == B.seqs[q][b-1])
							jjj = jj; 
						cost += modpair.hGammaCloses[ii-1][jjj];	
					}
				}
			}
		}
		return cost;
	}
	
	
	protected long subCost (int a, int b) {
		long cost = 0;
		for (int p = 0; p < K; p++) { 
			for (int q = 0; q < L; q++) {
				ConsistencyModifiers_Pair modpair =  mods.modifiers[p][q];						
				int ii = A.posInUgappedString[p][a==-1?aPos:a];
				int jj = B.posInUgappedString[q][b==-1?bPos:b];

				if (a<0) {//horiz	
					if (SequenceConverter.GAP_VAL != B.seqs[q][b-1]) 
						cost += modpair.hLambdas[ii][jj];			
				} else if (b<0) { //vert
					if (SequenceConverter.GAP_VAL != A.seqs[p][a-1]) 
						cost += modpair.vLambdas[ii][jj];
				} else {				
					if (SequenceConverter.GAP_VAL != B.seqs[q][b-1] && SequenceConverter.GAP_VAL != A.seqs[p][a-1]) {
						cost += modpair.subs[ii][jj];
					} else if (SequenceConverter.GAP_VAL != A.seqs[p][a-1] ) { //vert
						cost += modpair.vLambdas[ii][jj];
					} else if (SequenceConverter.GAP_VAL != B.seqs[q][b-1] ) {
						cost += modpair.hLambdas[ii][jj];
					}

				}
			}
		}
		return cost;
	}

	public void closeGaps () {
		//long cost = 0;
		for (int p = 0; p < K; p++) { 
			for (int q = 0; q < L; q++) {
				ConsistencyModifiers_Pair modpair =  mods.modifiers[p][q];
				int mm = A.posInUgappedString[p][M];
				int nn = B.posInUgappedString[q][N];

				if (seqBlocks[p] > seqBlocks[q+K]) { //vert
					shapeCost += modpair.vGammaCloses[mm][nn];
				} else if (seqBlocks[p] < seqBlocks[q+K]) { //horiz
					shapeCost += modpair.hGammaCloses[mm][nn];
				}					
			}
		}
	}
	
	public static void setMods (ConsistencyModifiers_AllPairs m) {
		mods = m;
	}
}

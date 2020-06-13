package opal.align.shapes;

import opal.IO.SequenceConverter;
import opal.align.Alignment;
import opal.align.ConsistencyModifiers_Pair;
import opal.IO.Configuration;

public class ConsistencyShapeTester extends ShapeTester {

	public ConsistencyShapeTester(long upperBound, long[][] lowerD, long[][] lowerH, long[][] lowerV, Configuration c) {
		super(upperBound, lowerD, lowerH, lowerV, c);
		// TODO Auto-generated constructor stub
	}
/*
	public boolean boundPrune (Shape s) {
		return false;
	}
	*/
	final protected void calcGapBounds(Shape s, long[] H_V_D) {
		//because the consistency approach uses both open and close costs, there are no
		//overcounts to be removed
		H_V_D[0] = H_V_D[1] = H_V_D[2] = 0;
		
		//in fact, it's possible that an open gap will be closed.  If
		// I can confirm it'll happen, I'll tighten the bound with a negative adjustment
		int K = s.K;
		int L = s.L;
		Alignment A = s.A;
		Alignment B = s.B;
		int M = A.M;
		int N = B.M;		
		for (int p = 0; p < K; p++) { 
			for (int q = K; q < K+L; q++) {
				ConsistencyModifiers_Pair modpair =  ConsistencyShape.mods.modifiers[p][q-K];						

				int ii = A.posInUgappedString[p][s.aPos];
				int jj = B.posInUgappedString[q-K][s.bPos];
				int mm = A.posInUgappedString[p][M];
				int nn = B.posInUgappedString[q-K][N];

				//horizontal
				if ( jj!= nn && SequenceConverter.GAP_VAL != B.seqs[q-K][s.bPos] && s.seqBlocks[p]>=s.seqBlocks[q]) {
					//a   or  a   -->    to    <--    -
					//a       -                       a
					//so close the gap coming from the right
					H_V_D[0] -= modpair.hGammaOpens[ii][jj+1];
					if (s.seqBlocks[p]>s.seqBlocks[q]) {
						H_V_D[0] -= modpair.vGammaCloses[ii][jj];
					}
				}
				//vertical
				if ( ii != mm && SequenceConverter.GAP_VAL != A.seqs[p][s.aPos] && s.seqBlocks[q]>=s.seqBlocks[p]) {
					//a   or  -   -->    to    <--    a
					//a       a                       -
					//so close the gap coming from the right
					H_V_D[1] -= modpair.vGammaOpens[ii+1][jj];
					if (s.seqBlocks[q]>s.seqBlocks[p]) {
						H_V_D[1] -= modpair.hGammaCloses[ii][jj];
					}
				}

				//diagonal
				if ( jj!=nn && ii!=mm ) {
					 
					if (SequenceConverter.GAP_VAL != B.seqs[q-K][s.bPos] && SequenceConverter.GAP_VAL == A.seqs[p][s.aPos]
					            && s.seqBlocks[p]>=s.seqBlocks[q]) {
						//a   or  a   -->    to    <--    -
						//a       -                       a
						//so close the gap coming from the right
						H_V_D[2] -= modpair.hGammaOpens[ii][jj+1];
					} else if ( SequenceConverter.GAP_VAL != A.seqs[p][s.aPos] && SequenceConverter.GAP_VAL == B.seqs[q-K][s.bPos]
					            && s.seqBlocks[q]>=s.seqBlocks[p]) {
						//a   or  -   -->    to    <--    a
						//a       a                       -
						//so close the gap coming from the right
						H_V_D[2] -= modpair.vGammaOpens[ii+1][jj];
					}

					if (SequenceConverter.GAP_VAL != B.seqs[q-K][s.bPos] && s.seqBlocks[p]>s.seqBlocks[q]) {
						//  a   -->    to    <--    -   or  a
						//  -                       a       a
						//so close the gap coming from the left
						H_V_D[2] -= modpair.vGammaCloses[ii][jj];
					}
					if (SequenceConverter.GAP_VAL != A.seqs[p][s.aPos] && s.seqBlocks[q]>s.seqBlocks[p]) {
						//  -   -->    to    <--    a   or  a
						//  a                       -       a
						//so close the gap coming from the left
						H_V_D[2] -= modpair.hGammaCloses[ii][jj];
					}				                                                                                                         				
				}
				
			}
		}
		
	}
	
	/* procedure: domGapBound
	 *   Input:  shapes s and t
	 *   Output: cost, a long representing an upperbound on the cost of
	 *           gaps that could possibly be opened or closed in s but not in t,
	 *           over all possible extensions.
	 *  This is a O(KL) version 
	 */
	final protected long domGapBound (Shape s, Shape t) {
		int K = s.K;
		int L = s.L;
		Alignment A = s.A;
		Alignment B = s.B;
		int M = A.M;
		int N = B.M;

		int cost = 0;
		for (int p = 0; p < K; p++) { 
			for (int q = K; q < K+L; q++) {
				ConsistencyModifiers_Pair modpair =  ConsistencyShape.mods.modifiers[p][q-K];						

				int ii = A.posInUgappedString[p][s.aPos];
				int jj = B.posInUgappedString[q-K][s.bPos];
				int mm = A.posInUgappedString[p][M];
				int nn = B.posInUgappedString[q-K][N];

				//opens
				if 	( t.seqBlocks[p]>t.seqBlocks[q] && s.seqBlocks[p]<=s.seqBlocks[q] && ii<mm ) {
					cost += modpair.vGammaOpens[ii+1][jj];
				} else if ( t.seqBlocks[p]<t.seqBlocks[q] && s.seqBlocks[p]>=s.seqBlocks[q] && jj<nn) {
					cost += modpair.hGammaOpens[ii][jj+1];		
				}
				//closes
				if 	( s.seqBlocks[p]>s.seqBlocks[q] && t.seqBlocks[p]<=t.seqBlocks[q] ) {
					cost += modpair.vGammaCloses[ii][jj];
				} else if ( s.seqBlocks[p]<s.seqBlocks[q] && t.seqBlocks[p]>=t.seqBlocks[q] ) {
					cost += modpair.hGammaCloses[ii][jj];		
				}
			}
		}
		return cost;
	}


}

package opal.align;

import com.traviswheeler.libs.LogWriter;
import opal.exceptions.GenericOpalException;
import opal.IO.Configuration;

public class PairSuboptimalityMatricesModified {

	private ConsistencyHeuristicAligner forward, backward; // good for exact pairwise alignments
	
	public Alignment A,B;
	//private int bigNumber = PairwiseAlignmentsContainer.bigBadVal;
	
	public long optCost;
	public long optCost2;
	private int M, P;
	ConsistencyModifiers_Pair modpair;
	LogWriter logger;
	Configuration config;
	
	public PairSuboptimalityMatricesModified (Alignment A, Alignment B, ConsistencyModifiers_Pair modpair, Configuration c) { // single sequence
		this.config = c;
		this.A = A;
		this.B = B;
		this.modpair = modpair;
		
		M = A.M;
		P = B.M;

		forward = new ConsistencyHeuristicAligner( new ConsistencyModifiers_AllPairs(modpair));
		forward.setAlignments(A, B);
		forward.align();
		optCost = forward.getEstimatedCost();
//LogWriter.stdErrLogln("ABC opt cost: " + optCost);
//		char[][] result3 = SequenceConverter.convertPathToCharAlignment(forward.getPath(), SequenceConverter.convertIntsToSeqs(A.seqs), SequenceConverter.convertIntsToSeqs(B.seqs));

		
		//flip modpair
		Alignment Arev = Alignment.buildNewAlignment(config.sc.buildReverseAlignment(A.seqs), A.seqIds, config, A.in);
		Alignment Brev = Alignment.buildNewAlignment(config.sc.buildReverseAlignment(B.seqs), B.seqIds, config, B.in);
		backward = new ConsistencyHeuristicAligner( new ConsistencyModifiers_AllPairs(modpair));  
		backward.setAlignments(Arev,Brev);
		backward.setReverseMode(true);
		backward.align();
		optCost2 = backward.getEstimatedCost();

//		char[][] result4 = SequenceConverter.convertPathToCharAlignment(backward.getPath(), SequenceConverter.convertIntsToSeqs(Arev.seqs), SequenceConverter.convertIntsToSeqs(Brev.seqs));

		// if I don't do this, I'll have double-counted gap opens for the subopt calcs
		forward.V[M][P] -= modpair.vGammaCloses[M][P];
		forward.H[M][P] -= modpair.hGammaCloses[M][P];
		backward.V[M][P] -= modpair.vGammaOpens[1][0];
		backward.H[M][P] -= modpair.hGammaOpens[0][1];
		
		if (optCost != optCost2) {
			LogWriter.stdErrLogln("modified matrix calculation: forward and revers disagree");
			throw new GenericOpalException("");
		}
	}

	
	
	/*
	 * for the reverse parts of the calcs below, the modifier lookups need to 
	 * be flipped, since they were calculated with alignments of 
	 * forward-facing strings. The flipping rule is:
	 *
	 * revSub(i,j) = Sub(M-(i-1), N-(j-1))
	 * revVOpen(i,j) = VClose(M-(i-1), N-j) 
	 * revVExt(i,j) = VExt(M-(i-1), N-j)
	 * revVClose(i,j) = VOpen(M-(i-1), N-j) 
	 * revHOpen(i,j) = HClose(M-i, N-(j-1)) 
	 * revHExt(i,j) = HExt(M-i, N-(j-1))
	 * revHClose(i,j) = HOpen(M-i, N-(j-1)) 
	 * 
	 */
	
		public int getSubsSubopt(int i,int j) {
		/*
		 *         AB 
		 * 
		 *         j  
		 *    -------------- 
		 *   |   o          |    
		 *   |    \         | 
		 * i |     o        | 
		 *   |              | 
		 *   |              | 
		 *   |              | 
		 *    -------------- 
		 * 
		 */	
		long best =  forward.D[i][j] + backward.D[M-(i-1)][P-(j-1)] - modpair.subs[i][j];
		return (int)(best - optCost);
	}



	public int getHExtSubopt(int i,int j) {
		/*
		 *         j
		 *    --------------
		 *   |              |    
		 *   |              | 
		 * i |   o-o        | 
		 *   |              | 
		 *   |              | 
		 *   |              | 
		 *    -------------- 
		 * 
		 */
		long best =  forward.H[i][j] + backward.H[M-i][P-(j-1)] - modpair.hLambdas[i][j];				
		return (int)(best - optCost);
	}

	public int getHOpenSubopt(int i,int j) {
		/*
		 *         j                  j=1
   		 *    --------------         --------------
		 *   |              |       |              |
		 *   |o             |       |              | 
		 *   | \            |       |              |
		 * i |   o-o        |       o-o            |
		 *   |              |       |              |
		 *   |              |       |              |
		 *    --------------         --------------
		 * 
		 */

		long best;
		if (j==1) {
			best =  forward.H[i][1] + backward.H[M-i][P] - modpair.hLambdas[i][j];
		} else {
			best =  forward.D[i][j-1] + backward.H[M-i][P-(j-1)] + modpair.hGammaOpens[i][j];
		}
		return (int)(best - optCost);
	}

	public int getHCloseSubopt(int i,int j) {
		/*
		 *         j                             j=P
   		 *    --------------         --------------
		 *   |              |       |              |
		 *   |              |       |              | 
		 *   |              |       |              |
		 * i |   o-o        |       |            o-o
		 *   |      \       |       |              |
		 *   |       o      |       |              |
		 *    --------------         --------------
		 * 
		 */

		long best;
		if (j==P) {
			best =  forward.H[i][P] + backward.H[M-i][1] - modpair.hLambdas[i][j];
		} else {
			best =  forward.H[i][j] + backward.D[M-i][P-j] + modpair.hGammaCloses[i][j];
		}
		return (int)(best - optCost);
	}

	
	public int getVExtSubopt(int i,int j) {
		/*
		 *         j
		 *    --------------
		 *   |              |    
		 *   |     o        | 
		 *   |     |        | 
		 * i |     o        | 
		 *   |              | 
		 *   |              | 
		 *    -------------- 
		 * 
		 */

		long best =  forward.V[i][j] + backward.V[M-(i-1)][P-j] - modpair.vLambdas[i][j];				
		return (int)(best - optCost);
	}

	public int getVOpenSubopt(int i,int j) {
		/*
		 *       j                        j
   		 *    --------------         -----o--------
		 *   |              |       |     |        |
		 *   | o            |   i=1 |     o        | 
		 *   |  \           |       |              |
		 *   |   o          |       |              |
		 *   |   |          |       |              |
		 * i |   o          |       |              |
		 *   |              |       |              |
		 *    --------------         --------------
		 * 
		 */

		long best;
		if (i==1) {
			best =  forward.V[1][j] + backward.V[M][P-j] - modpair.vLambdas[i][j];
		} else {
			best =  forward.D[i-1][j] + backward.V[M-(i-1)][P-j] + modpair.vGammaOpens[i][j];
		}
		return (int)(best - optCost);
	}

	public int getVCloseSubopt(int i,int j) {
		/*
		 *         j                       j      
   		 *    --------------         --------------
		 *   |              |       |              |
		 *   |     o        |       |              | 
		 *   |     |        |       |              |
		 * i |     o        |       |              |
		 *   |      \       |       |      o       |
		 *   |       o      |       |      |       |
		 *    --------------   i=M   ------o-------
		 * 
		 */

		long best;
		if (i==M) {
			best =  forward.V[M][j] + backward.V[1][P-j] - modpair.vLambdas[i][j];
		} else {
			best =  forward.V[i][j] + backward.D[M-i][P-j] + modpair.vGammaCloses[i][j];
		}
		return (int)(best - optCost);
	}
	
}

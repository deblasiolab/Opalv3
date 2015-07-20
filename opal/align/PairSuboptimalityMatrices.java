package opal.align;

import com.traviswheeler.libs.LogWriter;
import opal.IO.SequenceConverter;
import opal.exceptions.GenericOpalException;
import opal.IO.Configuration;

public class PairSuboptimalityMatrices {

	private CompressedAlignmentGraph forward, backward;
		
	public int[] nearlyOptimal_horizStart; //for each A_i, the first position B_j that is "good"
	public int[] nearlyOptimal_horizEnd; //for each A_i, the last position B_j that is "good"
	//public int[] nearlyOptimal_vertStart; //for each B_j, the first position A_i that is "good"
	//public int[] nearlyOptimal_vertEnd;
	public Alignment A,C;
	
	private long bigNumber = PairwiseAlignmentsContainer.bigBadVal;
	
	public long optCost;
	public long optCost2;
	private int M, P;
	public static int delta = 50;
		
	char[][] result;
	
	Configuration config;
	
	public PairSuboptimalityMatrices (Alignment A, Alignment B, Configuration c) { // single sequence
		this(A,B,true,c);
	}
	
	public PairSuboptimalityMatrices (Alignment A, Alignment B, boolean compress, Configuration c) { // single sequence
		this.A = A;
		this.C = B;
		this.config = c;
		
		M = A.M;
		P = B.M;	

		Alignment Arev = Alignment.buildNewAlignment(config.sc.buildReverseAlignment(A.seqs), A.seqIds, config);
		Alignment Brev = Alignment.buildNewAlignment(config.sc.buildReverseAlignment(B.seqs), B.seqIds, config);

		PairAligner_SplitGamma forwardAligner, backwardAligner;
		forwardAligner = new PairAligner_SplitGamma(A,B);
		forwardAligner.align();
		optCost = forwardAligner.getEstimatedCost();

		
		result = config.sc.convertPathToCharAlignment(forwardAligner.getPath(), config.sc.convertIntsToSeqs(A.seqs), config.sc.convertIntsToSeqs(B.seqs));
		
		backwardAligner = new PairAligner_SplitGamma(Arev,Brev);
		backwardAligner.align();


		nearlyOptimal_horizStart = new int[M+1];
		nearlyOptimal_horizEnd = new int[M+1];
		//nearlyOptimal_vertStart = new int[P+1];
		//nearlyOptimal_vertEnd = new int[P+1];


		for (int i=0; i<=M; i++) {
			nearlyOptimal_horizStart[i] = -1;
		}
//		for (int i=0; i<=P; i++) {
//			nearlyOptimal_vertStart[i] = -1;
//		}
		
		
		//These graphs all had the gap-close cost added during alignment. Because of how the graphs are used for recovery, I remove it here
		forwardAligner.H[M][P] -= config.gammaTerm/2;
		forwardAligner.V[M][P] -= config.gammaTerm/2;
		backwardAligner.H[M][P] -= config.gammaTerm/2;
		backwardAligner.V[M][P] -= config.gammaTerm/2;
		
		for (int i=0; i<M; i++) {
			for (int j=0; j<P; j++) {

				if (i==0 && j==0) {
					nearlyOptimal_horizStart[0] = 0;
//					nearlyOptimal_vertStart[0] = 0;
					continue;
				}
				
				//compare these ...
				long best = bigNumber ;
				if (i>0 && j>0) {
					int sigma = config.cost.costs[A.seqs[0][i-1]][C.seqs[0][j-1]]; 
					best =  forwardAligner.D[i][j] + backwardAligner.D[M-(i-1)][P-(j-1)] - sigma;
				}
				if (j>0)
					best = Math.min(best, forwardAligner.H[i][j] + backwardAligner.H[M-i][P-(j-1)] - ((i==0||i==M)?config.lambdaTerm:config.lambda));
				if (i>0)
					best = Math.min(best, forwardAligner.V[i][j] + backwardAligner.V[M-(i-1)][P-j] - ((j==0||j==P)?config.lambdaTerm:config.lambda));

				
				if (best <= optCost + delta) {
					nearlyOptimal_horizEnd[i] = j;
					if (nearlyOptimal_horizStart[i] == -1) nearlyOptimal_horizStart[i] = j;

//					nearlyOptimal_vertEnd[j] = i;					
//					if (nearlyOptimal_vertStart[j] == -1) nearlyOptimal_vertStart[j] = i;
				}
			}
		}

		for (int j=0; j<P; j++) {
			long best =            forwardAligner.D[M][j] + config.lambdaTerm * (P-j) + config.gammaTerm;
			best = Math.min(best, forwardAligner.H[M][j] + config.lambdaTerm * (P-j) + config.gammaTerm/2);
			best = Math.min(best, forwardAligner.V[M][j] + config.lambdaTerm * (P-j) + config.gammaTerm + config.gamma/2);
			if (best <= optCost + delta) {
				nearlyOptimal_horizEnd[M] = j;
				if (nearlyOptimal_horizStart[M] == -1) nearlyOptimal_horizStart[M] = j;

//				nearlyOptimal_vertEnd[j] = M;					
//				if (nearlyOptimal_vertStart[j] == -1) nearlyOptimal_vertStart[j] = M;
			}		
		}
		
		for (int i=0; i<M; i++) {
			long best =            forwardAligner.D[i][P] + config.lambdaTerm * (M-i) + config.gammaTerm;
			best = Math.min(best, forwardAligner.H[i][P] + config.lambdaTerm * (M-i) + config.gammaTerm + config.gamma/2);
			best = Math.min(best, forwardAligner.V[i][P] + config.lambdaTerm * (M-i) + config.gammaTerm/2);
			if (best <= optCost + delta) {
				nearlyOptimal_horizEnd[i] = P;
				if (nearlyOptimal_horizStart[i] == -1) nearlyOptimal_horizStart[i] = P;

//				nearlyOptimal_vertEnd[P] = i;					
//				if (nearlyOptimal_vertStart[P] == -1) nearlyOptimal_vertStart[P] = i;
			}		
		}

		nearlyOptimal_horizEnd[M] = P;
		if (nearlyOptimal_horizStart[M] == -1) nearlyOptimal_horizStart[M] = P;

		forward = new CompressedAlignmentGraph(forwardAligner, nearlyOptimal_horizStart, nearlyOptimal_horizEnd, false, compress ); 
		backward = new CompressedAlignmentGraph(backwardAligner, nearlyOptimal_horizStart, nearlyOptimal_horizEnd, true, compress );
		
//		nearlyOptimal_vertEnd[P] = M;					
//		if (nearlyOptimal_vertStart[P] == -1) nearlyOptimal_vertStart[P] = M;

	}

	
	public int getSubsSubopt(int i,int k) {
		/*
		 *         AC                     BC
		 * 
		 *         k                      k
		 *    --------------         --------------
		 *   |   o          |       |   o          |    
		 *   |    \         |       |    \         | 
		 * i |     o        |     j |     o        | 
		 *   |              |       |              | 
		 *   |              |       |              | 
		 *   |              |       |              | 
		 *    --------------         --------------
		 * 
		 */	
		int sigma = config.cost.costs[A.seqs[0][i-1]][C.seqs[0][k-1]]; 
		long best =  forward.D[i][k-forward.shift[i]] + backward.D[M-(i-1)][P-(k-1)-backward.shift[M-(i-1)]] - sigma;
		
		return (int)(best - optCost);
	}


	public int getExtSubopt(int i,int k, boolean isLeftGraph) {
		
		int minval = -1;
		if ( k > 0 ) {
			minval =  Math.min( 
					getExtSubopt1(i, k, isLeftGraph), 
					getExtSubopt2(i, k, isLeftGraph)
					);
		} else { //k==0.  
			//This is only valid if j==0, then only in one case:
			if ( isLeftGraph ||  i==0 ) {
				minval = getExtSubopt2(i, k, isLeftGraph);
			} else {	
				LogWriter.stdErrLogln("suprising attempt to calc ExtSubopt with k==0, i>0 on right graph");
				throw new GenericOpalException("");
			}
		}
		return minval;
		
	}
	
	
	private int getExtSubopt1(int i,int k, boolean isLeftGraph) {
		/*
		 *         k                      k
		 *    --------------         --------------
		 *   |              |       |              |    
		 *   |   o          |       |              | 
		 *   |    \         |     j |   o-o        | 
		 * i |     o        |       |              | 
		 *   |              |       |              | 
		 *   |              |       |              | 
		 *    --------------         --------------
		 * 
		 */

		long best=0;

		if (isLeftGraph) {
			//i and k are both > 0  ... from function call
			int sigma = config.cost.costs[A.seqs[0][i-1]][C.seqs[0][k-1]]; 
			best =  forward.D[i][k-forward.shift[i]] + backward.D[M-(i-1)][P-(k-1)-backward.shift[M-(i-1)]] - sigma;			
		} else {
			int j = i; int N=M; // just to sync the code below with the graph above. In other words: read "i" as "j" in the right graph above
			//k > 0  ,   j can be 0
			
			best =  forward.H[j][k-forward.shift[j]] + backward.H[N-j][P-(k-1)-backward.shift[N-j]] - ((j==0||j==N)?config.lambdaTerm:config.lambda);
		}
						
		return (int)(best - optCost);
	}

	
	
	private int getExtSubopt2(int i,int k, boolean isLeftGraph) {
		/*Standard case: i>0, 0<k<P;    if j==0, just ignore first diag
		 *         k                      k
		 *    --------------         --------------
		 *   |              |       |   o          |    
		 *   |     o        |       |    \         | 
		 *   |     |        |     j |   o-o-o      | 
		 * i |     o        |       |      \       | 
		 *   |              |       |       o      | 
		 *   |              |       |              | 
		 *    --------------         --------------
		 * 
		 */

		 /* Special case: k=0;  j must be 0  - this needs to be guaranteed by calling function.
		  *                                 if j>0, cost should be delta*mult
		 *   k                      
		 *    --------------        o-o------------
		 *   |              |       |\             |    
		 *   o              |       | o            | 
		 *   |              |       |              | 
		 * i o              |       |              | 
		 *   |              |       |              | 
		 *   |              |       |              | 
		 *    --------------         --------------
		 * 
		 */

		 /* Special case: k=P;  j must be N  - this needs to be guaranteed by calling function.
		  *                                 if j<N, cost should be delta*mult
		 *                  k                      
		 *    --------------         --------------
		 *   |              |       |              |    
		 *   |              o       |              | 
		 *   |              |       |              | 
		 * i |              o       |              | 
		 *   |              |       |            o | 
		 *   |              |       |             \| 
		 *    --------------         ------------o-o
		 * 
		 */
		
		long best=0;
		
		if (isLeftGraph) {
			//i > 0   ... from function call  (except in case where it's coopted for alignment of A-B
			best =  forward.V[i][k-forward.shift[i]] + backward.V[M-(i-1)][P-k-backward.shift[M-(i-1)]] - ((k==0||k==P)?config.lambdaTerm:config.lambda);
		} else {
			int j = i; int N=M; // just to sync the code below with the graph above. In other words: read "i" as "j" in the right graph above
			//k > 0   ... from function call
			if (j==N && k==P) {
				best = Math.min(forward.D[N][P-forward.shift[N]], forward.H[N][P-forward.shift[N]] + config.gammaTerm/2);
			} else if (j==0 && k==0){
				best = Math.min(backward.D[N][P-backward.shift[N]], backward.H[N][P-backward.shift[N]] + config.gammaTerm/2);
			} else if (k==0 || k==P) {
				/*
				 *  This is a case where seq C isn't able to inform the alignment of positions A_i and B_j
				 *           i
				 *    A  XXXYYZZ 
				 *    C  XXX----
				 *    B  XXXY--Z
				 *          j
				 *          
				 *    what to do here? I think it makes sense to require that the k-range 
				 *    always include P-1 if it includes P. Then this won't be the largest number.
				 *    I do that in the calling code ... and then make sure this scenario never returns a value
				 */
				LogWriter.stdErrLogln("Calling getExtSubopt2 with k =" + k + " and j = " + j + ".  Should never happen");
				throw new GenericOpalException("");
				//return bigNumber;
			} else { // 0<k<P
				if (j<N ) {					
					best = forward.H[j][k-forward.shift[j]] + backward.H[N-j][P-k-backward.shift[N-j]];					
					best = Math.min(best,forward.H[j][k-forward.shift[j]] + (i==0?config.gammaTerm:config.gamma)/2 + backward.D[N-j][P-k-backward.shift[N-j]]);
					if (j>0) {
						best = Math.min(best, forward.D[j][k-forward.shift[j]] + backward.D[N-j][P-k-backward.shift[N-j]]);
						best = Math.min(best, forward.D[j][k-forward.shift[j]] + config.gamma/2 + backward.H[N-j][P-k-backward.shift[N-j]]);
					}
				} else {  //j==N , k<P
					//diag:horiz  at bottom boundary			
					best = forward.D[N][k-forward.shift[N]] + config.gammaTerm/2 + backward.H[0][P-k-backward.shift[0]];
					//	horiz:horiz
					best = Math.min(best, forward.H[N][k-forward.shift[N]] + backward.H[0][P-k-backward.shift[0]]);
				}
			}
		}
		return (int)(best - optCost);
	}

	

	public int getOpenSubopt(int i,int k, boolean isLeftGraph) {
		int minval = Math.min( getOpenSubopt1(i, k, isLeftGraph), getOpenSubopt2(i, k, isLeftGraph));
		return minval;
	}
	
	
	public int getOpenSubopt1(int i,int k, boolean isLeftGraph) {
		/*
		 *         k                        k
		 *    --------------         --------------
		 *   |              |       |              |    
		 *   |   o          |       |   o          | 
		 *   |    \         |       |    \         |
		 * i |     o        |     j |     o-o      | 
		 *   |              |       |              | 
		 *   |              |       |              | 
		 *    --------------         --------------
		 * 
		 */
		
		/*
		 * Special case when j==0,  value will only be used if i == 1
		 *         k                        k
		 *    --------------      j  -----o-o------
		 *   |              |       |              |    
		 *   |   o          |       |              | 
		 *   |    \         |       |              |
		 * i |     o        |       |              | 
		 *   |              |       |              | 
		 *   |              |       |              | 
		 *    --------------         --------------
		 * 
		 */
		
		long best=0;
		if (isLeftGraph) {
			//i>0  and k>1 from function call (unless j==0, then i==1 && k>0)
			int sigma = config.cost.costs[A.seqs[0][i-1]][C.seqs[0][k-1]]; 
			best =  forward.D[i][k-forward.shift[i]] + backward.D[M-(i-1)][P-(k-1)-backward.shift[M-(i-1)]] - sigma;
		} else {
			int j = i; int N=M; // just to sync the code below with the graph above. In other words: read "i" as "j" in the right graph above
			//k > 1  from function call, unless i==1 && j==0;
			if (j>0) 
				best =  forward.D[j][k-1-forward.shift[j]] + (j==N?config.gammaTerm:config.gamma)/2 + backward.H[N-j][P-(k-1)-backward.shift[N-j]] ;
			else 
				best =  forward.H[0][k-forward.shift[0]] + backward.H[N][P-(k-1)-backward.shift[N]] - config.lambdaTerm ;
		}
						
		return (int)(best - optCost);
	}

	public int getOpenSubopt2(int i,int k, boolean isLeftGraph) {
		/*
		 * when i>1 (j>0, 0<k<P)
		 *         k                       k
		 *    --------------         --------------
		 *   |              |       |              |    
		 *   |   o          |       |    o         |    
		 *   |    \         |       |     \        |    
		 *   |     o        |     j |      o-o     | 
		 *   |     |        |       |       \      | 
		 * i |     o        |       |        o     | 
		 *   |              |       |              |  
		 *    --------------         --------------
		 * 
		 */

		/* Special case with j==N, 0<k<=P.  (k==P allowed)
		 *         k                         k
		 *    --------------         --------------
		 *   |              |       |              |    
		 *   |   same       |       |              |    
		 *   |    as        |       |              |    
		 *   |   above      |       |              | 
		 *   |    or        |       |              | 
		 *   |   below      |       |      o       | 
		 *   |              |       |       \      |  
		 *    --------------       j --------o-----
		 */

		/* Special case with i == 1; j>0, 0<k<P.
		 *         k                       k
		 *    -----o--------         --------------
		 *   |     |        |       |              |    
		 * 1 |     o        |       |    o         |    
		 *   |              |       |     \        |    
		 *   |              |     j |      o-o     | 
		 *   |              |       |       \      | 
		 *   |              |       |        o     | 
		 *   |              |       |              |  
		 *    --------------         --------------
		 * 
		 */

		/* Special case with i == 1, j==0, 0=<k<P  (k==0 allowed)
		 * 
		 *         k                    k       
		 *    -----o--------        ----o-o--------
		 *   |     |        |       |    \         |    
		 * 1 |     o        |       |     o        |    
		 *   |              |       |              |    
		 *   |              |       |              | 
		 *   |              |       |              | 
		 *   |              |       |              | 
		 *   |              |       |              |  
		 *    --------------         --------------
		 * 
		 */
		
		long best;
		
		if (isLeftGraph) {
			int gamma_open = (k==P ? config.gammaTerm :config.gamma)/2;
			if (i>1) {
				//check this at boundaries
				best =  forward.D[i-1][k-forward.shift[i-1]] + gamma_open + backward.V[M-(i-1)][P-k-backward.shift[M-(i-1)]];
			} else { //if ( i==1 ) {
				best = forward.V[1][k-forward.shift[1]] + backward.V[M][P-k-backward.shift[M]] - (k==0||k==P?config.lambdaTerm:config.lambda);
			}
		} else {
			int j = i; int N=M; // just to sync the code below with the graph above. In other words: read "i" as "j" in the right graph above
			best = bigNumber;

			if (j==0 /*k<P from calling function*/) { // should really only be run if i==1
				int sigma = config.cost.costs[A.seqs[0][1-1]][C.seqs[0][(k+1)-1]];
				best = forward.D[1][k+1-forward.shift[1]] + backward.D[N][P-k-backward.shift[N]] -  sigma;
				best = Math.min(best, forward.H[0][k+1-forward.shift[0]] + backward.H[N][P-k-backward.shift[N]] - config.lambdaTerm);
			} else if (j==N /*k>0 from calling function*/) {
				int sigma = config.cost.costs[A.seqs[0][N-1]][C.seqs[0][k-1]];
				best = forward.D[N][k-forward.shift[N]] + backward.D[0][P-(k-1)-backward.shift[0]] - sigma;
			} else { // (k<P from calling function)
				best = forward.D[j][k-forward.shift[j]] + config.gamma/2 + backward.H[N-j][P-k-backward.shift[N-j]] ;
				best = Math.min(best, forward.D[j][k-forward.shift[j]] +  backward.D[N-j][P-k-backward.shift[N-j]]);
			}
		}
						
		return (int)(best - optCost);
	}
	
	
	public int getOpenSubopt_special(int i, int k_prime, int k, boolean isLeftGraph) {
		/*
		 *          k'        k                  k'      k  
		 *      -------------------         -------------------
		 *     |                   |       |                   |    
		 *     |                   |       |   o               |    
		 * i-1 |    o-o-...-o      |       |    \              |    
		 *     |             \     |     j |     o-o-...-o     | 
		 * i   |              o    |       |                   | 
		 *     |                   |       |                   | 
		 *     |                   |       |                   |  
		 *      -------------------         -------------------
		 * 
		 */
		long best;
		if (isLeftGraph) {
			//k_prime > 0  and  i > 0 ... from function call
			int ext_cnt = k - k_prime - 2;
			int ext =  (i==1?config.lambdaTerm : config.lambda) * ext_cnt;
			int close =  (i==1?config.gammaTerm : config.gamma)/2;
				
			best = forward.H[i-1][k_prime+1-forward.shift[i-1]] + ext + close  + backward.D[M-(i-1)][P-(k-1)-backward.shift[M-(i-1)]];
		} else {
			int j = i; int N=M; // just to sync the code below with the graph above. In other words: read "i" as "j" in the right graph above
			//j and k_prime > 0 ... from function call
			int ext_cnt = k - k_prime - 1;
			int ext =  (j==N ? config.lambdaTerm : config.lambda) * ext_cnt;
			int open = (j==N ? config.gammaTerm :config.gamma)/2;
		
			best = forward.D[j][k_prime-forward.shift[j]] + ext + open  + backward.H[N-j][P-(k-1)-backward.shift[N-j]];
		}
		return (int)(best - optCost);
	}


	
	
	
	
	public int getCloseSubopt(int i,int k, boolean isLeftGraph) {
		int minval = Math.min(getCloseSubopt1(i, k, isLeftGraph),  getCloseSubopt2(i, k, isLeftGraph));
		return minval;
	}
	
	
	public int getCloseSubopt1(int i,int k, boolean isLeftGraph) {
		/*
		 *         k                      k
		 *    --------------         --------------
		 *   |              |       |              |    
		 *   |   o          |     j |   o-o        | 
		 *   |    \         |       |      \       |
		 * i |     o        |       |       o      | 
		 *   |              |       |              | 
		 *   |              |       |              | 
		 *    --------------         --------------
		 * 
		 */

		/*  special case, j==N.  Only actually used when i==M
		 *         k                         k
		 *    --------------         --------------
		 *   |              |       |              |    
		 *   |              |       |              | 
		 *   |              |       |              |
		 *   |              |       |              | 
		 *   |   o          |       |              | 
		 *   |    \         |       |              | 
		 *    -----o--------         ------o-o-----
		 * 
		 */		
			
		long best=0;
		if (isLeftGraph) {
			//i and k are both > 0   ... from function call
			int sigma = config.cost.costs[A.seqs[0][i-1]][C.seqs[0][k-1]]; 
			best =  forward.D[i][k-forward.shift[i]] + backward.D[M-(i-1)][P-(k-1)-backward.shift[M-(i-1)]] - sigma;
		} else {
			int j = i; int N=M; // just to sync the code below with the graph above. In other words: read "i" as "j" in the right graph above
			//k<P from function call, unless j==N && i==M
			if (j<N ) {
				best =  forward.H[j][k-forward.shift[j]] + (j==0?config.gammaTerm:config.gamma)/2 + backward.D[N-j][P-k-backward.shift[N-j]] ;
			} else { // j==N; only makes sense when i==M, so calling function filters for that.
				best = forward.H[N][k-forward.shift[N]] + backward.H[0][P-(k-1)-backward.shift[0]] - config.lambdaTerm ;
			}			
		
		}
						
		return (int)(best - optCost);
	}

	public int getCloseSubopt2(int i,int k, boolean isLeftGraph) {


		/*
		 * when i>1 (j>0, 0<k<P)
		 *       k                         k  
		 *    --------------         --------------
		 *   |              |       |              |    
		 *   |   o          |       |    o         |    
		 *   |   |          |       |     \        |    
		 * i |   o          |     j |    o-o       | 
		 *   |    \         |       |       \      | 
		 *   |     o        |       |        o     | 
		 *   |              |       |              |  
		 *    --------------         --------------
		 * 
		 */

		/* Special case with j==0. i>0 (i==M ok),  0<=k<P (k==0 ok)
		 * 
		 *       k                       k   
		 *    --------------         ----o---------
		 *   |              |       |     \        |    
		 *   |   same       |       |      o       |    
		 *   |    as        |       |              |    
		 * i |   above      |       |              | 
		 *   |    or        |       |              | 
		 *   |   below      |       |              | 
		 *   |              |       |              |  
		 *    --------------         --------------
		 * 
		 */

		
		/* Special case with i==M; j<N, k<P. (j,k>0)
		 *         k                       k
		 *    --------------         --------------
		 *   |              |       |              |    
		 *   |              |       |    o         |    
		 *   |              |       |     \        |    
		 *   |              |     j |    o-o       | 
		 *   |              |       |       \      | 
		 *   |     o        |       |        o     | 
		 *   |     |        |       |              |  
		 *    -----o--------         --------------
		 * 
		 */

		/* Special case with i == M, j==N, 0<k<=P (k==P allowed)
		 * 
		 *            k                        k
		 *    --------------         --------------
		 *   |              |       |              |    
		 *   |              |       |              |    
		 *   |              |       |              |    
		 *   |              |       |              | 
		 *   |              |       |              | 
		 *   |        o     |       |        o     | 
		 *   |        |     |       |         \    |  
		 *    --------o-----         --------o-o---
		 * 
		 */
		
		long best=0;
		
		if (isLeftGraph) {
			if (i<M) {
				best =  forward.V[i][k-forward.shift[i]] + (k==0?config.gammaTerm:config.gamma)/2 + backward.D[M-i][P-k-backward.shift[M-i]];
			} else { //if ( i==M ) {     k==P only if i==M
				best =  forward.V[M][k-forward.shift[M]] + backward.V[1][P-k-backward.shift[1]] - (k==0||k==P?config.lambdaTerm:config.lambda);
			}
		} else {
			int j = i; int N=M; // just to sync the code below with the graph above. In other words: read "i" as "j" in the right graph above
			best = bigNumber;
			if (j==0 ) {  //ok for k to be 0
				int sigma = config.cost.costs[A.seqs[0][1-1]][C.seqs[0][(k+1)-1]];
				best = forward.D[1][k+1-forward.shift[1]] + backward.D[N][P-k-backward.shift[N]] - sigma ;
			} else if (j==N ) { // ok for k==P, but not 0  
				int sigma = config.cost.costs[A.seqs[0][N-1]][C.seqs[0][k-1]];
				best = forward.D[N][k-forward.shift[N]] + backward.D[1][P-(k-1)-backward.shift[1]] - sigma ;
				best = Math.min(best, forward.H[N][k-forward.shift[N]] +  backward.H[0][P-(k-1)-backward.shift[0]] - config.lambdaTerm);
			} else { // 0<j<N , 0<k<P
				best = forward.D[j][k-forward.shift[j]] + backward.D[N-j][P-k-backward.shift[N-j]] ;
				best = Math.min(best, forward.H[j][k-forward.shift[j]] + config.gamma/2 + backward.D[N-j][P-k-backward.shift[N-j]]);
			}	
		}
		
		return (int)(best - optCost);
	}
	
	
	public int getCloseSubopt_special(int i, int k_prime, int k, boolean isLeftGraph) {
		/*
		 *           k'     k-1                  k'      k  
		 *      -------------------         -------------------
		 *     |                   |       |                   |    
		 *     |                   |       |                   |    
		 *     |   o               |       |                   |    
		 *     |    \              |     j |   o-o-...-o       | 
		 * i   |     o-o-...-o     |       |            \      | 
		 *     |                   |       |             o     | 
		 *     |                   |       |                   |  
		 *      -------------------         -------------------
		 * 
		 */
		long best = 0;
		if (isLeftGraph) {
			//k_prime > 0  and  i > 0 ... from function call
			int ext_cnt = k - k_prime - 2;
			int ext =  (i==M?config.lambdaTerm : config.lambda) * ext_cnt;
			int close =  (i==M?config.gammaTerm : config.gamma)/2;
			
			best = forward.D[i][k_prime-forward.shift[i]] + ext + close  + backward.H[M-i][P-(k-2)-backward.shift[M-i]];
		} else {
			int j = i; int N=M; // just to sync the code below with the graph above. In other words: read "i" as "j" in the right graph above
			//j and k_prime both > 0 ... from function call
			int ext_cnt = k - k_prime - 1;
			int ext =  (j==0 ? config.lambdaTerm : config.lambda) * ext_cnt;
			int open = (j==0 ? config.gammaTerm :config.gamma)/2;
						
			best = forward.H[j][k_prime-forward.shift[j]] + ext + open  + backward.D[N-j][P-(k-1)-backward.shift[N-j]];
		}
		
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

		long best =  forward.H[i][j-forward.shift[i]] + backward.H[M-i][P-(j-1)-backward.shift[M-i]] - ((i==0||i==M)?config.lambdaTerm:config.lambda);				
		
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

		long best=0;
		if (j==1) {
			best =  forward.H[i][1-forward.shift[i]] + backward.H[M-i][P-backward.shift[M-i]] - ((i==0||i==M)?config.lambdaTerm:config.lambda);
		} else {
			best =  forward.D[i][j-1-forward.shift[i]] + backward.H[M-i][P-(j-1)-backward.shift[M-i]] + (i==M?config.gammaTerm:config.gamma);;
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

		long best=0;
		if (j==P) {
			best =  forward.H[i][P-forward.shift[i]] + backward.H[M-i][1-backward.shift[M-i]] - ((i==0||i==M)?config.lambdaTerm:config.lambda);
		} else {
			best =  forward.H[i][j-forward.shift[i]] + backward.D[M-i][P-j-backward.shift[M-i]] + (i==0?config.gammaTerm:config.gamma);;
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

		long best =  forward.V[i][j-forward.shift[i]] + backward.V[M-(i-1)][P-j-backward.shift[M-(i-1)]] - ((j==0||j==P)?config.lambdaTerm:config.lambda);				

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
			best =  forward.V[1][j-forward.shift[1]] + backward.V[M][P-j-backward.shift[M]] - ((j==0||j==P)?config.lambdaTerm:config.lambda);
		} else {
			best =  forward.D[i-1][j-forward.shift[i-1]] + backward.V[M-(i-1)][P-j-backward.shift[M-(i-1)]] + (j==P?config.gammaTerm:config.gamma);
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
			best =  forward.V[M][j-forward.shift[M]] + backward.V[1][P-j-backward.shift[1]] - ((j==0||j==P)?config.lambdaTerm:config.lambda);
		} else {
			best =  forward.V[i][j-forward.shift[i]] + backward.D[M-i][P-j-backward.shift[M-i]] + (j==0?config.gammaTerm:config.gamma);
		}

		return (int)(best - optCost);
	}
	
	
	public PairAligner_SplitGamma getAlignment() {
		return null;
	}
	
	public static void setDelta (int d) {
		delta = d;
	}


}

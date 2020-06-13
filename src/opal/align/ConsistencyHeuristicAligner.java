package opal.align;

import java.util.ArrayList;
import java.util.Collections;

import com.traviswheeler.libs.LogWriter;
import opal.IO.SequenceConverter;
import opal.exceptions.GenericOpalException;

public class ConsistencyHeuristicAligner extends Aligner {

	boolean testing = false;
	//boolean testing = true;

	
	PairwiseAlignmentsContainer alignmentContainer;	
	ConsistencyModifiers_AllPairs mods;
	boolean reverseMode = false;
	/*
	 * If reverseMode == true, then all the modifier lookups need to 
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
	
	
	public ConsistencyHeuristicAligner ( Aligner al) {
		super(al);
		if (al instanceof ConsistencyAligner) {
			alignmentContainer = ((ConsistencyAligner)al).getAlignmentContainer();
			mods =  ((ConsistencyAligner)al).getConsistencyModifiers();
		} else {
			LogWriter.stdErrLogln("Unexpected attempt to create a ConsistencyHeuristicAligner");
			throw new GenericOpalException("");
		}
	}

	public ConsistencyHeuristicAligner( ConsistencyModifiers_AllPairs mods) {
		super(true);
		this.mods = mods;
	}

	public ConsistencyHeuristicAligner( LogWriter logger ) {
		super(true);
	}
	
	public ConsistencyHeuristicAligner( Alignment A, Alignment B) {
		this(A,B,true);
	}

	public ConsistencyHeuristicAligner( boolean pess) {
		super(pess);
	}

	public ConsistencyHeuristicAligner( Alignment alA, Alignment alB, boolean pess) {
		super(alA,alB,pess);
	}
	
	protected void initialize () {
		//the +1 thing here shifts over the characters to a 1-based counting method
		H = new long[M+1][N+1];
		V = new long[M+1][N+1];
		D = new long[M+1][N+1];				
	}		

	public void cleanup () {
		//the +1 thing here shifts over the characters to a 1-based counting method
		H = null;
		V = null;
		D = null;				
	}		
	
	protected void fillBoundary() {
		int bigNumber = Integer.MAX_VALUE/2 ; // half, to ensure no overflow from the next row/column 

		V[0][0] = H[0][0] = D[0][0] = 0; 
		for (int j=1; j<=N; j++) {
			int tmp=0;
			for (int p=0; p<K; p++) {
				int mm = A.posInUgappedString[p][M];
				for (int q=0; q<L; q++) {
					int jj = B.posInUgappedString[q][j];
					int nn = B.posInUgappedString[q][N];
					
					if (B.seqs[q][j-1] != SequenceConverter.GAP_VAL) {
						if (reverseMode) tmp += mods.modifiers[p][q].hLambdas[mm][nn-(jj-1)]; 
						else tmp += mods.modifiers[p][q].hLambdas[0][jj]; 
						
						if  ( jj == 1) {
							if (reverseMode)  tmp += mods.modifiers[p][q].hGammaCloses[mm][nn];
							else tmp += mods.modifiers[p][q].hGammaOpens[0][1];
						}
				    } 					
				}
			}
			H[0][j] = H[0][j-1] + tmp; 
			D[0][j] = V[0][j] = bigNumber;
		}
				
		
		
		for (int i=1; i<=M; i++) {
			int tmp=0;
			for (int p=0; p<K; p++) {
				int ii = A.posInUgappedString[p][i];
				int mm = A.posInUgappedString[p][M];
				for (int q=0; q<L; q++) {					
					int nn = B.posInUgappedString[q][N];
					
					if (A.seqs[p][i-1] != SequenceConverter.GAP_VAL) {
						if (reverseMode) tmp += mods.modifiers[p][q].vLambdas[mm-(ii-1)][nn]; 
						else tmp += mods.modifiers[p][q].vLambdas[ii][0]; 						
						if  ( ii == 1) { 
							if (reverseMode)  tmp += mods.modifiers[p][q].vGammaCloses[mm][nn];
							else tmp += mods.modifiers[p][q].vGammaOpens[1][0];
						}
				    } 					
				}
			}
			V[i][0] = V[i-1][0] + tmp; 
			D[i][0] = H[i][0] = bigNumber;			
		}		
	}
	
	
	protected void fillTable() {
		
		
		for (int i=1; i<=M; i++){
			for (int j=1; j<=N; j++){
				
				int tmpHH=0, tmpHD=0, tmpHV=0;
				int tmpDH=0, tmpDD=0, tmpDV=0;
				int tmpVH=0, tmpVD=0, tmpVV=0;
				int tmpH=0, tmpD=0, tmpV=0;

				for (int p=0; p<K; p++) {
					int ii = A.posInUgappedString[p][i];
					int mm = A.posInUgappedString[p][M];
					
					for (int q=0; q<L; q++) {
						ConsistencyModifiers_Pair modpair =  mods.modifiers[p][q];
						int jj = B.posInUgappedString[q][j];
						int nn = B.posInUgappedString[q][N];
	
						
						/* fill V[i][j] */
						if (A.seqs[p][i-1] != SequenceConverter.GAP_VAL) {
							//   x
 							//   -
							if (reverseMode) tmpV += modpair.vLambdas[mm-(ii-1)][nn-jj];
							else tmpV += modpair.vLambdas[ii][jj];   							
							
							if (i>1) { // there's no legit cost for D(0,j) or V(0,j)
								//from V[i-1][j]
								//    -  x
								//    -  -
								if  (isPessimistic && A.seqs[p][i-2] == SequenceConverter.GAP_VAL) {
									
									if (reverseMode) tmpVV += modpair.vGammaCloses[mm-(ii-1)][nn-jj];
									else tmpVV += modpair.vGammaOpens[ii][jj];   
										
									if ( jj>0 ) {
										if (reverseMode) tmpVV += modpair.hGammaOpens[mm-(ii-1)][nn-(jj-1)];
										else tmpVV += modpair.hGammaCloses[ii-1][jj];   
									}
								}
							
								//D[i-1][j]
								//    o x    or    - x
								//    x -          - -							
								if  ( B.seqs[q][j-1] != SequenceConverter.GAP_VAL || 
									(isPessimistic  && A.seqs[p][i-2] == SequenceConverter.GAP_VAL) ) { 
									if (reverseMode) tmpDV += modpair.vGammaCloses[mm-(ii-1)][nn-jj];
									else  tmpDV += modpair.vGammaOpens[ii][jj];
								}
								//    - x 
								//    o - 
								if ( A.seqs[p][i-2] == SequenceConverter.GAP_VAL && 
										(B.seqs[q][j-1] != SequenceConverter.GAP_VAL || (jj>0 && isPessimistic) ) ) {
									if (reverseMode) tmpDV += modpair.hGammaOpens[mm-(ii-1)][nn-(jj-1)];
									else tmpDV += modpair.hGammaCloses[ii-1][jj];
								}
							}
							
							//H[i-1][j]
							//      - x     (A[i-1] must be "-", b/c from horiz)
							//      o - 
							if  (B.seqs[q][j-1] != SequenceConverter.GAP_VAL || (jj>0 && isPessimistic) ) { 
								if (reverseMode) tmpHV += modpair.vGammaCloses[mm-(ii-1)][nn-jj] + modpair.hGammaOpens[mm-(ii-1)][nn-(jj-1)];
								else tmpHV += modpair.vGammaOpens[ii][jj] + modpair.hGammaCloses[ii-1][jj];
							} else if ( jj==0 && i==1 ) {
								if (reverseMode) tmpHV += modpair.vGammaCloses[mm-(ii-1)][nn-jj];
								else tmpHV += modpair.vGammaOpens[ii][jj];
							}
							
						}
						if (i==M && j==N){
							if (A.seqs[p][M-1] != SequenceConverter.GAP_VAL) {
								if (reverseMode) tmpV += modpair.vGammaOpens[mm-(ii-1)][nn-jj];
								else  tmpV += modpair.vGammaCloses[ii][jj];
							} else if (isPessimistic) {
								if (reverseMode) tmpV += Math.max(modpair.vGammaOpens[mm-(ii-1)][nn-jj],modpair.hGammaOpens[mm-ii][nn-(jj-1)]);
								else  tmpV += Math.max(modpair.vGammaCloses[ii][jj],modpair.hGammaCloses[ii][jj]);
							}
						}
						
						
						/* fill H[i][j] */
						if (B.seqs[q][j-1] != SequenceConverter.GAP_VAL) {
							//  -
							//  x
							if (reverseMode) tmpH += modpair.hLambdas[mm-ii][nn-(jj-1)];
							else tmpH += modpair.hLambdas[ii][jj];
							
							//from V[i][j-1]
							//  o -    
							//  - x    (B[j-1] will be "-" b/c from vert)
							if  (A.seqs[p][i-1] != SequenceConverter.GAP_VAL || (ii>0 && isPessimistic) ) { 
								if (reverseMode) tmpVH += modpair.hGammaCloses[mm-ii][nn-(jj-1)] + modpair.vGammaOpens[mm-(ii-1)][nn-(jj-1)];
								else tmpVH += modpair.hGammaOpens[ii][jj] + modpair.vGammaCloses[ii][jj-1];
							} else if ( ii==0 && j==1 ) {
								if (reverseMode) tmpVH += modpair.hGammaCloses[mm-ii][nn-(jj-1)];
								else tmpVH += modpair.hGammaOpens[ii][jj];
							}
							
							if (j>1) { // there's no legit cost from H(i,0) or D(i,0)
								//from D[i][j-1]
								//   x -    or   - -
								//   o x         - x
								if  (A.seqs[p][i-1] != SequenceConverter.GAP_VAL || 
										(isPessimistic && B.seqs[q][j-2] == SequenceConverter.GAP_VAL) ) { 
									if (reverseMode) tmpDH += modpair.hGammaCloses[mm-ii][nn-(jj-1)];
									else tmpDH += modpair.hGammaOpens[ii][jj];
								}
								//   o -   
								//   - x
								if (B.seqs[q][j-2] == SequenceConverter.GAP_VAL && 
										(A.seqs[p][i-1] != SequenceConverter.GAP_VAL || (ii>0 && isPessimistic)) ) {
									if (reverseMode) tmpDH += modpair.vGammaOpens[mm-(ii-1)][nn-(jj-1)];
									else tmpDH += modpair.vGammaCloses[ii][jj-1];
								}
							
								//from H[i-1][j]
								//  - -
								//  - x
								if  (isPessimistic && B.seqs[q][j-2] == SequenceConverter.GAP_VAL) {
									if (reverseMode) tmpHH += modpair.hGammaCloses[mm-ii][nn-(jj-1)];
									else tmpHH += modpair.hGammaOpens[ii][jj];
									if (ii>0) {
										if (reverseMode) tmpHH += modpair.vGammaOpens[mm-(ii-1)][nn-(jj-1)];
										else tmpHH += modpair.vGammaCloses[ii][jj-1];
									}
								}
							}
						}						
						if (i==M && j==N){
							if (B.seqs[q][N-1] != SequenceConverter.GAP_VAL) {
								if (reverseMode) tmpH += modpair.hGammaOpens[mm-ii][nn-(jj-1)];
								else tmpH += modpair.hGammaCloses[ii][jj];
							} else if ( isPessimistic) {
								if (reverseMode) tmpH += Math.max(modpair.hGammaOpens[mm-ii][nn-(jj-1)],modpair.vGammaOpens[mm-(ii-1)][nn-jj]);
								else tmpH += Math.max(modpair.hGammaCloses[ii][jj], modpair.vGammaCloses[ii][jj]);
							}

						}
						
						
						/* fill D[i][j] */
						//universal sub/ext cost
						if (A.seqs[p][i-1] != SequenceConverter.GAP_VAL && B.seqs[q][j-1] != SequenceConverter.GAP_VAL) {
							if (reverseMode) tmpD += modpair.subs[mm-(ii-1)][nn-(jj-1)];
							else tmpD += modpair.subs[ii][jj];
						} else if (A.seqs[p][i-1] != SequenceConverter.GAP_VAL) { // and B is a gap
							if (reverseMode) tmpD +=  modpair.vLambdas[mm-(ii-1)][nn-jj];
							else tmpD +=  modpair.vLambdas[ii][jj];
						} else if (B.seqs[q][j-1] != SequenceConverter.GAP_VAL) { // and A is a gap
							if (reverseMode) tmpD +=  modpair.hLambdas[mm-ii][nn-(jj-1)];
							else tmpD +=  modpair.hLambdas[ii][jj];
						}


						//from V[i-1][j-1]
						if (i>1) { // V(0,j) isn't valid
							if (A.seqs[p][i-1] != SequenceConverter.GAP_VAL && B.seqs[q][j-1] == SequenceConverter.GAP_VAL 
									&& (isPessimistic && A.seqs[p][i-2] == SequenceConverter.GAP_VAL)) {
								//  - x
								//  - -   B[j-1] is "-" b/c/ from vert
								if (reverseMode) tmpVD += modpair.vGammaCloses[mm-(ii-1)][nn-jj];
								else tmpVD += modpair.vGammaOpens[ii][jj];
								if (jj>0) {
									if (reverseMode) tmpVD += modpair.hGammaOpens[mm-(ii-1)][nn-(jj-1)];
									else tmpVD += modpair.hGammaCloses[ii-1][jj];
								}
						    } else if (A.seqs[p][i-1] == SequenceConverter.GAP_VAL && B.seqs[q][j-1] != SequenceConverter.GAP_VAL) {
						    	//  o -
						    	//  - x  B[j-1] is "-" b/c/ from vert
						    	if ( A.seqs[p][i-2] != SequenceConverter.GAP_VAL || isPessimistic) { 
						    		if (reverseMode) tmpVD += modpair.hGammaCloses[mm-ii][nn-(jj-1)];
						    		else tmpVD += modpair.hGammaOpens[ii][jj];
						    		if (ii>0){
						    			if (reverseMode) tmpVD += modpair.vGammaOpens[mm-(ii-1)][nn-(jj-1)];
						    			else tmpVD += modpair.vGammaCloses[ii][jj-1];
						    		}
						    	}
							} else if (A.seqs[p][i-1] != SequenceConverter.GAP_VAL && B.seqs[q][j-1] != SequenceConverter.GAP_VAL) {
								if (A.seqs[p][i-2] != SequenceConverter.GAP_VAL ) {	
									//  x x
									//  - x
									if (reverseMode) tmpVD +=  modpair.vGammaOpens[mm-(ii-1-1)][nn-(jj-1)];
									else tmpVD +=  modpair.vGammaCloses[ii-1][jj-1];
								} else if (isPessimistic && (ii>1 || jj>1)) {
									//  - x
									//  - x
									if (ii>1 && jj>1){
										if (reverseMode) tmpVD += Math.max(modpair.vGammaOpens[mm-(ii-1-1)][nn-(jj-1)], modpair.hGammaOpens[mm-(ii-1)][nn-(jj-1-1)]);
										else tmpVD += Math.max(modpair.vGammaCloses[ii-1][jj-1], modpair.hGammaCloses[ii-1][jj-1]);
									} else if (ii==1 && jj>1) { 
										if (reverseMode) tmpVD += modpair.hGammaOpens[mm-(ii-1)][nn-(jj-1-1)];
										else tmpVD += modpair.hGammaCloses[ii-1][jj-1];
									} else { // ii>1, jj==1 
										if (reverseMode) tmpVD += modpair.vGammaOpens[mm-(ii-1-1)][nn-(jj-1)];
										else tmpVD += modpair.vGammaCloses[ii-1][jj-1];
									}
								}
							}		
						}
						
						//from D[i-1][j-1]
						if (i>1 && j>1) {  //D(0,j) and D(i,0) are not legit
							if (A.seqs[p][i-1] != SequenceConverter.GAP_VAL && B.seqs[q][j-1] == SequenceConverter.GAP_VAL) { 
									if ( B.seqs[q][j-2] != SequenceConverter.GAP_VAL || 
											(isPessimistic && A.seqs[p][i-2] == SequenceConverter.GAP_VAL)) {
										//  o x    or   - x
										//  x -         - -
										if (reverseMode) tmpDD += modpair.vGammaCloses[mm-(ii-1)][nn-jj] ;
										else tmpDD += modpair.vGammaOpens[ii][jj] ;
									}
									if ( (A.seqs[p][i-2] == SequenceConverter.GAP_VAL && B.seqs[q][j-2] != SequenceConverter.GAP_VAL)
										|| (jj>0 && isPessimistic && A.seqs[p][i-2] == SequenceConverter.GAP_VAL && B.seqs[q][j-2] == SequenceConverter.GAP_VAL)) {
										//  - x    or   - x
										//  x -         - -
										if (reverseMode) tmpDD += modpair.hGammaOpens[mm-(ii-1)][nn-(jj-1)];
										else tmpDD += modpair.hGammaCloses[ii-1][jj];
									}
						    } else if (A.seqs[p][i-1] == SequenceConverter.GAP_VAL && B.seqs[q][j-1] != SequenceConverter.GAP_VAL) {
						    	if ( A.seqs[p][i-2] != SequenceConverter.GAP_VAL || 
										(isPessimistic && B.seqs[q][j-2] == SequenceConverter.GAP_VAL)) {
									//  x -    or   - -
									//  o x         - x
						    		if (reverseMode) tmpDD += modpair.hGammaCloses[mm-ii][nn-(jj-1)] ;
						    		else tmpDD += modpair.hGammaOpens[ii][jj] ;
								}
								if ( (A.seqs[p][i-2] != SequenceConverter.GAP_VAL && B.seqs[q][j-2] == SequenceConverter.GAP_VAL)
									|| (ii>0 && isPessimistic && A.seqs[p][i-2] == SequenceConverter.GAP_VAL && B.seqs[q][j-2] == SequenceConverter.GAP_VAL)) {
									//  x -    or   - -
									//  - x         - x
									if (reverseMode) tmpDD += modpair.vGammaOpens[mm-(ii-1)][nn-(jj-1)];
									else tmpDD += modpair.vGammaCloses[ii][jj-1];
								}
						    } else if (A.seqs[p][i-1] != SequenceConverter.GAP_VAL && B.seqs[q][j-1] != SequenceConverter.GAP_VAL) {
						    	if (A.seqs[p][i-2] != SequenceConverter.GAP_VAL && B.seqs[q][j-2] == SequenceConverter.GAP_VAL) { 
									//  x x    
									//  - x    
						    		if (reverseMode) tmpDD +=  modpair.vGammaOpens[mm-(ii-1-1)][nn-(jj-1)];
						    		else tmpDD +=  modpair.vGammaCloses[ii-1][jj-1];
						    	} else if (A.seqs[p][i-2] == SequenceConverter.GAP_VAL && B.seqs[q][j-2] != SequenceConverter.GAP_VAL) {
						    		//  - x    
									//  x x    
						    		if (reverseMode) tmpDD +=  modpair.hGammaOpens[mm-(ii-1)][nn-(jj-1-1)];
						    		else tmpDD +=  modpair.hGammaCloses[ii-1][jj-1];
						    	} else if (isPessimistic && (ii>1 || jj>1) && A.seqs[p][i-2] == SequenceConverter.GAP_VAL && B.seqs[q][j-2] == SequenceConverter.GAP_VAL) {
						    		//  - x  
									//  - x  
									if (ii>1 && jj>1) {
										if (reverseMode) tmpDD += Math.max(modpair.vGammaOpens[mm-(ii-1-1)][nn-(jj-1)], modpair.hGammaOpens[mm-(ii-1)][nn-(jj-1-1)]);
										else tmpDD += Math.max(modpair.vGammaCloses[ii-1][jj-1], modpair.hGammaCloses[ii-1][jj-1]);
									} else if (ii==1 && jj>1){ 
										if (reverseMode) tmpDD += modpair.hGammaOpens[mm-(ii-1)][nn-(jj-1-1)];
										else tmpDD += modpair.hGammaCloses[ii-1][jj-1];
									}else { // ii>1, jj==1 
										if (reverseMode) tmpDD += modpair.vGammaOpens[mm-(ii-1-1)][nn-(jj-1)];
										else tmpDD += modpair.vGammaCloses[ii-1][jj-1];
									}
						    	}

						    }
						}						
						
						//from H[i-1][j-1]
						if (j>1) { // H(0,j) isn't valid

							if (A.seqs[p][i-1] != SequenceConverter.GAP_VAL && B.seqs[q][j-1] == SequenceConverter.GAP_VAL) {
								//  - x   
								//  o -   
								if (B.seqs[q][j-2] != SequenceConverter.GAP_VAL || isPessimistic) { 
									if (reverseMode) tmpHD += modpair.vGammaCloses[mm-(ii-1)][nn-jj];
									else tmpHD += modpair.vGammaOpens[ii][jj];
									if (jj>1) {
										if (reverseMode) tmpHD += modpair.hGammaOpens[mm-(ii-1)][nn-(jj-1)];
										else tmpHD += modpair.hGammaCloses[ii-1][jj];
									}
								}
						    } else if (A.seqs[p][i-1] == SequenceConverter.GAP_VAL && B.seqs[q][j-1] != SequenceConverter.GAP_VAL
						    		&& isPessimistic && B.seqs[q][j-2] == SequenceConverter.GAP_VAL) {
									//  - -   
									//  - x   
						    		if (reverseMode) tmpHD += modpair.hGammaCloses[mm-ii][nn-(jj-1)];
						    		else tmpHD += modpair.hGammaOpens[ii][jj];
							    	if (ii>1) {
							    		if (reverseMode) tmpHD += modpair.vGammaOpens[mm-(ii-1)][nn-(jj-1)];
							    		else tmpHD += modpair.vGammaCloses[ii][jj-1];
							    	}
							} else if (A.seqs[p][i-1] != SequenceConverter.GAP_VAL && B.seqs[q][j-1] != SequenceConverter.GAP_VAL) {
								if (B.seqs[q][j-2] != SequenceConverter.GAP_VAL) {
									//  - x   
									//  x x   
									if (reverseMode) tmpHD +=  modpair.hGammaOpens[mm-(ii-1)][nn-(jj-1-1)];
									else tmpHD +=  modpair.hGammaCloses[ii-1][jj-1];
								} else if (isPessimistic && (ii>1 || jj>1) ) {
									//  - x   
									//  - x   
									if (ii>1 && jj>1) {
										if (reverseMode) tmpHD += Math.max(modpair.vGammaOpens[mm-(ii-1-1)][nn-(jj-1)], modpair.hGammaOpens[mm-(ii-1)][nn-(jj-1-1)]);
										else tmpHD += Math.max(modpair.vGammaCloses[ii-1][jj-1], modpair.hGammaCloses[ii-1][jj-1]);
									} else if (ii==1 && jj>1) { 
										
										if (reverseMode) tmpHD += modpair.hGammaOpens[mm-(ii-1)][nn-(jj-1-1)];
										else tmpHD += modpair.hGammaCloses[ii-1][jj-1];
									} else { // ii>1, jj==1 
										if (reverseMode) tmpHD += modpair.vGammaOpens[mm-(ii-1-1)][nn-(jj-1)];
										else tmpHD += modpair.vGammaCloses[ii-1][jj-1];
									}
								}
							}
						}
						
						if (i==1 && j==1) { // special case for beginning of matrix
							int tmp = 0;
							if (A.seqs[p][i-1] != SequenceConverter.GAP_VAL && B.seqs[q][j-1] == SequenceConverter.GAP_VAL) {
								//opening a gap
								if (reverseMode) tmp += modpair.vGammaCloses[mm-(ii-1)][nn-jj] ;
								else tmp += modpair.vGammaOpens[ii][jj] ;
							} else if (A.seqs[p][i-1] == SequenceConverter.GAP_VAL && B.seqs[q][j-1] != SequenceConverter.GAP_VAL) {
								if (reverseMode) tmp += modpair.hGammaCloses[mm-ii][nn-(jj-1)] ;
					    		else tmp += modpair.hGammaOpens[ii][jj] ;
							}
							tmpDD += tmp;
							tmpVD += tmp;
							tmpHD += tmp;
						} 
						
						if (i==M && j==N){
							if (A.seqs[p][M-1] != SequenceConverter.GAP_VAL && B.seqs[q][N-1] == SequenceConverter.GAP_VAL) { 
								if (reverseMode) tmpD += modpair.vGammaOpens[1][0];
								else tmpD += modpair.vGammaCloses[mm][nn];
							} else if (A.seqs[p][M-1] == SequenceConverter.GAP_VAL && B.seqs[q][N-1] != SequenceConverter.GAP_VAL) { 
								if (reverseMode) tmpD += modpair.hGammaOpens[0][1];
								else tmpD += modpair.hGammaCloses[mm][nn];
							} else if (isPessimistic && A.seqs[p][M-1] == SequenceConverter.GAP_VAL && B.seqs[q][N-1] == SequenceConverter.GAP_VAL) {
								if (reverseMode) tmpD += Math.max( modpair.hGammaOpens[0][1], modpair.vGammaOpens[1][0]);
								else tmpD += Math.max( modpair.hGammaCloses[mm][nn], modpair.vGammaCloses[mm][nn]);
							}
						}

					}
				}
				V[i][j] = tmpV + Math.min(V[i-1][j] + tmpVV,
						Math.min(D[i-1][j]+tmpDV, H[i-1][j]+tmpHV)) ;
				H[i][j] = tmpH + Math.min(V[i][j-1] + tmpVH,
						Math.min(D[i][j-1]+tmpDH, H[i][j-1]+tmpHH)) ;
				D[i][j] = tmpD + Math.min(V[i-1][j-1] + tmpVD,
						Math.min(D[i-1][j-1]+tmpDD, H[i-1][j-1]+tmpHD)) ;
				
			}
		}
	}	

	
	protected boolean recover () {
		
		path = new ArrayList<Direction>(2*Math.max(M,N));
		int i=M, j=N;
		Direction dir, nextdir=null;
		long cost = min(H[i][j], V[i][j], D[i][j]);
		estimatedCost = cost;
		
		dir =  (cost == H[i][j]) ?  Direction.horiz : (cost == V[i][j]) ?  Direction.vert :  Direction.diag;

		while (i>0 && j>0) {			
			if (Direction.vert == dir){			
				float tmpHV=0, tmpDV=0, tmpVV=0, tmpV=0;				
				for (int p=0; p<K; p++) {
					int ii = A.posInUgappedString[p][i];
					int mm = A.posInUgappedString[p][M];

					for (int q=0; q<L; q++) {
						ConsistencyModifiers_Pair modpair =  mods.modifiers[p][q];						
						int jj = B.posInUgappedString[q][j];
						int nn = B.posInUgappedString[q][N];
						
						if (A.seqs[p][i-1] != SequenceConverter.GAP_VAL) {

							if (reverseMode) tmpV += modpair.vLambdas[mm-(ii-1)][nn-jj];
							else tmpV += modpair.vLambdas[ii][jj];   							
							if (i>1) { // there's no legit cost for D(0,j) or V(0,j)
								if  (isPessimistic && A.seqs[p][i-2] == SequenceConverter.GAP_VAL) {
									if (reverseMode) tmpVV += modpair.vGammaCloses[mm-(ii-1)][nn-jj];
									else tmpVV += modpair.vGammaOpens[ii][jj];   
									if ( jj>0 ) {
										if (reverseMode) tmpVV += modpair.hGammaOpens[mm-(ii-1)][nn-(jj-1)];
										else tmpVV += modpair.hGammaCloses[ii-1][jj];   
									}
								}								
								if  ( B.seqs[q][j-1] != SequenceConverter.GAP_VAL || 
									(isPessimistic  && A.seqs[p][i-2] == SequenceConverter.GAP_VAL) ) { 
									if (reverseMode) tmpDV += modpair.vGammaCloses[mm-(ii-1)][nn-jj];
									else  tmpDV += modpair.vGammaOpens[ii][jj];
								}
								if ( A.seqs[p][i-2] == SequenceConverter.GAP_VAL && 
										(B.seqs[q][j-1] != SequenceConverter.GAP_VAL || (jj>0 && isPessimistic) ) ) {
									if (reverseMode) tmpDV += modpair.hGammaOpens[mm-(ii-1)][nn-(jj-1)];
									else tmpDV += modpair.hGammaCloses[ii-1][jj];
								}
							}
							
							if  (B.seqs[q][j-1] != SequenceConverter.GAP_VAL || (jj>0 && isPessimistic) ) { 
								if (reverseMode) tmpHV += modpair.vGammaCloses[mm-(ii-1)][nn-jj] + modpair.hGammaOpens[mm-(ii-1)][nn-(jj-1)];
								else tmpHV += modpair.vGammaOpens[ii][jj] + modpair.hGammaCloses[ii-1][jj];
							} else if ( jj==0 && i==1 ) {
								if (reverseMode) tmpHV += modpair.vGammaCloses[mm-(ii-1)][nn-jj];
								else tmpHV += modpair.vGammaOpens[ii][jj];
							}
						}
						if (i==M && j==N){
							if (A.seqs[p][M-1] != SequenceConverter.GAP_VAL) {
								if (reverseMode) tmpV += modpair.vGammaOpens[mm-(ii-1)][nn-jj];
								else  tmpV += modpair.vGammaCloses[ii][jj];
							} else if (isPessimistic) {
								if (reverseMode) tmpV += Math.max(modpair.vGammaOpens[mm-(ii-1)][nn-jj],modpair.hGammaOpens[mm-ii][nn-(jj-1)]);
								else  tmpV += Math.max(modpair.vGammaCloses[ii][jj],modpair.hGammaCloses[ii][jj]);
							}
						}
					}
				}
				
				if ( cost == H[i-1][j] + tmpV + tmpHV) {
					nextdir = Direction.horiz;
					if (testing) LogWriter.stdOutLogln("hv(" + i + "," + j+ "): " + (tmpV+tmpHV));
				} else if ( cost == V[i-1][j] + tmpV + tmpVV)	 {			
					nextdir = Direction.vert;
					if (testing) LogWriter.stdOutLogln("vv (" + i + "," + j+ "): "  + (tmpV+tmpVV));
				} else if ( cost == D[i-1][j] + tmpV + tmpDV) {					
					nextdir = Direction.diag;
					if (testing) LogWriter.stdOutLogln("dv(" + i + "," + j+ "): " + (tmpV+tmpDV));
				} else {
					LogWriter.stdErrLogln("no cost source found in dir == Direction.vert");
					return false;					
				}
				i--;
			} else if (Direction.horiz == dir){			
				float tmpHH=0, tmpDH=0, tmpVH=0, tmpH=0;				
				for (int p=0; p<K; p++) {
					int ii = A.posInUgappedString[p][i];
					int mm = A.posInUgappedString[p][M];

					for (int q=0; q<L; q++) {
						ConsistencyModifiers_Pair modpair =  mods.modifiers[p][q];						
						int jj = B.posInUgappedString[q][j];
						int nn = B.posInUgappedString[q][N];
						
						if (B.seqs[q][j-1] != SequenceConverter.GAP_VAL) {
							if (reverseMode) tmpH += modpair.hLambdas[mm-ii][nn-(jj-1)];
							else tmpH += modpair.hLambdas[ii][jj];
							if  (A.seqs[p][i-1] != SequenceConverter.GAP_VAL || (ii>0 && isPessimistic) ) { 
								if (reverseMode) tmpVH += modpair.hGammaCloses[mm-ii][nn-(jj-1)] + modpair.vGammaOpens[mm-(ii-1)][nn-(jj-1)];
								else tmpVH += modpair.hGammaOpens[ii][jj] + modpair.vGammaCloses[ii][jj-1];
							} else if ( ii==0 && j==1 ) {
								if (reverseMode) tmpVH += modpair.hGammaCloses[mm-ii][nn-(jj-1)];
								else tmpVH += modpair.hGammaOpens[ii][jj];
							}
							if (j>1) { // there's no legit cost from H(i,0) or D(i,0)
								if  (A.seqs[p][i-1] != SequenceConverter.GAP_VAL || 
										(isPessimistic && B.seqs[q][j-2] == SequenceConverter.GAP_VAL) ) {
									if (reverseMode) tmpDH += modpair.hGammaCloses[mm-ii][nn-(jj-1)];
									else tmpDH += modpair.hGammaOpens[ii][jj];
								}
								if (B.seqs[q][j-2] == SequenceConverter.GAP_VAL && 
										(A.seqs[p][i-1] != SequenceConverter.GAP_VAL || (ii>0 && isPessimistic)) ){
									if (reverseMode) tmpDH += modpair.vGammaOpens[mm-(ii-1)][nn-(jj-1)];
									else tmpDH += modpair.vGammaCloses[ii][jj-1];
								}
								if  (isPessimistic && B.seqs[q][j-2] == SequenceConverter.GAP_VAL) {
									if (reverseMode) tmpHH += modpair.hGammaCloses[mm-ii][nn-(jj-1)];
									else tmpHH += modpair.hGammaOpens[ii][jj];
									if (ii>0) {
										if (reverseMode) tmpHH += modpair.vGammaOpens[mm-(ii-1)][nn-(jj-1)];
										else tmpHH += modpair.vGammaCloses[ii][jj-1];
									}
								}
							}
						}						
						if (i==M && j==N){
							if (B.seqs[q][N-1] != SequenceConverter.GAP_VAL) {
								if (reverseMode) tmpH += modpair.hGammaOpens[mm-ii][nn-(jj-1)];
								else tmpH += modpair.hGammaCloses[ii][jj];
							} else if ( isPessimistic) {
								if (reverseMode) tmpH += Math.max(modpair.hGammaOpens[mm-ii][nn-(jj-1)],modpair.vGammaOpens[mm-(ii-1)][nn-jj]);
								else tmpH += Math.max(modpair.hGammaCloses[ii][jj], modpair.vGammaCloses[ii][jj]);
							}

						}

					}
				}
												
				if ( cost == H[i][j-1] + tmpH + tmpHH ) {
					nextdir = Direction.horiz;
					if (testing) LogWriter.stdOutLogln("hh(" + i + "," + j+ "): " + (tmpH+tmpHH));
				} else if ( cost == V[i][j-1] + tmpH + tmpVH ) {
					nextdir = Direction.vert;
					if (testing) LogWriter.stdOutLogln("vh(" + i + "," + j+ "): " + (tmpH+tmpVH));
				} else if ( cost == D[i][j-1] + tmpH + tmpDH) {		
					nextdir = Direction.diag;
					if (testing) LogWriter.stdOutLogln("dh(" + i + "," + j+ "): " + (tmpH+tmpDH));
				} else {
					LogWriter.stdErrLogln("no cost source found in dir == Direction.horiz");
					return false;					
				}
				j--;

			} else if (Direction.diag == dir ) {		
				float tmpHD=0, tmpDD=0, tmpVD=0, tmpD=0;				
				for (int p=0; p<K; p++) {
					int ii = A.posInUgappedString[p][i];
					int mm = A.posInUgappedString[p][M];

					for (int q=0; q<L; q++) {
						ConsistencyModifiers_Pair modpair =  mods.modifiers[p][q];						
						int jj = B.posInUgappedString[q][j];
						int nn = B.posInUgappedString[q][N];
						

						if (A.seqs[p][i-1] != SequenceConverter.GAP_VAL && B.seqs[q][j-1] != SequenceConverter.GAP_VAL) {
							if (reverseMode) tmpD += modpair.subs[mm-(ii-1)][nn-(jj-1)];
							else tmpD += modpair.subs[ii][jj];
						} else if (A.seqs[p][i-1] != SequenceConverter.GAP_VAL) { // and B is a gap
							if (reverseMode) tmpD +=  modpair.vLambdas[mm-(ii-1)][nn-jj];
							else tmpD +=  modpair.vLambdas[ii][jj];
						} else if (B.seqs[q][j-1] != SequenceConverter.GAP_VAL) { // and A is a gap
							if (reverseMode) tmpD +=  modpair.hLambdas[mm-ii][nn-(jj-1)];
							else tmpD +=  modpair.hLambdas[ii][jj];
						}
						if (i>1) { // V(0,j) isn't valid
							if (A.seqs[p][i-1] != SequenceConverter.GAP_VAL && B.seqs[q][j-1] == SequenceConverter.GAP_VAL 
									&& (isPessimistic && A.seqs[p][i-2] == SequenceConverter.GAP_VAL)) {
								if (reverseMode) tmpVD += modpair.vGammaCloses[mm-(ii-1)][nn-jj];
								else tmpVD += modpair.vGammaOpens[ii][jj];
								if (jj>0) {
									if (reverseMode) tmpVD += modpair.hGammaOpens[mm-(ii-1)][nn-(jj-1)];
									else tmpVD += modpair.hGammaCloses[ii-1][jj];
								}
						    } else if (A.seqs[p][i-1] == SequenceConverter.GAP_VAL && B.seqs[q][j-1] != SequenceConverter.GAP_VAL) {
						    	if ( A.seqs[p][i-2] != SequenceConverter.GAP_VAL || isPessimistic) { 
						    		if (reverseMode) tmpVD += modpair.hGammaCloses[mm-ii][nn-(jj-1)];
						    		else tmpVD += modpair.hGammaOpens[ii][jj];
						    		if (ii>0){
						    			if (reverseMode) tmpVD += modpair.vGammaOpens[mm-(ii-1)][nn-(jj-1)];
						    			else tmpVD += modpair.vGammaCloses[ii][jj-1];
						    		}
						    	}
							} else if (A.seqs[p][i-1] != SequenceConverter.GAP_VAL && B.seqs[q][j-1] != SequenceConverter.GAP_VAL) {
								if (A.seqs[p][i-2] != SequenceConverter.GAP_VAL ) {	
									if (reverseMode) tmpVD +=  modpair.vGammaOpens[mm-(ii-1-1)][nn-(jj-1)];
									else tmpVD +=  modpair.vGammaCloses[ii-1][jj-1];
								} else if (isPessimistic && (ii>1 || jj>1)) {
									if (ii>1 && jj>1){
										if (reverseMode) tmpVD += Math.max(modpair.vGammaOpens[mm-(ii-1-1)][nn-(jj-1)], modpair.hGammaOpens[mm-(ii-1)][nn-(jj-1-1)]);
										else tmpVD += Math.max(modpair.vGammaCloses[ii-1][jj-1], modpair.hGammaCloses[ii-1][jj-1]);
									} else if (ii==1 && jj>1) { 
										if (reverseMode) tmpVD += modpair.hGammaOpens[mm-(ii-1)][nn-(jj-1-1)];
										else tmpVD += modpair.hGammaCloses[ii-1][jj-1];
									} else { // ii>1, jj==1 
										if (reverseMode) tmpVD += modpair.vGammaOpens[mm-(ii-1-1)][nn-(jj-1)];
										else tmpVD += modpair.vGammaCloses[ii-1][jj-1];
									}
								}
							}		
						}
						
						if (i>1 && j>1) {  //D(0,j) and D(i,0) are not legit
							if (A.seqs[p][i-1] != SequenceConverter.GAP_VAL && B.seqs[q][j-1] == SequenceConverter.GAP_VAL) { 
									if ( B.seqs[q][j-2] != SequenceConverter.GAP_VAL || 
											(isPessimistic && A.seqs[p][i-2] == SequenceConverter.GAP_VAL)) {
										if (reverseMode) tmpDD += modpair.vGammaCloses[mm-(ii-1)][nn-jj] ;
										else tmpDD += modpair.vGammaOpens[ii][jj] ;
									}
									if ( (A.seqs[p][i-2] == SequenceConverter.GAP_VAL && B.seqs[q][j-2] != SequenceConverter.GAP_VAL)
										|| (jj>0 && isPessimistic && A.seqs[p][i-2] == SequenceConverter.GAP_VAL && B.seqs[q][j-2] == SequenceConverter.GAP_VAL)) {
										if (reverseMode) tmpDD += modpair.hGammaOpens[mm-(ii-1)][nn-(jj-1)];
										else tmpDD += modpair.hGammaCloses[ii-1][jj];
									}
						    } else if (A.seqs[p][i-1] == SequenceConverter.GAP_VAL && B.seqs[q][j-1] != SequenceConverter.GAP_VAL) {
						    	if ( A.seqs[p][i-2] != SequenceConverter.GAP_VAL || 
										(isPessimistic && B.seqs[q][j-2] == SequenceConverter.GAP_VAL)) {
						    		if (reverseMode) tmpDD += modpair.hGammaCloses[mm-ii][nn-(jj-1)] ;
						    		else tmpDD += modpair.hGammaOpens[ii][jj] ;
								}
								if ( (A.seqs[p][i-2] != SequenceConverter.GAP_VAL && B.seqs[q][j-2] == SequenceConverter.GAP_VAL)
									|| (ii>0 && isPessimistic && A.seqs[p][i-2] == SequenceConverter.GAP_VAL && B.seqs[q][j-2] == SequenceConverter.GAP_VAL)) {
									if (reverseMode) tmpDD += modpair.vGammaOpens[mm-(ii-1)][nn-(jj-1)];
									else tmpDD += modpair.vGammaCloses[ii][jj-1];
								}
						    } else if (A.seqs[p][i-1] != SequenceConverter.GAP_VAL && B.seqs[q][j-1] != SequenceConverter.GAP_VAL) {
						    	if (A.seqs[p][i-2] != SequenceConverter.GAP_VAL && B.seqs[q][j-2] == SequenceConverter.GAP_VAL) { 
						    		if (reverseMode) tmpDD +=  modpair.vGammaOpens[mm-(ii-1-1)][nn-(jj-1)];
						    		else tmpDD +=  modpair.vGammaCloses[ii-1][jj-1];
						    	} else if (A.seqs[p][i-2] == SequenceConverter.GAP_VAL && B.seqs[q][j-2] != SequenceConverter.GAP_VAL) {
						    		if (reverseMode) tmpDD +=  modpair.hGammaOpens[mm-(ii-1)][nn-(jj-1-1)];
						    		else tmpDD +=  modpair.hGammaCloses[ii-1][jj-1];
						    	} else if (isPessimistic && (ii>1 || jj>1) && A.seqs[p][i-2] == SequenceConverter.GAP_VAL && B.seqs[q][j-2] == SequenceConverter.GAP_VAL) {
									if (ii>1 && jj>1) {
										if (reverseMode) tmpDD += Math.max(modpair.vGammaOpens[mm-(ii-1-1)][nn-(jj-1)], modpair.hGammaOpens[mm-(ii-1)][nn-(jj-1-1)]);
										else tmpDD += Math.max(modpair.vGammaCloses[ii-1][jj-1], modpair.hGammaCloses[ii-1][jj-1]);
									} else if (ii==1 && jj>1){ 
										if (reverseMode) tmpDD += modpair.hGammaOpens[mm-(ii-1)][nn-(jj-1-1)];
										else tmpDD += modpair.hGammaCloses[ii-1][jj-1];
									}else { // ii>1, jj==1 
										if (reverseMode) tmpDD += modpair.vGammaOpens[mm-(ii-1-1)][nn-(jj-1)];
										else tmpDD += modpair.vGammaCloses[ii-1][jj-1];
									}
						    	}
						    }
						}						
						if (j>1) { // H(0,j) isn't valid
							if (A.seqs[p][i-1] != SequenceConverter.GAP_VAL && B.seqs[q][j-1] == SequenceConverter.GAP_VAL) {
								if (B.seqs[q][j-2] != SequenceConverter.GAP_VAL || isPessimistic) { 
									if (reverseMode) tmpHD += modpair.vGammaCloses[mm-(ii-1)][nn-jj];
									else tmpHD += modpair.vGammaOpens[ii][jj];
									if (jj>1) {
										if (reverseMode) tmpHD += modpair.hGammaOpens[mm-(ii-1)][nn-(jj-1)];
										else tmpHD += modpair.hGammaCloses[ii-1][jj];
									}
								}
						    } else if (A.seqs[p][i-1] == SequenceConverter.GAP_VAL && B.seqs[q][j-1] != SequenceConverter.GAP_VAL
						    		&& isPessimistic && B.seqs[q][j-2] == SequenceConverter.GAP_VAL) {
						    		if (reverseMode) tmpHD += modpair.hGammaCloses[mm-ii][nn-(jj-1)];
						    		else tmpHD += modpair.hGammaOpens[ii][jj];
							    	if (ii>1) {
							    		if (reverseMode) tmpHD += modpair.vGammaOpens[mm-(ii-1)][nn-(jj-1)];
							    		else tmpHD += modpair.vGammaCloses[ii][jj-1];
							    	}
							} else if (A.seqs[p][i-1] != SequenceConverter.GAP_VAL && B.seqs[q][j-1] != SequenceConverter.GAP_VAL) {
								if (B.seqs[q][j-2] != SequenceConverter.GAP_VAL) {
									if (reverseMode) tmpHD +=  modpair.hGammaOpens[mm-(ii-1)][nn-(jj-1-1)];
									else tmpHD +=  modpair.hGammaCloses[ii-1][jj-1];
								} else if (isPessimistic && (ii>1 || jj>1) ) {
									if (ii>1 && jj>1) {
										if (reverseMode) tmpHD += Math.max(modpair.vGammaOpens[mm-(ii-1-1)][nn-(jj-1)], modpair.hGammaOpens[mm-(ii-1)][nn-(jj-1-1)]);
										else tmpHD += Math.max(modpair.vGammaCloses[ii-1][jj-1], modpair.hGammaCloses[ii-1][jj-1]);
									} else if (ii==1 && jj>1) { 
										
										if (reverseMode) tmpHD += modpair.hGammaOpens[mm-(ii-1)][nn-(jj-1-1)];
										else tmpHD += modpair.hGammaCloses[ii-1][jj-1];
									} else { // ii>1, jj==1 
										if (reverseMode) tmpHD += modpair.vGammaOpens[mm-(ii-1-1)][nn-(jj-1)];
										else tmpHD += modpair.vGammaCloses[ii-1][jj-1];
									}
								}
							}
						}

						if (i==1 && j==1) { // special case for beginning of matrix
							int tmp = 0;
							if (A.seqs[p][i-1] != SequenceConverter.GAP_VAL && B.seqs[q][j-1] == SequenceConverter.GAP_VAL) {
								//opening a gap
								if (reverseMode) tmp += modpair.vGammaCloses[mm-(ii-1)][nn-jj] ;
								else tmp += modpair.vGammaOpens[ii][jj] ;
							} else if (A.seqs[p][i-1] == SequenceConverter.GAP_VAL && B.seqs[q][j-1] != SequenceConverter.GAP_VAL) {
								if (reverseMode) tmp += modpair.hGammaCloses[mm-ii][nn-(jj-1)] ;
					    		else tmp += modpair.hGammaOpens[ii][jj] ;
							}
							tmpDD += tmp;
							tmpVD += tmp;
							tmpHD += tmp;
						} 
						
						if (i==M && j==N){
							if (A.seqs[p][M-1] != SequenceConverter.GAP_VAL && B.seqs[q][N-1] == SequenceConverter.GAP_VAL) { 
								if (reverseMode) tmpD += modpair.vGammaOpens[1][0];
								else tmpD += modpair.vGammaCloses[mm][nn];
							} else if (A.seqs[p][M-1] == SequenceConverter.GAP_VAL && B.seqs[q][N-1] != SequenceConverter.GAP_VAL) { 
								if (reverseMode) tmpD += modpair.hGammaOpens[0][1];
								else tmpD += modpair.hGammaCloses[mm][nn];
							} else if (isPessimistic && A.seqs[p][M-1] == SequenceConverter.GAP_VAL && B.seqs[q][N-1] == SequenceConverter.GAP_VAL) {
								if (reverseMode) tmpD += Math.max( modpair.hGammaOpens[0][1], modpair.vGammaOpens[1][0]);
								else tmpD += Math.max( modpair.hGammaCloses[mm][nn], modpair.vGammaCloses[mm][nn]);
							}
						}
					}
				}
				
				
				
				if ( cost == H[i-1][j-1] + tmpD + tmpHD )	{			
					nextdir = Direction.horiz;
					if (testing) LogWriter.stdOutLogln("hd(" + i + "," + j+ "): " + (tmpD+tmpHD));
				} else if ( cost == V[i-1][j-1] + tmpD + tmpVD) {				
					nextdir = Direction.vert;		
					if (testing) LogWriter.stdOutLogln("vd(" + i + "," + j+ "): " + (tmpD+tmpVD));
				} else if ( cost == D[i-1][j-1] + tmpD + tmpDD) {									 
					nextdir = Direction.diag;
					if (testing) LogWriter.stdOutLogln("dd(" + i + "," + j+ "): " + (tmpD+tmpDD));
				} else {
					LogWriter.stdErrLogln("no cost source found in dir == Direction.diag");
					return false;					
				}
				i--;
				j--;
			} else {
				LogWriter.stdErrLogln("Encountered an unknown direction in the DP table");
				return false;
			}
			
			path.add(dir);
			dir = nextdir;
			cost = (Direction.horiz == dir ) ? H[i][j] : (Direction.vert == dir ) ?  V[i][j] : D[i][j] ;
			
		}
		
		//I'm on a boundary. Just add the remaining slots
		while (i>0) { //then j==0
			path.add(Direction.vert);
			if (testing) LogWriter.stdOutLogln("v: (cost not calculated)");
			i--;
		}
		while (j>0) { //then i==0  ... only one of these ever happens
			path.add(Direction.horiz);
			if (testing) LogWriter.stdOutLogln("h: (cost not calculated)");
			j--;
		}

		Collections.reverse(path);
		return true;

	}
	
	public long[][] getD (){
		return D;
	}
	public long[][] getH (){
		return H;
	}
	public long[][] getV (){
		return V;
	}
	public void setReverseMode (boolean mode) {
		reverseMode = mode;
	}
}

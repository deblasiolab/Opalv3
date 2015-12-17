package opal.align;

import java.util.ArrayList;
import java.util.Collections;

import com.traviswheeler.libs.LogWriter;

public class ProfileAligner extends Aligner {
		
//	boolean testing = true;
	boolean testing = false;
	
	public static boolean doReverse = false;
	
	//vars in the profile subclasses
	public ProfileAligner ( Aligner al) {
		super(al);
	}
	public ProfileAligner() {
		this(true);
	}
	
	public ProfileAligner( boolean pess) {
		super(pess);
	}

	public ProfileAligner( Alignment A, Alignment B) {
		this(A,B,true);
	}

	public ProfileAligner( Alignment A, Alignment B, boolean pess) {
		super(A,B,pess);
	}
	
	public void align () {

		int[] idsAB = new int[A.seqIds.length + B.seqIds.length];
		System.arraycopy(A.seqIds, 0, idsAB, 0, A.seqIds.length);
		System.arraycopy(B.seqIds, 0, idsAB, A.seqIds.length, B.seqIds.length);

		
		int[][] result2 = null;			
		long cost2 = Long.MAX_VALUE;
		if (isPessimistic && doReverse) { 
			// I also want to do a reverse pass, to see if that gives a better alignment. 
			// If so, use that one instead
			Alignment alA = A; 
			Alignment alB = B; 
			int[][] revA = config.sc.buildReverseAlignment(A.seqs);
			int[][] revB = config.sc.buildReverseAlignment(B.seqs);

			A = Alignment.buildNewAlignment(revA, A.seqIds, true, config);
			B = Alignment.buildNewAlignment(revB, B.seqIds, true, config);
	
			super.align();
			
			result2 = config.sc.convertPathToIntAlignment (getPath(), A, B);

			if (config.useStructure) {
				// structure modifiers are stored forwards. Easiest code to write involves reversing this alignment
				int[][] rev = config.sc.buildReverseAlignment(result2);
				result2 = rev;		
			}
			
			cost2  = Aligner.calcCost(result2, A.K, B.K, idsAB, config);
		
			A = alA;
			B = alB;			
		}

		super.align(); // forward pass
		int[][] result1 = config.sc.convertPathToIntAlignment (getPath(), A, B);
		
		long cost1 = 0;
		if (isPessimistic && doReverse) 
			cost1 = Aligner.calcCost(result1, A.K, B.K, idsAB, config);
		
		
		int[] ids = new int[A.seqIds.length + B.seqIds.length];
		for (int i=0; i<A.seqIds.length; i++) ids[i] = A.seqIds[i];
		for (int i=0; i<B.seqIds.length; i++) ids[A.seqIds.length+i] = B.seqIds[i];
			
		if (cost2 < cost1) {
			if (config.useStructure) {
				resultAlignment = Alignment.buildNewAlignment(result2, ids, config);
			} else {
				resultAlignment = Alignment.buildNewAlignment(config.sc.buildReverseAlignment(result2), ids, config);
			}
		} else { 
			resultAlignment = Alignment.buildNewAlignment(result1, ids, config);
		}
	}
	
	protected void initialize () {
		// +1 to shift over the characters to a 1-based counting method
		H = new long[M+1][N+1];
		V = new long[M+1][N+1];
		D = new long[M+1][N+1];				
	}		

	public void cleanup () {
		H = null;
		V = null;
		D = null;
	}
	
	protected void fillBoundary() {
		int bigNumber = Integer.MAX_VALUE/2 ; // half, to ensure no overflow from the next row/column 

		V[0][0] = H[0][0] = D[0][1] = V[0][1] = D[1][0] = H[1][0] = bigNumber;
		H[0][1] = config.leftLambdaTerm()*K*B.f1[1] + config.leftGammaTerm()*K*B.firstLetterCount[1];
		V[1][0] = config.leftLambdaTerm()*L*A.f1[1] + config.leftGammaTerm()*L*A.firstLetterCount[1];
		D[0][0] = 0; 
		for (int j=2; j<=N; j++) {
			H[0][j] = H[0][j-1] + config.leftLambdaTerm()*K*B.f1[j] + config.leftGammaTerm()*K*B.firstLetterCount[j];
			D[0][j] = V[0][j] = bigNumber;
		}
		for (int i=2; i<=M; i++) {
			V[i][0] = V[i-1][0] + config.leftLambdaTerm()*L*A.f1[i] + config.leftGammaTerm()*L*A.firstLetterCount[i];
			D[i][0] = H[i][0] = bigNumber;
		}
	}
	
	protected void fillTable() {
		
		int begGapOpenCount_H, termGapExtCountLeft_H, termGapExtCountRight_H; 
		int begGapOpenCount_V, termGapExtCountLeft_V, termGapExtCountRight_V;
		int begGapOpenCount_D, termGapExtCountLeft_D, termGapExtCountRight_D;
		int fA1fB1;
		long x,y,z;
		for (int i=1; i<=M; i++){
			for (int j=1; j<=N; j++){				
				
				begGapOpenCount_H = A.gapsBeforeFirst[i] * B.firstLetterCount[j];				
				termGapExtCountLeft_H = B.f1[j] * (A.gapsBeforeFirst[i]) ;
				termGapExtCountRight_H = B.f1[j] * (A.gapsAfterLast[i]+A.lastLetterCount[i]) ;

				begGapOpenCount_V = B.gapsBeforeFirst[j] * A.firstLetterCount[i];
				termGapExtCountLeft_V = A.f1[i] * (B.gapsBeforeFirst[j]) ;				
				termGapExtCountRight_V = A.f1[i] * (B.gapsAfterLast[j]+B.lastLetterCount[j]) ;				

				begGapOpenCount_D = begGapOpenCount_H + begGapOpenCount_V;
				termGapExtCountLeft_D = B.f1[j] * (A.gapsBeforeFirst[i]) 
						+A.f1[i] * (B.gapsBeforeFirst[j]); 
				termGapExtCountRight_D = B.f1[j] * (A.gapsAfterLast[i]) 
						+A.f1[i] * (B.gapsAfterLast[j]); 
				
				fA1fB1 = A.f1[i] * B.f1[j];
				
				if (isPessimistic) {
					int endGapOpenCount_H_pess = B.f01[j] * (A.gapsAfterLast[i]+A.lastLetterCount[i]) ;
					int endGapOpenCount_Hv_pess = B.f1[j] * (A.gapsAfterLast[i]+A.lastLetterCount[i]);
					int endGapOpenCount_Hd_pess = B.f1[j] * A.lastLetterCount[i] +
												B.f01[j] * A.gapsAfterLast[i];
					
					int endGapOpenCount_V_pess = A.f01[i] * (B.gapsAfterLast[j]+B.lastLetterCount[j]) ;
					int endGapOpenCount_Vh_pess = A.f1[i] * (B.gapsAfterLast[j]+B.lastLetterCount[j]);
					int endGapOpenCount_Vd_pess = A.f1[i] * B.lastLetterCount[j] +
												  A.f01[i] * B.gapsAfterLast[j];

					int endGapOpenCount_Dv_pess = B.f1[j] * A.gapsAfterLast[i] + A.f01[i] * B.gapsAfterLast[j];  
					int endGapOpenCount_Dh_pess = A.f1[i] * B.gapsAfterLast[j] + B.f01[j] * A.gapsAfterLast[i];  
					int endGapOpenCount_Dd_pess = B.f1[j] * A.lastLetterCount[i-1] + B.f01[j] * (A.gapsAfterLast[i]-A.lastLetterCount[i-1])
												+ A.f1[i] * B.lastLetterCount[j-1] + A.f01[i] * (B.gapsAfterLast[j]-B.lastLetterCount[j-1]);
				
					
					x = H[i][j-1] + config.gamma*( (K * B.f01[j])- (begGapOpenCount_H+endGapOpenCount_H_pess) ) 
							+ config.leftGammaTerm()*(begGapOpenCount_H)
							+ config.rightGammaTerm()*(endGapOpenCount_H_pess);
					y = V[i][j-1] + config.gamma*( (K * B.f1[j]) - (begGapOpenCount_H+endGapOpenCount_Hv_pess) ) 
							+ config.leftGammaTerm()*(begGapOpenCount_H)
							+ config.rightGammaTerm()*(endGapOpenCount_Hv_pess);
					z = D[i][j-1] + config.gamma*( (fA1fB1 + A.f0[i]*B.f01[j]) -(begGapOpenCount_H+endGapOpenCount_Hd_pess) ) 
							+ config.leftGammaTerm()*(begGapOpenCount_H)
							+ config.rightGammaTerm()*(endGapOpenCount_Hd_pess);
					if (config.useStructure) {
						x += getStructGapOpenModifer((StructureAlignment)A, (StructureAlignment)B, i, j, Direction.horiz, Direction.horiz, isPessimistic); 
						y += getStructGapOpenModifer((StructureAlignment)A, (StructureAlignment)B, i, j, Direction.horiz, Direction.vert, isPessimistic); 
						z += getStructGapOpenModifer((StructureAlignment)A, (StructureAlignment)B, i, j, Direction.horiz, Direction.diag, isPessimistic); 
					}
					H[i][j] = x<y ? x<z?x:z  :  y<z?y:z  ;

					x = H[i-1][j] + config.gamma*( (L * A.f1[i])-(begGapOpenCount_V+endGapOpenCount_Vh_pess) ) 
							+ config.leftGammaTerm()*(begGapOpenCount_V)
							+ config.rightGammaTerm()*(endGapOpenCount_Vh_pess);
					y = V[i-1][j] + config.gamma*( (L * A.f01[i])-(begGapOpenCount_V+endGapOpenCount_V_pess) ) 
							+ config.leftGammaTerm()*(begGapOpenCount_V)
							+ config.rightGammaTerm()*(endGapOpenCount_V_pess);
					z = D[i-1][j] + config.gamma*( (fA1fB1 + A.f01[i]*B.f0[j])-(begGapOpenCount_V+endGapOpenCount_Vd_pess)) 
							+ config.leftGammaTerm()*(begGapOpenCount_V)
							+ config.rightGammaTerm()*(endGapOpenCount_Vd_pess);
					if (config.useStructure) {
						x += getStructGapOpenModifer((StructureAlignment)A, (StructureAlignment)B, i, j, Direction.vert, Direction.horiz, isPessimistic); 
						y += getStructGapOpenModifer((StructureAlignment)A, (StructureAlignment)B, i, j, Direction.vert, Direction.vert, isPessimistic); 
						z += getStructGapOpenModifer((StructureAlignment)A, (StructureAlignment)B, i, j, Direction.vert, Direction.diag, isPessimistic); 
					}					
					V[i][j] = x<y ? x<z?x:z  :  y<z?y:z  ;

					x = H[i-1][j-1] + config.gamma*( (A.f1[i]*B.f0[j] + A.f0[i]*B.f01[j])-(begGapOpenCount_D+endGapOpenCount_Dh_pess)) 
							+ config.leftGammaTerm()*(begGapOpenCount_D)
							+ config.rightGammaTerm()*(endGapOpenCount_Dh_pess);
					y = V[i-1][j-1] + config.gamma*( (B.f1[j]*A.f0[i] + B.f0[j]*A.f01[i])-(begGapOpenCount_D+endGapOpenCount_Dv_pess)) 
							+ config.leftGammaTerm()*(begGapOpenCount_D)
							+ config.rightGammaTerm()*(endGapOpenCount_Dv_pess);
					z = D[i-1][j-1] + config.gamma*( (A.f1[i]*B.f10[j] + A.f10[i]*B.f1[j] + A.f00[i]*B.f01[j] + A.f01[i]*B.f00[j] ) - (begGapOpenCount_D+endGapOpenCount_Dd_pess))
						+ config.leftGammaTerm()*(begGapOpenCount_D)
						+ config.rightGammaTerm()*(endGapOpenCount_Dd_pess);
					if (config.useStructure) {
						x += getStructGapOpenModifer((StructureAlignment)A, (StructureAlignment)B, i, j, Direction.diag, Direction.horiz, isPessimistic); 
						y += getStructGapOpenModifer((StructureAlignment)A, (StructureAlignment)B, i, j, Direction.diag, Direction.vert, isPessimistic); 
						z += getStructGapOpenModifer((StructureAlignment)A, (StructureAlignment)B, i, j, Direction.diag, Direction.diag, isPessimistic); 
					}
					D[i][j] = x<y ? x<z?x:z  :  y<z?y:z  ;
					
					
				} else {
						

					int endGapOpenCount_H_opt = B.f1[j] * A.lastLetterCount[i] ; // only for V and D
					int endGapOpenCount_V_opt = A.f1[i] * B.lastLetterCount[j] ; // only for H and D
					
					int endGapOpenCount_Dh_opt = A.f1[i] * B.lastLetterCount[j-1];
					int endGapOpenCount_Dv_opt = B.f1[j] * A.lastLetterCount[i-1];
					int endGapOpenCount_Dd_opt = endGapOpenCount_Dh_opt + endGapOpenCount_Dv_opt;
					
					//note: the optimistic counts don't include the case found by begGapOpenCounts,
					//  ... so no need to subtract them out.  (in other words, this tightens the optimistic bound a bit)

					x = H[i][j-1] + config.leftGammaTerm()*begGapOpenCount_H;
					y = V[i][j-1] + config.gamma * (fA1fB1 - (begGapOpenCount_H+endGapOpenCount_H_opt)) 
							+ config.leftGammaTerm()*(begGapOpenCount_H)
							+ config.rightGammaTerm()*(endGapOpenCount_H_opt);
					z = D[i][j-1] + config.gamma * (fA1fB1 - (begGapOpenCount_H+endGapOpenCount_H_opt)) 
							+ config.leftGammaTerm()*(begGapOpenCount_H)
							+ config.rightGammaTerm()*(endGapOpenCount_H_opt);
					if (config.useStructure) {
						x += getStructGapOpenModifer((StructureAlignment)A, (StructureAlignment)B, i, j, Direction.horiz, Direction.horiz, isPessimistic); 
						y += getStructGapOpenModifer((StructureAlignment)A, (StructureAlignment)B, i, j, Direction.horiz, Direction.vert, isPessimistic); 
						z += getStructGapOpenModifer((StructureAlignment)A, (StructureAlignment)B, i, j, Direction.horiz, Direction.diag, isPessimistic); 
					}
					H[i][j] = x<y ? x<z?x:z  :  y<z?y:z  ;

					x = H[i-1][j] + config.gamma * (fA1fB1 - (begGapOpenCount_V+endGapOpenCount_V_opt)) 
							+ config.leftGammaTerm()*(begGapOpenCount_V)
							+ config.rightGammaTerm()*(endGapOpenCount_V_opt);
					y = V[i-1][j] + config.leftGammaTerm()*begGapOpenCount_V;
					z = D[i-1][j] + config.gamma * (fA1fB1 - (begGapOpenCount_V+endGapOpenCount_V_opt)) 
							+ config.leftGammaTerm()*(begGapOpenCount_V)
							+ config.rightGammaTerm()*(endGapOpenCount_V_opt);
					if (config.useStructure) {
						x += getStructGapOpenModifer((StructureAlignment)A, (StructureAlignment)B, i, j, Direction.vert, Direction.horiz, isPessimistic); 
						y += getStructGapOpenModifer((StructureAlignment)A, (StructureAlignment)B, i, j, Direction.vert, Direction.vert, isPessimistic); 
						z += getStructGapOpenModifer((StructureAlignment)A, (StructureAlignment)B, i, j, Direction.vert, Direction.diag, isPessimistic); 
					}					
					V[i][j] = x<y ? x<z?x:z  :  y<z?y:z  ;

					int vTermGapOpens_notFirstCol = A.firstLetterCount[i]*(B.gapsBeforeFirst[j]-B.f10[j]);
					int hTermGapOpens_notFirstCol = B.firstLetterCount[j]*(A.gapsBeforeFirst[i]-A.f10[i]);
					
					x = H[i-1][j-1] + config.gamma*(A.f1[i]*B.f10[j] + vTermGapOpens_notFirstCol - (begGapOpenCount_V + endGapOpenCount_Dh_opt)) 
							+ config.leftGammaTerm()*(begGapOpenCount_D)
							+ config.rightGammaTerm()*(endGapOpenCount_Dh_opt);
					y = V[i-1][j-1] + config.gamma*(B.f1[j]*A.f10[i] + hTermGapOpens_notFirstCol - (begGapOpenCount_H + endGapOpenCount_Dv_opt)) 
							+ config.leftGammaTerm()*(begGapOpenCount_D)
							+ config.rightGammaTerm()*(endGapOpenCount_Dv_opt);
					z = D[i-1][j-1] + config.gamma*(A.f1[i]*B.f10[j] + A.f10[i]*B.f1[j] + vTermGapOpens_notFirstCol + hTermGapOpens_notFirstCol- (begGapOpenCount_D + endGapOpenCount_Dd_opt)) 
							+ config.leftGammaTerm()*(begGapOpenCount_D)
							+ config.rightGammaTerm()*(endGapOpenCount_Dd_opt);
					if (config.useStructure) {
						x += getStructGapOpenModifer((StructureAlignment)A, (StructureAlignment)B, i, j, Direction.diag, Direction.horiz, isPessimistic); 
						y += getStructGapOpenModifer((StructureAlignment)A, (StructureAlignment)B, i, j, Direction.diag, Direction.vert, isPessimistic); 
						z += getStructGapOpenModifer((StructureAlignment)A, (StructureAlignment)B, i, j, Direction.diag, Direction.diag, isPessimistic); 
					}			
					D[i][j] = x<y ? x<z?x:z  :  y<z?y:z  ;

				
				}

				H[i][j] += config.lambda * (K*B.f1[j] - termGapExtCountLeft_H - termGapExtCountRight_H)
						+ config.rightLambdaTerm() * termGapExtCountRight_H
						+ config.leftLambdaTerm() * termGapExtCountLeft_H;
				V[i][j] += config.lambda * (L*A.f1[i] - termGapExtCountLeft_V - termGapExtCountRight_V)
						+ config.rightLambdaTerm() * termGapExtCountRight_V
						+ config.leftLambdaTerm() * termGapExtCountLeft_V;

//				all substitution letter combos
				for (int m=0; m<A.chars[i].length; m++){
					for (int n=0; n<B.chars[j].length; n++){
						int s = config.cost.costs[A.chars[i][m]][B.chars[j][n]]; 
						D[i][j] += A.freqs[i][m]*B.freqs[j][n]*s;
					}
				}
				
				// plus extension costs				
				D[i][j] += config.lambda * ( (A.f0[i] * B.f1[j] + A.f1[i] * B.f0[j])-termGapExtCountLeft_D-termGapExtCountRight_D) 
						+ config.leftLambdaTerm() * termGapExtCountLeft_D
						+ config.rightLambdaTerm() * termGapExtCountRight_D;
				
				if (config.useStructure) {
					D[i][j] += getStructSubModifier((StructureAlignment)A, (StructureAlignment)B, i,j, config) +
							getStructGapExtModifer((StructureAlignment)A, (StructureAlignment)B, i, j, Direction.diag);
					V[i][j] += getStructGapExtModifer((StructureAlignment)A, (StructureAlignment)B, i, j, Direction.vert);
					H[i][j] += getStructGapExtModifer((StructureAlignment)A, (StructureAlignment)B, i, j, Direction.horiz);
				}				
//			    LogWriter.stdOutLogln("[" + i + "," + j + "]: H = " + H[i][j] + " , V = " + V[i][j] + " , D = " + D[i][j]);
			}
		}
	}	
	
	protected boolean recover () {
		
		path = new ArrayList<Direction>(2*Math.max(M,N));
		int i=M, j=N;
		Direction dir, nextdir=null;
		long cost = min(H[i][j], V[i][j], D[i][j]);
		estimatedCost = cost;

		long x=0,y=0,z=0;
		
		dir =  (cost == H[i][j]) ?  Direction.horiz : (cost == V[i][j]) ?  Direction.vert :  Direction.diag;
		int fA1fB1;
		
		//the calculations below are a bit hard to follow. Go up to the fillTable function to see them with white space		
		while (i>0 && j>0) {			
			fA1fB1 = A.f1[i] * B.f1[j];
									
			if (Direction.horiz == dir){
				int begGapOpenCount_H = A.gapsBeforeFirst[i] * B.firstLetterCount[j];				
				int termGapExtCountLeft_H = B.f1[j] * (A.gapsBeforeFirst[i]) ;
				int termGapExtCountRight_H = B.f1[j] * (A.gapsAfterLast[i]+A.lastLetterCount[i]) ;
				int endGapOpenCount_H_pess = B.f01[j] * (A.gapsAfterLast[i]+A.lastLetterCount[i]) ;
				int endGapOpenCount_Hv_pess = B.f1[j] * (A.gapsAfterLast[i]+A.lastLetterCount[i]);
				int endGapOpenCount_Hd_pess = B.f1[j] * A.lastLetterCount[i] + B.f01[j] * A.gapsAfterLast[i];
				int endGapOpenCount_H_opt = B.f1[j] * A.lastLetterCount[i] ; // only for V and D
								
				int base = config.lambda * (K*B.f1[j] - termGapExtCountLeft_H - termGapExtCountRight_H) 
						+ config.leftLambdaTerm() * termGapExtCountLeft_H
						+ config.rightLambdaTerm() * termGapExtCountRight_H;
				if (config.useStructure) {
					base += getStructGapExtModifer((StructureAlignment)A, (StructureAlignment)B, i, j, Direction.horiz);
					x = getStructGapOpenModifer((StructureAlignment)A, (StructureAlignment)B, i, j, Direction.horiz, Direction.horiz, isPessimistic); 
					y = getStructGapOpenModifer((StructureAlignment)A, (StructureAlignment)B, i, j, Direction.horiz, Direction.vert, isPessimistic); 
					z = getStructGapOpenModifer((StructureAlignment)A, (StructureAlignment)B, i, j, Direction.horiz, Direction.diag, isPessimistic); 
				}
				
				if ( cost == H[i][j-1] + base + x + ( isPessimistic ?  
						config.gamma * ( (K * B.f01[j])-(begGapOpenCount_H+endGapOpenCount_H_pess)) 
						+ config.leftGammaTerm()*(begGapOpenCount_H) 
						+ config.rightGammaTerm() * (endGapOpenCount_H_pess) 
						: config.leftGammaTerm()*begGapOpenCount_H )) {				
					nextdir = Direction.horiz;
					if (testing) LogWriter.stdOutLogln("hh(" + i + "," + j+ "): " + (cost));
				} else if ( cost == V[i][j-1] + base + y + (isPessimistic ? 
						config.gamma * ( (K * B.f1[j])-(begGapOpenCount_H+endGapOpenCount_Hv_pess) ) 
						+ config.leftGammaTerm()*(begGapOpenCount_H)
						+ config.rightGammaTerm() * (endGapOpenCount_Hv_pess)
						: config.gamma * (fA1fB1 - (begGapOpenCount_H+endGapOpenCount_H_opt)) 
						+ config.leftGammaTerm()*(begGapOpenCount_H)
						+ config.rightGammaTerm()*(endGapOpenCount_H_opt) ) ) {
					nextdir = Direction.vert;
					if (testing) LogWriter.stdOutLogln("vh(" + i + "," + j+ "): " + (cost));
				} else if ( cost == D[i][j-1] + base + z + (isPessimistic ? 
						config.gamma * ( (A.f1[i]*B.f1[j] + A.f0[i]*B.f01[j])-(begGapOpenCount_H+endGapOpenCount_Hd_pess)) 
						+ config.leftGammaTerm()*(begGapOpenCount_H)
						+ config.rightGammaTerm()*(endGapOpenCount_Hd_pess) 
						: config.gamma * (fA1fB1 - (begGapOpenCount_H+endGapOpenCount_H_opt)) 
						+ config.leftGammaTerm()*(begGapOpenCount_H)
						+ config.rightGammaTerm()*(endGapOpenCount_H_opt)  )) {				
					nextdir = Direction.diag;
					if (testing) LogWriter.stdOutLogln("dh(" + i + "," + j+ "): " + (cost));
				} else {
					if (testing) LogWriter.stdErrLogln("no cost source found in dir == Direction.horiz");
					return false;					
				}
				j--;
			} else if (Direction.vert == dir) {
				int begGapOpenCount_V = B.gapsBeforeFirst[j] * A.firstLetterCount[i];
				int termGapExtCountLeft_V = A.f1[i] * (B.gapsBeforeFirst[j]) ;	
				int termGapExtCountRight_V = A.f1[i] * (B.gapsAfterLast[j]+B.lastLetterCount[j]) ;				
				int endGapOpenCount_V_pess = A.f01[i] * (B.gapsAfterLast[j]+B.lastLetterCount[j]) ;
				int endGapOpenCount_Vh_pess = A.f1[i] * (B.gapsAfterLast[j]+B.lastLetterCount[j]);
				int endGapOpenCount_Vd_pess = A.f1[i] * B.lastLetterCount[j] +  A.f01[i] * B.gapsAfterLast[j];
				int endGapOpenCount_V_opt = A.f1[i] * B.lastLetterCount[j] ; // only for H and D
				
				int base = config.lambda * (L*A.f1[i] - termGapExtCountLeft_V - termGapExtCountRight_V) 
						+ config.leftLambdaTerm() * termGapExtCountLeft_V
						+ config.rightLambdaTerm() * termGapExtCountRight_V;
				if (config.useStructure) {
					base += getStructGapExtModifer((StructureAlignment)A, (StructureAlignment)B, i, j, Direction.vert);
					x = getStructGapOpenModifer((StructureAlignment)A, (StructureAlignment)B, i, j, Direction.vert, Direction.horiz, isPessimistic); 
					y = getStructGapOpenModifer((StructureAlignment)A, (StructureAlignment)B, i, j, Direction.vert, Direction.vert, isPessimistic); 
					z = getStructGapOpenModifer((StructureAlignment)A, (StructureAlignment)B, i, j, Direction.vert, Direction.diag, isPessimistic); 
				}
				
				if ( cost == H[i-1][j] + base + x +  ( isPessimistic? 
						config.gamma*( (L * A.f1[i])-(begGapOpenCount_V+endGapOpenCount_Vh_pess) ) 
						+ config.leftGammaTerm()*(begGapOpenCount_V)
						+ config.rightGammaTerm()*(endGapOpenCount_Vh_pess) 
						: config.gamma * (fA1fB1 - (begGapOpenCount_V+endGapOpenCount_V_opt)) 
						+ config.leftGammaTerm()*(begGapOpenCount_V)
						+ config.rightGammaTerm()*(endGapOpenCount_V_opt) )) {
					nextdir = Direction.horiz;
					if (testing) LogWriter.stdOutLogln("hv(" + i + "," + j+ "): " + (cost));
				} else if ( cost == V[i-1][j] + base + y + ( isPessimistic ? 
						config.gamma*( (L * A.f01[i])-(begGapOpenCount_V+endGapOpenCount_V_pess) ) 
						+ config.leftGammaTerm()*(begGapOpenCount_V)
						+ config.rightGammaTerm()*(endGapOpenCount_V_pess) 
						: config.leftGammaTerm()*begGapOpenCount_V )) {				
					nextdir = Direction.vert;
					if (testing) LogWriter.stdOutLogln("vv(" + i + "," + j+ "): " + (cost));
			 	} else if ( cost == D[i-1][j] + base + z +  ( isPessimistic ? 
			 			config.gamma*( (fA1fB1 + A.f01[i]*B.f0[j])-(begGapOpenCount_V+endGapOpenCount_Vd_pess) ) 
			 			+ config.leftGammaTerm()*(begGapOpenCount_V)
			 			+ config.rightGammaTerm()*(endGapOpenCount_Vd_pess) 
			 			: config.gamma * (fA1fB1 - (begGapOpenCount_V+endGapOpenCount_V_opt)) 
			 			+ config.leftGammaTerm()*(begGapOpenCount_V)
			 			+ config.rightGammaTerm()*(endGapOpenCount_V_opt) )) {					
					nextdir = Direction.diag;
					if (testing) LogWriter.stdOutLogln("dv(" + i + "," + j+ "): " + (cost));
			 	} else {
			 		LogWriter.stdErrLogln("no cost source found in dir == Direction.vert");
					return false;					
				}
				i--;
			} else if (Direction.diag == dir ) {
				int begGapOpenCount_H = A.gapsBeforeFirst[i] * B.firstLetterCount[j];				
				int begGapOpenCount_V = B.gapsBeforeFirst[j] * A.firstLetterCount[i];
				int begGapOpenCount_D = begGapOpenCount_H + begGapOpenCount_V;
				int termGapExtCountLeft_D = B.f1[j] * (A.gapsBeforeFirst[i]) 
						+ A.f1[i] * (B.gapsBeforeFirst[j]); 

				int termGapExtCountRight_D = B.f1[j] * (A.gapsAfterLast[i]) 
									+ A.f1[i] * (B.gapsAfterLast[j]); 

				int endGapOpenCount_Dv_pess = B.f1[j] * A.gapsAfterLast[i] + A.f01[i] * B.gapsAfterLast[j];  
				int endGapOpenCount_Dh_pess = A.f1[i] * B.gapsAfterLast[j] + B.f01[j] * A.gapsAfterLast[i];  
				int endGapOpenCount_Dd_pess = B.f1[j] * A.lastLetterCount[i-1] + B.f01[j] * (A.gapsAfterLast[i]-A.lastLetterCount[i-1])
											+ A.f1[i] * B.lastLetterCount[j-1] + A.f01[i] * (B.gapsAfterLast[j]-B.lastLetterCount[j-1]);

				int endGapOpenCount_Dh_opt = A.f1[i] * B.lastLetterCount[j-1];
				int endGapOpenCount_Dv_opt = B.f1[j] * A.lastLetterCount[i-1];
				int endGapOpenCount_Dd_opt = endGapOpenCount_Dh_opt + endGapOpenCount_Dv_opt;

				int base = config.lambda * ( (A.f0[i] * B.f1[j] + A.f1[i] * B.f0[j])-termGapExtCountLeft_D-termGapExtCountRight_D) 
						+ config.leftLambdaTerm() * termGapExtCountLeft_D
						+ config.rightLambdaTerm() * termGapExtCountRight_D;

				if (config.useStructure) {
					base += getStructSubModifier((StructureAlignment)A, (StructureAlignment)B, i,j, config) +
						getStructGapExtModifer((StructureAlignment)A, (StructureAlignment)B, i, j, Direction.diag);
					x = getStructGapOpenModifer((StructureAlignment)A, (StructureAlignment)B, i, j, Direction.diag, Direction.horiz, isPessimistic); 
					y = getStructGapOpenModifer((StructureAlignment)A, (StructureAlignment)B, i, j, Direction.diag, Direction.vert, isPessimistic); 
					z = getStructGapOpenModifer((StructureAlignment)A, (StructureAlignment)B, i, j, Direction.diag, Direction.diag, isPessimistic); 
				}
								
				for (int p=0; p<A.chars[i].length; p++){
					for (int q=0; q<B.chars[j].length; q++){
						int s = config.cost.costs[A.chars[i][p]][B.chars[j][q]];
						base += A.freqs[i][p]*B.freqs[j][q]*s;
					}
				}			

				int vTermGapOpens_notFirstCol = A.firstLetterCount[i]*(B.gapsBeforeFirst[j]-B.f10[j]);
				int hTermGapOpens_notFirstCol = B.firstLetterCount[j]*(A.gapsBeforeFirst[i]-A.f10[i]);

				if ( cost == H[i-1][j-1] + base + x + ( isPessimistic ? 
						config.gamma * ((A.f1[i]*B.f0[j] + A.f0[i]*B.f01[j])-(begGapOpenCount_D+endGapOpenCount_Dh_pess)) 
						+ config.leftGammaTerm()*(begGapOpenCount_D)
						+ config.rightGammaTerm()*(endGapOpenCount_Dh_pess) 
						: config.gamma*(A.f1[i]*B.f10[j] + vTermGapOpens_notFirstCol - (begGapOpenCount_V+endGapOpenCount_Dh_opt)) 
						+ config.leftGammaTerm()*(begGapOpenCount_D) 
						+ config.rightGammaTerm()*(endGapOpenCount_Dh_opt) )) {				
					nextdir = Direction.horiz;
					if (testing) LogWriter.stdOutLogln("hd(" + i + "," + j+ "): " + (cost));
				} else if ( cost == V[i-1][j-1] + base + y + ( isPessimistic ? 
						config.gamma * ( (B.f1[j]*A.f0[i] + B.f0[j]*A.f01[i])-(begGapOpenCount_D+endGapOpenCount_Dv_pess)) 
						+ config.leftGammaTerm()*(begGapOpenCount_D)
						+ config.rightGammaTerm()*(endGapOpenCount_Dv_pess) 
						: config.gamma*(B.f1[j]*A.f10[i] + hTermGapOpens_notFirstCol - (begGapOpenCount_H+endGapOpenCount_Dv_opt)) 
						+ config.leftGammaTerm()*(begGapOpenCount_D)
						+ config.rightGammaTerm()*(endGapOpenCount_Dv_opt) )) {				
					nextdir = Direction.vert;		
					if (testing) LogWriter.stdOutLogln("vd(" + i + "," + j+ "): " + (cost));
				} else if ( cost == D[i-1][j-1] + base + z + ( isPessimistic ? 
						config.gamma * ( (A.f1[i]*B.f10[j] + A.f10[i]*B.f1[j] + A.f00[i]*B.f01[j] + A.f01[i]*B.f00[j])-(begGapOpenCount_D+endGapOpenCount_Dd_pess)) 
						+ config.leftGammaTerm()*(begGapOpenCount_D)
						+ config.rightGammaTerm()*(endGapOpenCount_Dd_pess) 
						: config.gamma*( A.f1[i]*B.f10[j] + A.f10[i]*B.f1[j] + vTermGapOpens_notFirstCol + hTermGapOpens_notFirstCol- (begGapOpenCount_D+endGapOpenCount_Dd_opt)) 
						+ config.leftGammaTerm()*(begGapOpenCount_D) 
						+ config.rightGammaTerm()*(endGapOpenCount_Dd_opt) )) {									 
					nextdir = Direction.diag;
					if (testing) LogWriter.stdOutLogln("dd(" + i + "," + j+ "): " + (cost));
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
			if (testing){
				cost = V[i][j] ;
				LogWriter.stdOutLogln("v(" + i + "," + j+ "): " + (cost));
			}
			i--;
		}
		while (j>0) { //then i==0  ... only one of these ever happens
			path.add(Direction.horiz);
			if (testing){
				cost = H[i][j] ;
				LogWriter.stdOutLogln("h(" + i + "," + j+ "): " + (cost));
			}
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

	
}

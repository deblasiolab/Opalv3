package opal.align;

import java.util.ArrayList;
import java.util.Collections;

import com.traviswheeler.libs.LogWriter;


public class PairAligner extends Aligner {
		
	public PairAligner ( Aligner al) {
		super(al);
	}
	
	public PairAligner( LogWriter logger ) {
		super(true);
	}
	public PairAligner( Alignment A, Alignment B) {
		super(A,B,true);
	}
	
	public PairAligner( boolean pess) {
		super(pess);
	}
	
	public PairAligner( Alignment A, Alignment B, boolean pess) {
		super(A,B,pess);
	}
	

	protected void initialize () {
		
		//the +1 thing here shifts over the characters to a 1-based counting method		
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
		V[0][0] = H[0][0] = D[0][0] = 0; 
		long bigNumber = Long.MAX_VALUE/2; // half, to ensure no overflow from the next row/column 
		for (int j=1; j<=N; j++) {
			H[0][j] = H[0][j-1] + config.leftLambdaTerm() + (j==1?config.leftGammaTerm():0);
			D[0][j] = V[0][j] = bigNumber; //H[0][j] + gamma;
		//	LogWriter.stdOutLogln("V[0][" + j + "] = " + V[0][j] + ";  H[0][" + j + "] = " + H[0][j] + ";  D[0][" + j + "] = " + D[0][j]);
		}
		for (int i=1; i<=M; i++) {
			V[i][0] = V[i-1][0] + config.leftLambdaTerm() + (i==1?config.leftGammaTerm():0);
			D[i][0] = H[i][0] = bigNumber;//= V[i][0] + gamma ;
		//	LogWriter.stdOutLogln("V[" + i + "][0] = " + V[i][0] + ";  H[" + i + "][0] = " + H[i][0] + ";  D[" + i + "][0] = " + D[i][0]);
		}
//		ShutdownHandler.shutdown(1);
	}
	
	protected void fillTable() {
		int gamma_tmp;
		int lambda_tmp;
		for (int i=1; i<=M; i++){
			for (int j=1; j<=N; j++){
					if (i==M || j==N) {
						gamma_tmp = config.rightGammaTerm();
						lambda_tmp = config.rightLambdaTerm();
					} else {
						gamma_tmp = config.gamma;
						lambda_tmp = config.lambda;						
					}

					long structGamma_tmpA=0;
					long structGamma_tmpB=0;
					if (config.useStructure) {
						structGamma_tmpA = Math.round(getStructGapOpenModiferPair(A.seqIds[0], i-1, config));
						structGamma_tmpB = Math.round(getStructGapOpenModiferPair(B.seqIds[0], j-1, config));
					}

					
					H[i][j] = min ( H[i][j-1], 
					                V[i][j-1] + gamma_tmp + structGamma_tmpA,
					                D[i][j-1] + gamma_tmp + structGamma_tmpA
							);
					
					V[i][j] = min ( 
							H[i-1][j] + gamma_tmp + structGamma_tmpB,
							V[i-1][j], 
							D[i-1][j] + gamma_tmp + structGamma_tmpB
							);
					
					D[i][j] = min ( 
							H[i-1][j-1],
							V[i-1][j-1],
							D[i-1][j-1] 
					       );
				
				H[i][j] += lambda_tmp;
				V[i][j] += lambda_tmp; 

				//substitution cost
				D[i][j] += config.cost.costs[A.seqs[0][i-1]][B.seqs[0][j-1]];

				if (config.useStructure) {
					D[i][j] += Math.round(getStructSubModifierPair(A.seqIds[0], B.seqIds[0], i-1, j-1, config));
					V[i][j] += Math.round(getStructGapExtModiferPair(A.seqIds[0], i-1, config));
					H[i][j] += Math.round(getStructGapExtModiferPair(B.seqIds[0], j-1, config));
				}
			}
		}
	}
		
	protected boolean recover () {
		
//		HeatMap heatmap = new HeatMap(A.M,B.M,"AB.png");
		
		path = new ArrayList<Direction>(2*Math.max(M,N));
		int i=M, j=N;
		Direction dir, nextdir=null;
		long cost = min(H[i][j], V[i][j], D[i][j]);
		estimatedCost = cost;
		
		dir =  (cost == H[i][j]) ?  Direction.horiz : (cost == V[i][j]) ?  Direction.vert :  Direction.diag;
		
		int gamma_tmp=0, base=0;
		
		//the calculations below are a bit hard to follow. Go up to the fillTable function to see them with white space		
		while (i>0 && j>0) {			
		
//			heatmap.setPointHeat(j, i, 4);
			
			if (Direction.diag != dir) {
				if (i==M || j==N) {
					gamma_tmp = config.rightGammaTerm();
					base = config.rightLambdaTerm();
				} else {
					gamma_tmp = config.gamma;
					base = config.lambda;						
				}
			}

			long structGamma_tmpA=0;
			long structGamma_tmpB=0;
			if (config.useStructure) {
				structGamma_tmpA = Math.round(getStructGapOpenModiferPair(A.seqIds[0], i-1, config));
				structGamma_tmpB = Math.round(getStructGapOpenModiferPair(B.seqIds[0], j-1, config));
			}

			
			if (Direction.horiz == dir){
				if (config.useStructure) base += Math.round(getStructGapExtModiferPair(B.seqIds[0], j-1, config));								
				if ( cost == H[i][j-1] + base )				
					nextdir = Direction.horiz;
				else if ( cost == V[i][j-1] + base + gamma_tmp + structGamma_tmpA)
					nextdir = Direction.vert;
				else if ( cost == D[i][j-1] + base + gamma_tmp + structGamma_tmpA)				
					nextdir = Direction.diag;
				else {
					LogWriter.stdErrLogln("no cost source found in dir == Direction.horiz");
					return false;					
				}
				j--;
			} else if (Direction.vert == dir) {
				if (config.useStructure) base += Math.round(getStructGapExtModiferPair(A.seqIds[0], i-1, config));
				if ( cost == H[i-1][j] + base + gamma_tmp + structGamma_tmpB)
					nextdir = Direction.horiz;
				else if ( cost == V[i-1][j] + base )				
					nextdir = Direction.vert;
				else if ( cost == D[i-1][j] + base + gamma_tmp + structGamma_tmpB)					
					nextdir = Direction.diag;
				else {
					LogWriter.stdErrLogln("no cost source found in dir == Direction.vert");
					return false;					
				}
				i--;
			} else if (Direction.diag == dir ) {
				base = config.cost.costs[A.seqs[0][i-1]][B.seqs[0][j-1]];
				if (config.useStructure) base += Math.round(getStructSubModifierPair(A.seqIds[0], B.seqIds[0], i-1, j-1, config));

				if ( cost == H[i-1][j-1] + base )				
					nextdir = Direction.horiz;
				else if ( cost == V[i-1][j-1] + base )				
					nextdir = Direction.vert;		
				else if ( cost == D[i-1][j-1] + base )									 
					nextdir = Direction.diag;
				else {
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
//			heatmap.setPointHeat(j, i, 4);
			path.add(Direction.vert);
			i--;
		}
		while (j>0) { //then i==0  ... only one of these ever happens
//			heatmap.setPointHeat(j, i, 4);
			path.add(Direction.horiz);
			j--;
		}
		
//		heatmap.saveImage();
		
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

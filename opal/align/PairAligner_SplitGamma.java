package opal.align;

import java.util.ArrayList;
import java.util.Collections;
import opal.IO.Configuration;

import com.traviswheeler.libs.LogWriter;

public class PairAligner_SplitGamma extends Aligner {
		
	public PairAligner_SplitGamma ( Aligner al) {
		super(al);
	}
	
	public PairAligner_SplitGamma( ) {
		super(true);
	}
	public PairAligner_SplitGamma( Alignment A, Alignment B) {
		super(A,B,true);
	}
	
	public PairAligner_SplitGamma( boolean pess) {
		super(pess);
	}
	
	public PairAligner_SplitGamma( Alignment A, Alignment B, boolean pess) {
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
			H[0][j] = H[0][j-1] + config.leftLambdaTerm() + (j==1?config.leftGammaTerm()/2:0);
			D[0][j] = V[0][j] = bigNumber; //H[0][j] + gamma;
		//	LogWriter.stdOutLogln("V[0][" + j + "] = " + V[0][j] + ";  H[0][" + j + "] = " + H[0][j] + ";  D[0][" + j + "] = " + D[0][j]);
		}
		for (int i=1; i<=M; i++) {
			V[i][0] = V[i-1][0] + config.leftLambdaTerm() + (i==1?config.leftGammaTerm()/2:0);
			D[i][0] = H[i][0] = bigNumber;//= V[i][0] + gamma ;
		//	LogWriter.stdOutLogln("V[" + i + "][0] = " + V[i][0] + ";  H[" + i + "][0] = " + H[i][0] + ";  D[" + i + "][0] = " + D[i][0]);
		}
//		ShutdownHandler.shutdown(1);
	}

	protected void fillTable() {
		int gamma_h_open_tmp, gamma_v_open_tmp;
		int gamma_close_tmp;
		for (int i=1; i<=M; i++){
			for (int j=1; j<=N; j++){
					if (i==1 || j==1) {
						gamma_close_tmp = config.leftGammaTerm()/2;
					} else {
						gamma_close_tmp = config.gamma/2;
					}	
					gamma_h_open_tmp = (i==M?config.rightGammaTerm():config.gamma)/2;
					gamma_v_open_tmp = (j==N?config.rightGammaTerm():config.gamma)/2;
					
					H[i][j] = min ( H[i][j-1], 
					                V[i][j-1] + gamma_h_open_tmp + gamma_close_tmp,
					                D[i][j-1] + gamma_h_open_tmp
							);
					
					V[i][j] = min ( 
							H[i-1][j] + gamma_close_tmp + gamma_v_open_tmp,
							V[i-1][j], 
							D[i-1][j] + gamma_v_open_tmp
							);
					
					D[i][j] = min ( 
							H[i-1][j-1] + gamma_close_tmp,
							V[i-1][j-1] + gamma_close_tmp,
							D[i-1][j-1] 
					       );
				

				H[i][j] += (i==M ? config.rightLambdaTerm() : config.lambda); 
				V[i][j] += (j==N ? config.rightLambdaTerm() : config.lambda); 

				//substitution cost
				D[i][j] += config.cost.costs[A.seqs[0][i-1]][B.seqs[0][j-1]];
			}
		}
		//close terminal gap
		H[M][N] += config.rightGammaTerm()/2;
		V[M][N] += config.rightGammaTerm()/2;		
	}	

	
	protected boolean recover () {
		
		path = new ArrayList<Direction>(2*Math.max(M,N));
		int i=M, j=N;
		Direction dir, nextdir=null;
		long cost = min(H[i][j], V[i][j], D[i][j]);
		estimatedCost = cost;
		
		dir =  (cost == H[i][j]) ?  Direction.horiz : (cost == V[i][j]) ?  Direction.vert :  Direction.diag;
		
		//clean up the extra cost added at end of terminal gaps
		if (Direction.horiz == dir || Direction.vert == dir) {
				cost -= config.rightGammaTerm()/2;
		}
		
//		int base=0;
		int gamma_open_tmp=0;
		int gamma_close_tmp=0;
		
		while (i>0 && j>0) {					
		
			if (Direction.diag != dir) {
				if ( (Direction.horiz == dir && i==M) || (Direction.vert == dir && j==N) ) {
					gamma_open_tmp = config.rightGammaTerm()/2;
				} else {
					gamma_open_tmp = config.gamma/2;		
				}
			}
			if (i==1 || j==1) {
				gamma_close_tmp = config.leftGammaTerm()/2;
			} else {
				gamma_close_tmp = config.gamma/2;
			}
			
			if (Direction.horiz == dir){
				int ext = (i==M?config.rightLambdaTerm() : config.lambda);
				if ( cost == H[i][j-1] + ext )				
					nextdir = Direction.horiz;
				else if ( cost == V[i][j-1] + ext + gamma_open_tmp + gamma_close_tmp)
					nextdir = Direction.vert;
				else if ( cost == D[i][j-1] + ext + gamma_open_tmp)				
					nextdir = Direction.diag;
				else {
					LogWriter.stdErrLogln("no cost source found in dir == Direction.horiz.  Pos " + i + "," + j);
					return false;					
				}
				j--;
			} else if (Direction.vert == dir) {
				int ext = (j==N?config.rightLambdaTerm() : config.lambda);
				if ( cost == H[i-1][j] + ext + gamma_open_tmp + gamma_close_tmp)
					nextdir = Direction.horiz;
				else if ( cost == V[i-1][j] + ext )				
					nextdir = Direction.vert;
				else if ( cost == D[i-1][j] + ext + gamma_open_tmp)					
					nextdir = Direction.diag;
				else {
					LogWriter.stdErrLogln("no cost source found in dir == Direction.vert.  Pos " + i + "," + j);
					return false;					
				}
				i--;
			} else if (Direction.diag == dir ) {
				int sub = config.cost.costs[A.seqs[0][i-1]][B.seqs[0][j-1]];

				if ( cost == H[i-1][j-1] + sub + gamma_close_tmp)				
					nextdir = Direction.horiz;
				else if ( cost == V[i-1][j-1] + sub + gamma_close_tmp )				
					nextdir = Direction.vert;		
				else if ( cost == D[i-1][j-1] + sub )									 
					nextdir = Direction.diag;
				else {
					LogWriter.stdErrLogln("no cost source found in dir == Direction.diag.  Pos " + i + "," + j);
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
			i--;
		}
		while (j>0) { //then i==0  ... only one of these ever happens
			path.add(Direction.horiz);
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

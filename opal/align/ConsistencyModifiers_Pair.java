package opal.align;


public class ConsistencyModifiers_Pair {

	public int[][] vGammaOpens;
	public int[][] hGammaOpens;
	public int[][] vGammaCloses;
	public int[][] hGammaCloses;
	public int[][] vLambdas;
	public int[][] hLambdas;
	public int[][] subs;
	public int[][] vGammaOpensAB;
	public int[][] hGammaOpensAB;
	public int[][] vGammaClosesAB;
	public int[][] hGammaClosesAB;
	public int[][] vLambdasAB;
	public int[][] hLambdasAB;
	public int[][] subsAB;
	int M;
	int N;
	
	
	public ConsistencyModifiers_Pair (int M, int N) {
		this.M = M;
		this.N = N;
		vGammaOpens = new int[M+1][N+1];
		hGammaOpens = new int[M+1][N+1];
		vGammaCloses = new int[M+1][N+1];
		hGammaCloses = new int[M+1][N+1];
		vLambdas = new int[M+1][N+1];
		hLambdas = new int[M+1][N+1];
		subs = new int[M+1][N+1];
	
		vGammaOpensAB = new int[M+1][N+1];
		hGammaOpensAB = new int[M+1][N+1];
		vGammaClosesAB = new int[M+1][N+1];
		hGammaClosesAB = new int[M+1][N+1];
		vLambdasAB = new int[M+1][N+1];
		hLambdasAB = new int[M+1][N+1];
		subsAB = new int[M+1][N+1];
		
	}
	/*
	public ConsistencyModifiers_Pair getReverse() {
		ConsistencyModifiers_Pair revpair = new ConsistencyModifiers_Pair(M,N);
		for (int i=0; i<=M; i++) {
			for (int j=0; j<=N; j++) {
				revpair.vGammaOpens[i][j] = vGammaCloses[M-i][N-j];
				revpair.hGammaOpens[i][j] = hGammaCloses[M-i][N-j];
				revpair.vGammaCloses[i][j] = vGammaOpens[M-i][N-j];
				revpair.hGammaCloses[i][j] = hGammaOpens[M-i][N-j];
				revpair.vLambdas[i][j] = vLambdas[M-i][N-j];
				revpair.hLambdas[i][j] = hLambdas[M-i][N-j];
				revpair.subs[i][j] = subs[M-i][N-j];
			}
		}
		return revpair;
	}*/

}

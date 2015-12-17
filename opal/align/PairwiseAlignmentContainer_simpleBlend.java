package opal.align;


import opal.tree.Tree;
import opal.IO.Configuration;

public class PairwiseAlignmentContainer_simpleBlend extends
		PairwiseAlignmentsContainer {

	public PairwiseAlignmentContainer_simpleBlend(Tree tree, float[][] distances,Configuration c) {
		super(tree, distances,c);
	}

	protected void setBlendParams (int neighborCnt) {
//		c_weight = neighborCnt==0? 0 : consistency_weight/neighborCnt;
		c_weight = consistency_weight;
		normalizer = neighborCnt==0? 0 : (float)1/ (neighborCnt * PairSuboptimalityMatrices.delta);
	}
	
	protected int calcSub(int a, int b, int i, int j, ConsistencyModifiers_Pair modpair, int neighborCnt) {
		int sigma = config.cost.costs[origSeqs[a][i-1]][origSeqs[b][j-1]];
		float mod = normalizer * modpair.subs[i][j];
		return Math.round( sigma * (1 + c_weight * mod));
	}
	
	protected int calcVLambda(int i, int j, ConsistencyModifiers_Pair modpair, int neighborCnt) {
		int v_ext = (j==0||j==N?((j==0)?config.leftLambdaTerm():config.rightLambdaTerm()) : config.lambda);
		float mod = normalizer * modpair.vLambdas[i][j];
		return Math.round( v_ext * (1 + c_weight * mod));
	}
	
	protected int calcVGammaOpen(int i, int j, ConsistencyModifiers_Pair modpair, int neighborCnt){
		int v_open = ( (j==0&&i==1)||j==N ? ((j==N)?config.rightLambdaTerm():config.leftLambdaTerm()) : config.gamma)/2;
		float mod = normalizer * modpair.vGammaOpens[i][j];
		return Math.round( v_open * (1 + c_weight * mod));
	}
	
	protected int calcVGammaClose(int i, int j, ConsistencyModifiers_Pair modpair, int neighborCnt){
		int v_close = ( j==0||(j==N&&i==M) ? ((j==0)?config.leftLambdaTerm():config.rightLambdaTerm()) : config.gamma)/2;
		float mod = normalizer * modpair.vGammaCloses[i][j];
		return Math.round( v_close * (1 + c_weight * mod));
	}
	
	
	protected int calcHLambda(int i, int j, ConsistencyModifiers_Pair modpair, int neighborCnt){
		int h_ext = (i==0||i==M?((i==0)?config.leftLambdaTerm():config.rightLambdaTerm()) : config.lambda);
		float mod = normalizer * modpair.hLambdas[i][j];
		return Math.round( h_ext * (1 + c_weight * mod));
	}
	
	protected int calcHGammaOpen(int i, int j, ConsistencyModifiers_Pair modpair, int neighborCnt){
		int h_open = ( (i==0&&j==1)||i==M ? ((i==M)?config.rightLambdaTerm():config.leftLambdaTerm()) : config.gamma)/2;
		float mod = normalizer * modpair.hGammaOpens[i][j];
		return Math.round( h_open * (1 + c_weight * mod));

	}
	
	protected int calcHGammaClose(int i, int j, ConsistencyModifiers_Pair modpair, int neighborCnt){
		int h_close = ( i==0||(j==N&&i==M) ? ((i==0)?config.leftLambdaTerm():config.rightLambdaTerm()): config.gamma)/2;
		float mod = normalizer * modpair.hGammaCloses[i][j];
		return Math.round( h_close * (1 + c_weight * mod));
	}
	
	
}

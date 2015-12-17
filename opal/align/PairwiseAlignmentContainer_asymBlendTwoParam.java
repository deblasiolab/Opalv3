package opal.align;


import opal.tree.Tree;
import opal.IO.Configuration;

public class PairwiseAlignmentContainer_asymBlendTwoParam extends
		PairwiseAlignmentsContainer {

	float alpha;
	
	public PairwiseAlignmentContainer_asymBlendTwoParam(Tree tree, float[][] distances, Configuration c) {
		super(tree, distances, c);
	}


	protected final void calcPrework () {
		alpha = consistency_other_seqs_weight / (1+ consistency_other_seqs_weight);
		normalizer = PairSuboptimalityMatrices.delta * 1; //Aligner.gammas.length;
	}
	
	protected final void setBlendParams (int neighborCnt) {
	}

	protected int calcSub(int a, int b, int i, int j, ConsistencyModifiers_Pair modpair, int neighborCnt) {
		int sigma = config.cost.costs[origSeqs[a][i-1]][origSeqs[b][j-1]];
		float mod;
		if (neighborCnt==0)
			mod = modpair.subsAB[i][j]/normalizer;
		else
			mod = ((1-alpha) * modpair.subsAB[i][j] + alpha * modpair.subs[i][j]/neighborCnt)/normalizer;
		return Math.round( sigma * (1 + consistency_weight * mod));
	}
	
	protected int calcVLambda(int i, int j, ConsistencyModifiers_Pair modpair, int neighborCnt) {
		int v_ext = (j==0||j==N?((j==0)?config.leftLambdaTerm():config.rightLambdaTerm()) : config.lambda);
		float mod;
		if (neighborCnt==0)
			mod = modpair.vLambdasAB[i][j]/normalizer;
		else
			mod = ((1-alpha) * modpair.vLambdasAB[i][j] + alpha * modpair.vLambdas[i][j]/neighborCnt)/normalizer;
		return Math.round( v_ext * (1 + consistency_weight * mod));
	}
	
	protected int calcVGammaOpen(int i, int j, ConsistencyModifiers_Pair modpair, int neighborCnt){
		int v_open = ( (j==0&&i==1)||j==N ?((j==0)?config.leftLambdaTerm():config.rightLambdaTerm()) : config.gamma)/2;
		float mod;
		if (neighborCnt==0)
			mod = modpair.vGammaOpensAB[i][j]/normalizer;
		else
			mod = ((1-alpha) * modpair.vGammaOpensAB[i][j] + alpha * modpair.vGammaOpens[i][j]/neighborCnt)/normalizer;
		return Math.round( v_open * (1 + consistency_weight * mod));

	}
	
	protected int calcVGammaClose(int i, int j, ConsistencyModifiers_Pair modpair, int neighborCnt){
		int v_close = ( j==0||(j==N&&i==M) ? ((j==0)?config.leftLambdaTerm():config.rightLambdaTerm()) : config.gamma)/2;
		float mod;
		if (neighborCnt==0)
			mod = modpair.vGammaClosesAB[i][j]/normalizer;
		else
			mod = ((1-alpha) * modpair.vGammaClosesAB[i][j] + alpha * modpair.vGammaCloses[i][j]/neighborCnt)/normalizer;
		return Math.round( v_close * (1 + consistency_weight * mod));
		
	}

	
	protected int calcHLambda(int i, int j, ConsistencyModifiers_Pair modpair, int neighborCnt){
		int h_ext = (i==0||i==M?((i==0)?config.leftLambdaTerm():config.rightLambdaTerm()) : config.lambda);
		float mod;
		if (neighborCnt==0)
			mod = modpair.hLambdasAB[i][j]/normalizer;
		else 
			mod = ((1-alpha) * modpair.hLambdasAB[i][j] + alpha * modpair.hLambdas[i][j]/neighborCnt)/normalizer;
		return Math.round( h_ext * (1 + consistency_weight * mod));
	}
	
	protected int calcHGammaOpen(int i, int j, ConsistencyModifiers_Pair modpair, int neighborCnt){
		int h_open = ( (i==0&&j==1)||i==M ? ((j==0)?config.leftLambdaTerm():config.rightLambdaTerm()) : config.gamma)/2;
		float mod;
		if (neighborCnt==0)
			mod = modpair.hGammaOpensAB[i][j]/normalizer;
		else 
			mod = ((1-alpha) * modpair.hGammaOpensAB[i][j] + alpha * modpair.hGammaOpens[i][j]/neighborCnt)/normalizer;
		return Math.round( h_open * (1 + consistency_weight * mod));
	}
	
	protected int calcHGammaClose(int i, int j, ConsistencyModifiers_Pair modpair, int neighborCnt){
		int h_close = ( i==0||(j==N&&i==M) ? ((i==0)?config.leftLambdaTerm():config.rightLambdaTerm()) : config.gamma)/2;
		float mod;
		if (neighborCnt==0)
			mod = modpair.hGammaClosesAB[i][j]/normalizer;
		else 
			mod = ((1-alpha) * modpair.hGammaClosesAB[i][j] + alpha * modpair.hGammaCloses[i][j]/neighborCnt)/normalizer;
		return Math.round( h_close * (1 + consistency_weight * mod));
	}
	
}

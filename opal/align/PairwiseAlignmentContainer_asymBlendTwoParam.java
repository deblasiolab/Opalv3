package opal.align;

import com.traviswheeler.libs.LogWriter;

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
		int v_ext = config.lambda;
		if(j==0) v_ext = config.leftLambdaTerm();
		if(j==N) v_ext = config.rightLambdaTerm();
		float mod;
		if (neighborCnt==0)
			mod = modpair.vLambdasAB[i][j]/normalizer;
		else
			mod = ((1-alpha) * modpair.vLambdasAB[i][j] + alpha * modpair.vLambdas[i][j]/neighborCnt)/normalizer;
		return Math.round( v_ext * (1 + consistency_weight * mod));
	}
	
	protected int calcVGammaOpen(int i, int j, ConsistencyModifiers_Pair modpair, int neighborCnt){
		int v_open = config.gamma/2;
		if(j==0 && i==1) v_open = config.leftGammaTerm()/2;
		if(j==N) v_open = config.rightGammaTerm()/2;
		float mod;
		if (neighborCnt==0)
			mod = modpair.vGammaOpensAB[i][j]/normalizer;
		else
			mod = ((1-alpha) * modpair.vGammaOpensAB[i][j] + alpha * modpair.vGammaOpens[i][j]/neighborCnt)/normalizer;
		return Math.round( v_open * (1 + consistency_weight * mod));

	}
	
	protected int calcVGammaClose(int i, int j, ConsistencyModifiers_Pair modpair, int neighborCnt){
		int v_close = config.gamma/2;
		if(j==0) v_close = config.leftGammaTerm()/2;
		if(j==N&&i==M) v_close = config.rightGammaTerm()/2;
		float mod;
		if (neighborCnt==0)
			mod = modpair.vGammaClosesAB[i][j]/normalizer;
		else
			mod = ((1-alpha) * modpair.vGammaClosesAB[i][j] + alpha * modpair.vGammaCloses[i][j]/neighborCnt)/normalizer;
		return Math.round( v_close * (1 + consistency_weight * mod));
		
	}

	
	protected int calcHLambda(int i, int j, ConsistencyModifiers_Pair modpair, int neighborCnt){
		int h_ext = config.lambda;
		if(i==0) h_ext = config.leftLambdaTerm();
		if(i==M) h_ext = config.rightLambdaTerm();
		float mod;
		if (neighborCnt==0)
			mod = modpair.hLambdasAB[i][j]/normalizer;
		else 
			mod = ((1-alpha) * modpair.hLambdasAB[i][j] + alpha * modpair.hLambdas[i][j]/neighborCnt)/normalizer;
		return Math.round( h_ext * (1 + consistency_weight * mod));
	}
	
	protected int calcHGammaOpen(int i, int j, ConsistencyModifiers_Pair modpair, int neighborCnt){
		int h_open = config.gamma/2;
		if(i==0&&j==1) h_open = config.leftGammaTerm()/2;
		if(i==M) h_open = config.rightGammaTerm()/2;
		float mod;
		if (neighborCnt==0)
			mod = modpair.hGammaOpensAB[i][j]/normalizer;
		else 
			mod = ((1-alpha) * modpair.hGammaOpensAB[i][j] + alpha * modpair.hGammaOpens[i][j]/neighborCnt)/normalizer;
		return Math.round( h_open * (1 + consistency_weight * mod));
	}
	
	protected int calcHGammaClose(int i, int j, ConsistencyModifiers_Pair modpair, int neighborCnt){
		int h_close = config.gamma/2;
		if(i==0) h_close = config.leftGammaTerm()/2;
		if(j==N&i==M) h_close = config.rightGammaTerm()/2;
		float mod;
		if (neighborCnt==0)
			mod = modpair.hGammaClosesAB[i][j]/normalizer;
		else 
			mod = ((1-alpha) * modpair.hGammaClosesAB[i][j] + alpha * modpair.hGammaCloses[i][j]/neighborCnt)/normalizer;
		return Math.round( h_close * (1 + consistency_weight * mod));
	}
	
}

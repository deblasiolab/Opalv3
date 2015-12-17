package opal.align;


import opal.tree.Tree;
import opal.IO.Configuration;

public class PairwiseAlignmentContainer_pureSuboptBlend extends
		PairwiseAlignmentsContainer {

	public PairwiseAlignmentContainer_pureSuboptBlend(Tree tree, float[][] distances, Configuration c) {
		super(tree, distances, c);
	}

	protected void setBlendParams (int neighborCnt) {
		c_weight = neighborCnt==0? 0 : consistency_weight/neighborCnt;
		normalizer = 1/(1+ neighborCnt==0? 0 : consistency_weight);
	}
	
	protected int calcSub(int a, int b, int i, int j, ConsistencyModifiers_Pair modpair, int neighborCnt) {
		float mod = modpair.subsAB[i][j] + c_weight * modpair.subs[i][j];
		return Math.round( normalizer * mod );	
	}
	
	protected int calcVLambda(int i, int j, ConsistencyModifiers_Pair modpair, int neighborCnt) {
		float mod = modpair.vLambdasAB[i][j] + c_weight * modpair.vLambdas[i][j];
		return Math.round( normalizer * mod );
	}
	
	protected int calcVGammaOpen(int i, int j, ConsistencyModifiers_Pair modpair, int neighborCnt){
		float mod = modpair.vGammaOpensAB[i][j] + c_weight * modpair.vGammaOpens[i][j];
		return Math.round( normalizer * mod );
	}
	
	protected int calcVGammaClose(int i, int j, ConsistencyModifiers_Pair modpair, int neighborCnt){
		float mod = modpair.vGammaClosesAB[i][j] + c_weight * modpair.vGammaCloses[i][j];
		return Math.round( normalizer * mod );
	}
	
	
	protected int calcHLambda(int i, int j, ConsistencyModifiers_Pair modpair, int neighborCnt){
		float mod = modpair.hLambdasAB[i][j] + c_weight * modpair.hLambdas[i][j];
		return Math.round( normalizer * mod );
	}
	
	protected int calcHGammaOpen(int i, int j, ConsistencyModifiers_Pair modpair, int neighborCnt){
		float mod = modpair.hGammaOpensAB[i][j] + c_weight * modpair.hGammaOpens[i][j];
		return Math.round( normalizer * mod );
	}
	
	protected int calcHGammaClose(int i, int j, ConsistencyModifiers_Pair modpair, int neighborCnt){
		float mod = modpair.hGammaClosesAB[i][j] + c_weight * modpair.hGammaCloses[i][j];
		return Math.round( normalizer * mod );
	}
		
}

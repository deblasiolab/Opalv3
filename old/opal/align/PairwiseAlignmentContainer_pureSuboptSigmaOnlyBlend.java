package opal.align;

import opal.tree.Tree;


public class PairwiseAlignmentContainer_pureSuboptSigmaOnlyBlend extends
		PairwiseAlignmentsContainer {

	int maxSub;
	
	public PairwiseAlignmentContainer_pureSuboptSigmaOnlyBlend(Tree tree, float[][] distances) {
		super(tree, distances);
	}


	protected void setBlendParams (int neighborCnt) {
		c_weight = neighborCnt==0? 0 : consistency_weight/neighborCnt;
		normalizer = 1/(1+ neighborCnt==0? 0 : consistency_weight);
		maxSub = 0;
	}

	protected int calcSub(int a, int b, int i, int j, ConsistencyModifiers_Pair modpair, int neighborCnt) {
		float mod = modpair.subsAB[i][j] + c_weight * modpair.subs[i][j];
		int ret = Math.round( normalizer * mod );
		if (ret > maxSub)
			maxSub = ret;
		
		return ret;
	}
	protected int calcVLambda(int i, int j, ConsistencyModifiers_Pair modpair, int neighborCnt) {
		return 0;
	}
	
	protected int calcVGammaOpen(int i, int j, ConsistencyModifiers_Pair modpair, int neighborCnt){
		return 0;
	}
	
	protected int calcVGammaClose(int i, int j, ConsistencyModifiers_Pair modpair, int neighborCnt){
		return 0;
	}
	
	protected int calcHLambda(int i, int j, ConsistencyModifiers_Pair modpair, int neighborCnt){
		return 0;
	}
	
	protected int calcHGammaOpen(int i, int j, ConsistencyModifiers_Pair modpair, int neighborCnt){
		return 0;
	}
	
	protected int calcHGammaClose(int i, int j, ConsistencyModifiers_Pair modpair, int neighborCnt){
		return 0;
	}
	
	protected void postProcessAB (int a, int b, ConsistencyModifiers_Pair modpair ) {
 			for (int i=0; i<=M; i++) {
				for (int j= (i==0?1:0); j<=N; j++) {																			
					if (i>0 && j>0) 		
						modpair.subs[i][j] -= maxSub+1;
				}
			}
	}
	
}

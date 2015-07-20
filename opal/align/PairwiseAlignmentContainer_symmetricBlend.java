package opal.align;

import com.traviswheeler.libs.LogWriter;

import opal.IO.SequenceConverter;
import opal.tree.Tree;
import opal.IO.Configuration;

public class PairwiseAlignmentContainer_symmetricBlend extends PairwiseAlignmentsContainer {

	float alpha;
	
	public PairwiseAlignmentContainer_symmetricBlend(Tree tree, float[][] distances, Configuration c) {
		super(tree, distances,c);
		capAtDelta = false;
	}

	protected final void calcPrework () {
		alpha = consistency_other_seqs_weight / (1+ consistency_other_seqs_weight); 
	}
	
	protected final void setBlendParams (int neighborCnt) {
		normalizer = neighborCnt * PairSuboptimalityMatrices.delta;
	}

	protected int calcSub(int a, int b, int i, int j, ConsistencyModifiers_Pair modpair, int neighborCnt) {		
		int sigma = config.cost.costs[origSeqs[a][i-1]][origSeqs[b][j-1]];
		float mod = modpair.subs[i][j]/normalizer;
		return Math.round( sigma * (1 + mod));
	}
	
	protected int calcVLambda(int i, int j, ConsistencyModifiers_Pair modpair, int neighborCnt) {
		int v_ext = (j==0||j==N?config.lambdaTerm : config.lambda);
		float mod =  modpair.vLambdas[i][j]/normalizer;
		return Math.round( v_ext * (1 + mod));
	}
	
	protected int calcVGammaOpen(int i, int j, ConsistencyModifiers_Pair modpair, int neighborCnt){
		int v_open = ( (j==0&&i==1)||j==N ? config.gammaTerm : config.gamma)/2;
		float mod = modpair.vGammaOpens[i][j]/normalizer;
		return Math.round( v_open * (1 + mod));
	}
	
	protected int calcVGammaClose(int i, int j, ConsistencyModifiers_Pair modpair, int neighborCnt){
		int v_close = ( j==0||(j==N&&i==M) ? config.gammaTerm : config.gamma)/2;
		float mod = modpair.vGammaCloses[i][j]/normalizer;
		return Math.round( v_close * (1 + mod));
	}
	
	protected int calcHLambda(int i, int j, ConsistencyModifiers_Pair modpair, int neighborCnt){
		int h_ext = (i==0||i==M?config.lambdaTerm : config.lambda);
		float mod = modpair.hLambdas[i][j]/normalizer;
		return Math.round( h_ext * (1 + mod));
	}
	
	protected int calcHGammaOpen(int i, int j, ConsistencyModifiers_Pair modpair, int neighborCnt){
		int h_open = ( (i==0&&j==1)||i==M ? config.gammaTerm : config.gamma)/2;
		float mod = modpair.hGammaOpens[i][j]/normalizer;
		return Math.round( h_open * (1 + mod));
	}
	
	protected int calcHGammaClose(int i, int j, ConsistencyModifiers_Pair modpair, int neighborCnt){
		int h_close = ( i==0||(j==N&&i==M) ? config.gammaTerm : config.gamma)/2;
		float mod =  modpair.hGammaCloses[i][j]/normalizer;
		return Math.round( h_close * (1 + mod));
	}
	

	protected void postProcessAB (int a, int b, ConsistencyModifiers_Pair modpair ) {

		
		int worstAB_Subopt = 0;
		for (int i=0; i<=M; i++) {
			for (int j= (i==0?1:0); j<=N; j++) {//it's not meaningful to modify 0,0
				if (i>0 && j>0) {
					worstAB_Subopt = Math.max(worstAB_Subopt, modpair.subsAB[i][j]);
				}
				if (i>0) {
					worstAB_Subopt = Math.max(worstAB_Subopt, modpair.vLambdasAB[i][j]);
					if (j>0 || i==1)
						worstAB_Subopt = Math.max(worstAB_Subopt, modpair.vGammaOpensAB[i][j]);
					if (j<N || i==M)
						worstAB_Subopt = Math.max(worstAB_Subopt, modpair.vGammaClosesAB[i][j]);
				}
				if (j>0) {
					worstAB_Subopt = Math.max(worstAB_Subopt, modpair.hLambdasAB[i][j]);
					if (i>0 || j==1)
						worstAB_Subopt = Math.max(worstAB_Subopt, modpair.hGammaOpensAB[i][j]);
					if (i<M || j==N)
						worstAB_Subopt = Math.max(worstAB_Subopt, modpair.hGammaClosesAB[i][j]);
				}
			}
		}

		
		// with the ABC subopt matrix, align AB
		// then get subopt for all positions of AB
		// scale those by (cost of unmodified algnt ;/ cost of mod alignt)
		// use below
		int worstABC_Subopt = calcModifiedModifiers(a, b, modpair);
		// now modpair has been changed ... use it
		
		float ratio_AB_ABC = (float)worstAB_Subopt/worstABC_Subopt; 
		
		
		for (int i=0; i<=M; i++) {
			for (int j= (i==0?1:0); j<=N; j++) {//it's not meaningful to modify 0,0														
				
				if (i>0 && j>0) {
					int sigma = config.cost.costs[origSeqs[a][i-1]][origSeqs[b][j-1]];
					float mod = ((1-alpha) * modpair.subsAB[i][j] + alpha * ratio_AB_ABC * modpair.subs[i][j])/PairSuboptimalityMatrices.delta;
					modpair.subs[i][j] = Math.round( sigma * (1 + consistency_weight * mod));
				}
				
				if (i>0) {
					int v_ext = (j==0||j==N?config.lambdaTerm : config.lambda);
					float mod = ((1-alpha) * modpair.vLambdasAB[i][j] + alpha * ratio_AB_ABC * modpair.vLambdas[i][j])/PairSuboptimalityMatrices.delta;
					modpair.vLambdas[i][j] = Math.round( v_ext * (1 + consistency_weight * mod));

					int v_open = ( (j==0&&i==1)||j==N ? config.gammaTerm : config.gamma)/2;
					mod = ((1-alpha) * modpair.vGammaOpensAB[i][j] + alpha * ratio_AB_ABC * modpair.vGammaOpens[i][j])/PairSuboptimalityMatrices.delta;
					modpair.vGammaOpens[i][j] = Math.round( v_open * (1 + consistency_weight * mod));

					int v_close = ( j==0||(j==N&&i==M) ? config.gammaTerm : config.gamma)/2;
					mod = ((1-alpha) * modpair.vGammaClosesAB[i][j] + alpha * ratio_AB_ABC * modpair.vGammaCloses[i][j])/PairSuboptimalityMatrices.delta;
					modpair.vGammaCloses[i][j] = Math.round( v_close * (1 + consistency_weight * mod));
				}
				
				if (j>0) {
					int h_ext = (i==0||i==M?config.lambdaTerm : config.lambda);
					float mod = ((1-alpha) * modpair.hLambdasAB[i][j] + alpha * ratio_AB_ABC * modpair.hLambdas[i][j])/PairSuboptimalityMatrices.delta;
					modpair.hLambdas[i][j] = Math.round( h_ext * (1 + consistency_weight * mod));

					int h_open = ( (i==0&&j==1)||i==M ? config.gammaTerm : config.gamma)/2;
					mod = ((1-alpha) * modpair.hGammaOpensAB[i][j] + alpha * ratio_AB_ABC * modpair.hGammaOpens[i][j])/PairSuboptimalityMatrices.delta;
					modpair.hGammaOpens[i][j] = Math.round( h_open * (1 + consistency_weight * mod));

					int h_close = ( i==0||(j==N&&i==M) ? config.gammaTerm : config.gamma)/2;
					mod = ((1-alpha) * modpair.hGammaClosesAB[i][j] + alpha * ratio_AB_ABC * modpair.hGammaCloses[i][j])/PairSuboptimalityMatrices.delta;
					modpair.hGammaCloses[i][j] = Math.round( h_close * (1 + consistency_weight * mod));
				}			
			}
		}
		
	}
	

	
	
	
	private int calcModifiedModifiers(int a, int b, ConsistencyModifiers_Pair modpair) {

		int M = origSeqs[a].length;
		int N = origSeqs[b].length;
		
		PairSuboptimalityMatricesModified mat = new PairSuboptimalityMatricesModified(
				Alignment.buildNewAlignment(origSeqs[a], a, config), 
				Alignment.buildNewAlignment(origSeqs[b], b, config),
				modpair, config);

		
		// because the ABC-modified alignmt is based on costs that are larger further from the 
		// optimal alignt, we need to scale the final costs down to keep away from the case that
		// all cells not on the optimal alignt path have maximum suboptimality
			

		int worstABC_Subopt = 0;

		for (int i=0; i<=M; i++) {
			for (int j= (i==0?1:0); j<=N; j++) {//it's not meaningful to modify 0,0
						
				if (i>0 && j>0) { 
					modpair.subs[i][j] = (int)(mat.getSubsSubopt(i,j));
					worstABC_Subopt = Math.max(worstABC_Subopt, modpair.subs[i][j]);
				} else { 
					modpair.subs[i][j] = bigBadVal;
				}
				
				if (i>0) {
					modpair.vLambdas[i][j] = (int)(mat.getVExtSubopt(i, j));
					modpair.vGammaOpens[i][j] = (int)(mat.getVOpenSubopt(i, j));
					modpair.vGammaCloses[i][j] = (int)(mat.getVCloseSubopt(i, j));
					worstABC_Subopt = Math.max(worstABC_Subopt, modpair.vLambdas[i][j]);
					if (j>0 || i==1)
						worstABC_Subopt = Math.max(worstABC_Subopt, modpair.vGammaOpens[i][j]);
					if (j<N || i==M)
						worstABC_Subopt = Math.max(worstABC_Subopt, modpair.vGammaCloses[i][j]);
				} else {
					modpair.vLambdas[i][j] = modpair.vGammaOpens[i][j] = modpair.vGammaCloses[i][j] = bigBadVal;  
				}
				
				if (j>0) {
					modpair.hLambdas[i][j] = (int)(mat.getHExtSubopt(i, j));
					modpair.hGammaOpens[i][j] = (int)(mat.getHOpenSubopt(i, j));					
					modpair.hGammaCloses[i][j] =  (int)(mat.getHCloseSubopt(i, j));
					worstABC_Subopt = Math.max(worstABC_Subopt, modpair.hLambdas[i][j]);
					if (i>0 || j==1)
						worstABC_Subopt = Math.max(worstABC_Subopt, modpair.hGammaOpens[i][j]);
					if (i<M || j==N)
						worstABC_Subopt = Math.max(worstABC_Subopt, modpair.hGammaCloses[i][j]);
				} else {
					modpair.hLambdas[i][j] = modpair.hGammaOpens[i][j] = modpair.hGammaCloses[i][j] = bigBadVal;  
				}
				
			}
		}	
		return worstABC_Subopt;
	}
	
}

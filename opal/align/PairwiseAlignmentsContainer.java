package opal.align;

import java.util.LinkedList;

import java.util.Hashtable;
import java.util.Stack;

import opal.exceptions.GenericOpalException;
import opal.tree.Tree;
import opal.tree.TreeNode;
import opal.IO.Configuration;

import com.traviswheeler.libs.LogWriter;

public abstract class PairwiseAlignmentsContainer {

	public static int neighborCount = 5;
	public static float consistency_weight = 1;
	public static float consistency_other_seqs_weight = 1;
	public static int maxSubtreeSize = 6;
	public static double badScoreMult = 1.05;
	public static boolean useMax = true;
	public static boolean useWeights = false;
	public static float flattenABC = -1;
	
	public static int bigBadVal = Integer.MAX_VALUE/1000;
	
	public static enum RangeOverlapType { intersect, union, union_fast, intersect_fast  };
	RangeOverlapType rangeMethod = RangeOverlapType.intersect_fast;
	
	public static enum BlendType { simple, asym_oneparam, asym_twoparam, symmetric  };
	public static BlendType blendMethod = BlendType.symmetric;

	public static enum NeighborDistType { max, avg  };
	public static NeighborDistType neighborDistMethod = NeighborDistType.max;

	public static float neighborDistThreshold = 1000;
	
	protected boolean capAtDelta = true;
	
	Configuration config;
	
	float[][] distances;
	long[][] pairCosts;
	int[][] origSeqs;

	float c_weight;
	float normalizer;
	int M, N;

	
	//for each sequencePair, store a list of neighbors and a list of the distances of those neighbors
	public Hashtable<String,LinkedList<Integer>> nearestNeighbors = 
								new Hashtable<String,LinkedList<Integer>>();
	public Hashtable<String,LinkedList<Float>> neighborDistances =
								new Hashtable<String,LinkedList<Float>>();
	public Hashtable<String,PairSuboptimalityMatrices[]> pairSuboptMatrices = 
								new Hashtable<String,PairSuboptimalityMatrices[]>();
	public Hashtable<String,Integer> pairUseCount = new Hashtable<String,Integer>();
	
	TreeNode root;
	
	public PairwiseAlignmentsContainer (Tree tree, float distances[][], Configuration c) {
		config = c;
		this.distances = distances;		
		root = tree.getRoots().get(0);
				
		//get the input (ungapped) sequences
		int K = root.lastLeaf - root.firstLeaf + 1;
		TreeNode[] orig2leaves = new TreeNode[K];  // for testing only
		TreeNode[] leaves = new TreeNode[K];
		origSeqs = new int[K][];
			
		root.fillLeafList(0, leaves);
		for (int p=0; p<leaves.length; p++) {
			int a = getInputSeqID(p);
			origSeqs[a] = leaves[p].alignment.seqs[0];
			orig2leaves[a] = leaves[p];
		}
		
		pairCosts = new long[distances.length][distances[0].length];
		for(int p=0; p<distances.length; p++){
			for(int q=p+1; q<distances[0].length; q++){
				pairCosts[p][q] = pairCosts[q][p] = (int)(0.5 + distances[p][q] * (origSeqs[p].length + origSeqs[q].length));
			}
		}

		
		//preprocess: get list of neighbors, and tally pair-counts
		for (TreeNode node : tree.getNodesToMerge()) {
			int first = node.firstLeaf;
			int size = node.lastLeaf - first + 1;
			if (size <= maxSubtreeSize ) {
				TreeNode left = node.leftChild;
				TreeNode right = node.rightChild;
				for (int i = left.firstLeaf; i<=left.lastLeaf; i++) {
					int a = getInputSeqID(i);
					for (int j = right.firstLeaf; j<=right.lastLeaf; j++) {
						int b = getInputSeqID(j);
						tallyNeighbors(a,b);		
					}
				}				
			}
		}
	
		for (String key: neighborDistances.keySet()) {
			LinkedList<Float> dists = neighborDistances.get(key);
			float sum = 0;

			if (dists.size() == 0) continue;
					
			int[] ab = SequenceIdPair.splitString(key);
			/*float ab_dist = distances[ab[0]][ab[1]];
			if (orig2leaves[ab[0]].parent == orig2leaves[ab[1]].parent ) {
				LogWriter.stdErrLog(ab[0] + "," + ab[1] + " (" + ab_dist + ") : " );
				LinkedList<Integer> neighbors = nearestNeighbors.get(key);
				for (int i=0; i<dists.size(); i++) {
					LogWriter.stdErrLog(neighbors.get(i) + "(" + dists.get(i) + "); " ); 
				}
				LogWriter.stdErrLogln("");
			}
			*/
						
			if (useWeights) {
//				int[] ab = SequenceIdPair.splitString(key);
				float self_dist = Math.max(distances[ab[0]][ab[0]], distances[ab[1]][ab[1]]);
	
				for (int i=0; i<dists.size(); i++) {
					float f = dists.get(i);
					float x = f - self_dist;
					sum += 1/x;
					dists.set(i,  x);
				}
				for (int i=0; i<dists.size(); i++) {
					float f = dists.get(i);
					f = (1/f)/sum;
					dists.set(i,  f);
				}

			} else {
				for (int i=0; i<dists.size(); i++) {
					dists.set(i, (float)1/neighborCount);
				}
			}
		}
		
		this.distances = null;
	}

	
	private void tallyNeighbors (int a, int b) {
		float neighborDist;
		String pair = SequenceIdPair.makeString(a, b);
		
		for (int k=0 ; k<root.leafOrderFromInput.size(); k++) {
			int c = root.leafOrderFromInput.get(k);
			if (c==a || c==b) continue;
			neighborDist = neighborDistMethod == NeighborDistType.max ? 
											Math.max(distances[a][c], distances[b][c]) : 
											(distances[a][c] + distances[b][c]) / 2 ;	 			
			placeNeighborInList(pair, c, neighborDist);
		}

		LinkedList<Float> neighborDists = neighborDistances.get(pair);
		LinkedList<Integer> neighbors = nearestNeighbors.get(pair);
		if (Float.MAX_VALUE == neighborDists.getLast()) { 
			// fewer than X seqs, total. Need to remove the placeholder I added at first.
			neighborDists.removeLast();
			neighbors.removeLast();
		}

		for (int c: neighbors) {
			pair = SequenceIdPair.makeString(a, c);
			Integer x = pairUseCount.get(pair);
			pairUseCount.put(pair, 1 + (x == null ? 0 : x.intValue()) );
			
			pair = SequenceIdPair.makeString(b, c);
			x = pairUseCount.get(pair);
			pairUseCount.put(pair, 1 + (x == null ? 0 : x.intValue()) );
		}
		
	}
	
	private void placeNeighborInList (String pair, int c, float dist) {
	
		
		LinkedList<Float> neighborDists = neighborDistances.get(pair);
		LinkedList<Integer> neighbors;
		if (neighborDists == null) {
			//initialize lists with bogus entries that should run off the end of the list soon
			neighborDists = new LinkedList<Float>();
			neighborDists.add(Float.MAX_VALUE);
			neighborDistances.put(pair, neighborDists);

			neighbors = new LinkedList<Integer>();
			neighbors.add(-1);
			nearestNeighbors.put(pair, neighbors);
		} else {
			neighbors = nearestNeighbors.get(pair);	
		}
		
		if ( dist < neighborDistThreshold) {
			//run through distance list. if reach end, don't add.
			//if find a distance that's > dist, insert dist and neighbor there, and possibly remove last entry on list.
			for (int i=0; i<neighborDists.size(); i++) {
				if (neighborDists.get(i)>dist) {
					//then insert before i
					neighborDists.add(i, dist);
					neighbors.add(i,c);
					//remove last entry, if we're over the threshold
					if (neighborDists.size() > neighborCount) {
						neighborDists.removeLast();
						neighbors.removeLast();
					}				
					break;
				}
			}
		}
	}

	protected void calcAB_Subopts (int a, int seq_i, int b, int seq_j, ConsistencyModifiers_AllPairs mods) {

		//build the A-B alignment, and find subopt ranges
		String pair = SequenceIdPair.makeString(a,b);
		PairSuboptimalityMatrices[] matAB = pairSuboptMatrices.get(pair);
		
		if ( matAB == null){
			Alignment A = Alignment.buildNewAlignment(origSeqs[a], a, config);
			Alignment B = Alignment.buildNewAlignment(origSeqs[b], b, config);
			
			/*
			 * 
			 * Not sure whats happeneing here
			matAB = new PairSuboptimalityMatrices[Aligner.lambdas.length];
			for (int z=0; z<Aligner.lambdas.length; z++) {
				Aligner.switchParamID(z);
				matAB[z] = new PairSuboptimalityMatrices(A,B,capAtDelta);
			}
			Aligner.switchParamID(0);
			*/
			pairSuboptMatrices.put(pair, matAB);
		}
		//Get the correct cell in the allpairs matrix
		ConsistencyModifiers_Pair modpair = mods.modifiers[seq_i][seq_j];
		
		int badScore;
		if (capAtDelta)
			badScore = (int)(badScoreMult * PairSuboptimalityMatrices.delta );
		else
			badScore = Integer.MAX_VALUE;
		//Date start = new Date();
		
		for (int i=0; i<=M; i++) {
			for (int j= (i==0?1:0); j<=N; j++) {//it's not meaningful to modify 0,0
				modpair.subsAB[i][j] = 
					modpair.vLambdasAB[i][j] = 
					modpair.vGammaOpensAB[i][j] = 
					modpair.vGammaClosesAB[i][j] = 
					modpair.hLambdasAB[i][j] = 
					modpair.hGammaOpensAB[i][j] = 					
					modpair.hGammaClosesAB[i][j] = 0;

				if (i==0) 
					modpair.subsAB[i][j] = modpair.vLambdasAB[i][j] = modpair.vGammaOpensAB[i][j] = modpair.vGammaClosesAB[i][j] = bigBadVal;  
				if (j==0)
					modpair.subsAB[i][j] = modpair.hLambdasAB[i][j] = modpair.hGammaOpensAB[i][j] = modpair.hGammaClosesAB[i][j] = bigBadVal;  
			
				
				for (int z=0; z<matAB.length; z++) {
					// Aligner.switchParamID(z);

					if (capAtDelta && (j < matAB[z].nearlyOptimal_horizStart[i] || j > matAB[z].nearlyOptimal_horizEnd[i])) { 
						modpair.subsAB[i][j] += badScore;
						modpair.vLambdasAB[i][j] += badScore;
						modpair.vGammaOpensAB[i][j] += badScore;
						modpair.vGammaClosesAB[i][j] += badScore;
						modpair.hLambdasAB[i][j] += badScore;
						modpair.hGammaOpensAB[i][j] += badScore;			
						modpair.hGammaClosesAB[i][j] += badScore;
					} else {
						if (i>0 && j>0) {
							modpair.subsAB[i][j] += Math.min(badScore, matAB[z].getSubsSubopt(i,j));
						}
						if (i>0) {
							modpair.vLambdasAB[i][j] += Math.min(badScore, matAB[z].getVExtSubopt(i, j));
							modpair.vGammaOpensAB[i][j] += Math.min(badScore, matAB[z].getVOpenSubopt(i, j));
							modpair.vGammaClosesAB[i][j] += Math.min(badScore, matAB[z].getVCloseSubopt(i, j));
						}
			
						if (j>0) {
							modpair.hLambdasAB[i][j] += Math.min(badScore, matAB[z].getHExtSubopt(i, j));
							modpair.hGammaOpensAB[i][j] += Math.min(badScore, matAB[z].getHOpenSubopt(i, j));					
							modpair.hGammaClosesAB[i][j] += Math.min(badScore, matAB[z].getHCloseSubopt(i, j));
						}
					}
				}	
				//Aligner.switchParamID(0);

			}
		}
	}


	protected void calcAB_Modifiers (int a, int seq_i, int b, int seq_j, int c, ConsistencyModifiers_AllPairs mods, float weight) {
		 //weight *= neighborCount; // avoid numerical errors
		
		//build the A-C & B-C alignments, and find subopt ranges
		String pair = SequenceIdPair.makeString(a,c);
			
		
		PairSuboptimalityMatrices[] matAC = pairSuboptMatrices.get(pair);
		if ( matAC == null){
			Alignment A = Alignment.buildNewAlignment(origSeqs[a], a, config);
			Alignment C = Alignment.buildNewAlignment(origSeqs[c], c, config);
			/*
			 * not sure whats happeneing here
			
			matAC = new PairSuboptimalityMatrices[Aligner.lambdas.length];
			for (int x=0; x<Aligner.lambdas.length; x++) {
				Aligner.switchParamID(x);
				matAC[x] = new PairSuboptimalityMatrices(A,C);
			}
			pairSuboptMatrices.put(pair, matAC);
			*/
		}

			
		pair = SequenceIdPair.makeString(b,c);
		PairSuboptimalityMatrices matBC[] = pairSuboptMatrices.get(pair);
		if ( matBC == null){
			Alignment B = Alignment.buildNewAlignment(origSeqs[b], b, config);
			Alignment C = Alignment.buildNewAlignment(origSeqs[c], c, config);
			
			/*
			 * Again this
			matBC = new PairSuboptimalityMatrices[Aligner.lambdas.length];
			for (int x=0; x<Aligner.lambdas.length; x++) {
				Aligner.switchParamID(x);
				matBC[x] = new PairSuboptimalityMatrices(B,C);
			}
			Aligner.switchParamID(0);
			pairSuboptMatrices.put(pair, matBC);
			*/
		}

		
		//Get the correct cell in the allpairs matrix
		ConsistencyModifiers_Pair modpair = mods.modifiers[seq_i][seq_j];
		int P = origSeqs[c].length;

		int badScore = Math.round((float)(badScoreMult * PairSuboptimalityMatrices.delta ));
				
		//Date start = new Date();	
		
		for (int i=0; i<=M; i++) {
			/*
			if (i%50 == 0 ) {
				Date now = new Date();
				long diff = now.getTime() - start.getTime();
				LogWriter.stdErrLogln("\tstart row " + i + " : " + diff + " ms");
				start = now;
			}
			*/
			for (int j= (i==0?1:0); j<=N; j++) {//it's not meaningful to modify 0,0?
								
				for (int z=0; z<matAC.length; z++) {
					// turned off Aligner.switchParamID(z);

//for each matAC,matBC pair (one for each parameter choice)				
				
					//pick positions k to test against i and j
					Stack<Integer> positions = new Stack<Integer>();
					PairSuboptimalityMatrices tmpMatA, tmpMatB;
					int ii,jj;
					if (matAC[z].nearlyOptimal_horizStart[i] <= matBC[z].nearlyOptimal_horizStart[j]) {
						tmpMatA = matAC[z]; 
						tmpMatB = matBC[z];
						ii = i; jj = j;
					} else {
						tmpMatA = matBC[z];
						tmpMatB = matAC[z];
						ii = j; jj = i;
					}
	
					if ( (rangeMethod == RangeOverlapType.intersect || rangeMethod == RangeOverlapType.intersect_fast) &&
							tmpMatA.nearlyOptimal_horizEnd[ii] >= tmpMatB.nearlyOptimal_horizStart[jj]) {
						//intersection isn't empty
						for (int k=tmpMatB.nearlyOptimal_horizStart[jj]; k<=tmpMatA.nearlyOptimal_horizEnd[ii] && k<=tmpMatB.nearlyOptimal_horizEnd[jj]; k++)
							positions.add(k);
					} else if (rangeMethod == RangeOverlapType.union || 
							(rangeMethod == RangeOverlapType.union_fast  
													&& tmpMatA.nearlyOptimal_horizEnd[ii] >= tmpMatB.nearlyOptimal_horizStart[jj])) {					
						if (tmpMatA.nearlyOptimal_horizEnd[ii] >= tmpMatB.nearlyOptimal_horizStart[jj]) {
							//union, and they overlap
							for (int k=tmpMatA.nearlyOptimal_horizStart[ii]; k<=tmpMatB.nearlyOptimal_horizEnd[jj]; k++)
								positions.add(k);				
						} else { // union, not overlapping
							for (int k=tmpMatA.nearlyOptimal_horizStart[ii]; k<=tmpMatA.nearlyOptimal_horizEnd[ii]; k++)
								positions.add(k);
							for (int k=tmpMatB.nearlyOptimal_horizStart[jj]; k<=tmpMatB.nearlyOptimal_horizEnd[jj]; k++)
								positions.add(k);
						}
					} else { // union_fast or intersect_fast, not overlapping
						// don't even bother with detailed calculations; just fill in default (high) values
						int bad;
						if (tmpMatB.nearlyOptimal_horizStart[jj] - tmpMatA.nearlyOptimal_horizEnd[ii] <= 10)
							bad = badScore;
						else
						//					else if (tmpMatB.nearlyOptimal_horizStart[jj] - tmpMatA.nearlyOptimal_horizEnd[ii] <= 100)
							bad = (int)( 1.1* badScoreMult * PairSuboptimalityMatrices.delta );
	
						//else 
							//bad = bigBadVal; // that far apart? no way this alignment option makes sense
						
						modpair.subs[i][j] += weight * bad;
						modpair.vLambdas[i][j] += weight * bad;
						modpair.hLambdas[i][j] += weight * bad;
						modpair.vGammaOpens[i][j] += weight * bad;
						modpair.hGammaOpens[i][j] += weight * bad;
						modpair.vGammaCloses[i][j] += weight * bad;				
						modpair.hGammaCloses[i][j] += weight * bad;
						
						continue;
					}
									
	
					if (positions.size()==1){
						int k=positions.peek();
						if (k==0)
							positions.add(1);
						else if (k==P)
							positions.add(0, P-1); 
					}
									
					int bestSub = bigBadVal;
					int bestExtHoriz = bigBadVal;
					int bestExtVert = bigBadVal;
					int bestOpenHoriz = bigBadVal;
					int bestOpenVert = bigBadVal;
					int bestCloseHoriz = bigBadVal;
					int bestCloseVert = bigBadVal;
					
					int best_k_prime = -1;
					//for the special case of gamma where we avoid quadratic time through keeping track of best preceding value
					//Reminder: this ony matters for k>2, and with runs of positions at least 3 long
					if ( i>0 && j>0  ) {  // the case where i or j == 0 doesn't require the k/k' solution, so it's handled later
						int mink = bigBadVal;
						int lastk = -1;
						
						for (int k : positions){
							if (k>lastk+1){ //disconnect between the two ranges
								mink=bigBadVal;
								best_k_prime = -1;
							}
							if (mink==bigBadVal) {
								lastk = mink=k;
								continue;
							}
							if (k==mink+1) {
								lastk = k;
								continue;
							}
							lastk=k;
	
							int x,y,xy1,xy2;
							
							/*Gamma opens*/
							//vertical  (a,-)
							xy1 = bigBadVal;
							if (best_k_prime>0 /*best_k_prime!=-1*/) {
								x = matAC[z].getOpenSubopt_special(i,best_k_prime,k,true);
								y = matBC[z].getOpenSubopt_special(j,best_k_prime,k,false);
								xy1 =  (useMax ? Math.max(x,y) : x + y );
								if ( xy1 < bestOpenVert) {
									bestOpenVert = xy1;
								}
							}
							
							if (k>2) {
								x = matAC[z].getOpenSubopt_special(i,k-2,k,true);
								y = matBC[z].getOpenSubopt_special(j,k-2,k,false);
								xy2 =  (useMax ? Math.max(x,y) : x + y );
								if (xy2 < xy1) {
									best_k_prime = k-2;
									if ( xy2 < bestOpenVert) {
										bestOpenVert = xy2;
									}
								}
							}
							
							//horiz (-,b)
							xy1 = bigBadVal;
							if (best_k_prime>0) {
								x = matBC[z].getOpenSubopt_special(j,best_k_prime,k,true);
								y = matAC[z].getOpenSubopt_special(i,best_k_prime,k,false);
								xy1 =  (useMax ? Math.max(x,y) : x + y );
								if ( xy1 < bestOpenHoriz) {
									bestOpenHoriz = xy1;
								}
							}
							
							if (k>2) {
								x = matBC[z].getOpenSubopt_special(j,k-2,k,true);
								y = matAC[z].getOpenSubopt_special(i,k-2,k,false);
								xy2 =  (useMax ? Math.max(x,y) : x + y );
								if (xy2 < xy1) {
									best_k_prime = k-2;
									if ( xy2 < bestOpenHoriz) {
										bestOpenHoriz = xy2;
									}
								}
							}
						
							/*Gamma closes*/
							//vertical  (a,-)
							xy1 = bigBadVal;
							if (best_k_prime>0 /*best_k_prime!=-1*/) {
								x = matAC[z].getCloseSubopt_special(i,best_k_prime,k,true);
								y = matBC[z].getCloseSubopt_special(j,best_k_prime,k,false);
								xy1 =  (useMax ? Math.max(x,y) : x + y );
								if ( xy1 < bestCloseVert) {
									bestCloseVert = xy1;
								}
							}
	
							if (k>2) {
								x = matAC[z].getCloseSubopt_special(i,k-2,k,true);
								y = matBC[z].getCloseSubopt_special(j,k-2,k,false);
								xy2 =  (useMax ? Math.max(x,y) : x + y );
								if (xy2 < xy1) {
									best_k_prime = k-2;
									if ( xy2 < bestCloseVert) {
										bestCloseVert = xy2;
									}
								}
							}
							
							//horiz (-,b)
							xy1 = bigBadVal;
							if (best_k_prime>0 /*best_k_prime!=-1*/) {
								x = matBC[z].getCloseSubopt_special(j,best_k_prime,k,true);
								y = matAC[z].getCloseSubopt_special(i,best_k_prime,k,false);
								xy1 =  (useMax ? Math.max(x,y) : x + y );
								if ( xy1 < bestCloseHoriz) {
									bestCloseHoriz = xy1;
								}
							}
							
							if (k>2) {
								x = matBC[z].getCloseSubopt_special(j,k-2,k,true);
								y = matAC[z].getCloseSubopt_special(i,k-2,k,false);
								xy2 =  (useMax ? Math.max(x,y) : x + y );
								if (xy2 < xy1) {
									best_k_prime = k-2;
									if ( xy2 < bestCloseHoriz) {
										bestCloseHoriz = xy2;
									}
								}
							}
						}
	
					}
					
					
					for (int k : positions){			
						int x,y,xy;
	
						if (i>0 && j>0 && k>0) {
							x = matAC[z].getSubsSubopt(i,k);
							y = matBC[z].getSubsSubopt(j,k);
							xy =  (useMax ? Math.max(x,y) : x + y );
							if ( xy < bestSub) {
								bestSub = xy;
							}
						}					
						
						/*Lambdas*/
						//vertical  (a,-)
						if (i>0) {
							if ( (k==P && j<N) || (k==0 && j>0) ) {
									bestExtVert = badScore * (useMax ? 1 : 2);
							} else {
								//either k<P or k==P&&j==N ...   and k>0 or k==0&&j==0
								x = matAC[z].getExtSubopt(i,k,true);
								y = matBC[z].getExtSubopt(j,k,false);
								xy =  (useMax ? Math.max(x,y) : x + y );
								if ( xy < bestExtVert) {
									bestExtVert = xy;
								}
							}
						}
					
						//horiz (-,b)
						if (j>0) {
							if ( (k==P && i<M) || (k==0 && i>0) ) {
								bestExtHoriz = badScore * (useMax ? 1 : 2);
							} else {
								x = matBC[z].getExtSubopt(j,k,true);
								y = matAC[z].getExtSubopt(i,k,false);
								xy =  (useMax ? Math.max(x,y) : x + y );
								if ( xy < bestExtHoriz) {
									bestExtHoriz = xy;
								}
							}
						}
						
						
						/*Gamma opens*/
						//vertical  (a,-)
						if ( i>0 ) {
							if ( j==0 && (i>1 || k==P) ) { // this value won't be used, because you can't open a gap at i,0 with i>1
								x = y = badScore;
							} else if (j==0&&i==1) { // special case: k=0 is ok for subopt2, but not subopt1;  k>=1 valid for both  
								if (k==0) {
									x = matAC[z].getOpenSubopt2(i,k,true);
									y = matBC[z].getOpenSubopt2(j,k,false);
								} else {
									x = matAC[z].getOpenSubopt(i,k,true);
									y = matBC[z].getOpenSubopt(j,k,false);								
								}
							} else {  // j>0
								if (k==0) {
									x = y = badScore;
								} else if (k==P && j<N) {
									x = matAC[z].getOpenSubopt1(i,k,true);
									y = matBC[z].getOpenSubopt1(j,k,false);
								} else {
									x = matAC[z].getOpenSubopt(i,k,true);
									y = matBC[z].getOpenSubopt(j,k,false);
								}
							}
							xy =  (useMax ? Math.max(x,y) : x + y );
							if ( xy < bestOpenVert) {
								bestOpenVert = xy;
							}
						}
						
						//horiz (-,b)
						if ( j>0 ) {
							if ( i==0 && (j>1 || k==P)) { // this value won't be used, because you can't open a gap at 0,j with j>1
								x = y = badScore;
							} else if (i==0&&j==1) { // special case: k=0 is ok for subopt2, but not subopt1;  k>=1 valid for both  
								if (k==0) {
									x = matBC[z].getOpenSubopt2(j,k,true);
									y = matAC[z].getOpenSubopt2(i,k,false);
								} else {
									x = matBC[z].getOpenSubopt(j,k,true);
									y = matAC[z].getOpenSubopt(i,k,false);								
								}
							} else {  // i>0
								if (k==0) {
									x = y = badScore;
								} else if (k==P && i<M) {
									x = matBC[z].getOpenSubopt1(j,k,true);
									y = matAC[z].getOpenSubopt1(i,k,false);
								} else {
									x = matBC[z].getOpenSubopt(j,k,true);
									y = matAC[z].getOpenSubopt(i,k,false);
								}
							}
							xy =  (useMax ? Math.max(x,y) : x + y );
							if ( xy < bestOpenHoriz) {
								bestOpenHoriz = xy;
							}
						}
	
	
						/*Gamma closes*/
						//vertical  (a,-)
						if ( i>0 ) {
							if (k==0 && j>0) {
								x = y = badScore;
							} else if (j==N && i<M) {
								x = y = badScore;
							} else if (j==N && i==M) {
								x = matAC[z].getCloseSubopt(i,k,true);
								y = matBC[z].getCloseSubopt(j,k,false);
							} else { // j<N,  and k==0 only if j==0, 
								if (k==P) {
									x = y = badScore;
								} else if (k==0 /* and j==0*/){
									x = matAC[z].getCloseSubopt2(i,k,true);
									y = matBC[z].getCloseSubopt2(j,k,false);								
								} else {
									x = matAC[z].getCloseSubopt(i,k,true);
									y = matBC[z].getCloseSubopt(j,k,false);								
								}
							}
							xy =  (useMax ? Math.max(x,y) : x + y );
							if ( xy < bestCloseVert) {
								bestCloseVert = xy;
							}
						}
						
						
						
						//horiz (-,b)
						if ( j>0 ) {
							if (k==0 && i>0) {
								x = y = badScore;
							} else if (i==M && j<N) {
								x = y = badScore;
							} else if (i==M && j==N) {
								x = matBC[z].getCloseSubopt(j,k,true);
								y = matAC[z].getCloseSubopt(i,k,false);
							} else { // i<M,  and k==0 only if i==0, 
								if (k==P) {
									x = y = badScore;
								} else if (k==0 /* and i==0*/){
									x = matBC[z].getCloseSubopt2(j,k,true);
									y = matAC[z].getCloseSubopt2(i,k,false);								
								} else {
									x = matBC[z].getCloseSubopt(j,k,true);
									y = matAC[z].getCloseSubopt(i,k,false);								
								}
							}
							xy =  (useMax ? Math.max(x,y) : x + y );
							if ( xy < bestCloseHoriz) {
								bestCloseHoriz = xy;
							}
						}				
					}
								
	
					if ( (bestSub == bigBadVal && i>0 && j>0) ||
						(bestExtVert == bigBadVal && i>0) ||
						(bestExtHoriz == bigBadVal && j>0) ||
						(bestOpenVert == bigBadVal && i>0) ||
						(bestOpenHoriz == bigBadVal &&j>0) ||
						(bestCloseVert == bigBadVal && i>0) ||
						(bestCloseHoriz == bigBadVal && j>0) ) {
							LogWriter.stdErrLogln("unexpected high value in pw-alignt-cont");
							throw new GenericOpalException("");
					}
	
	
					if ( bestSub < 0 ||
						bestExtVert < 0 ||
						bestExtHoriz < 0 ||
						bestOpenVert < 0 ||
						bestOpenHoriz < 0 ||
						bestCloseVert < 0 ||
						bestCloseHoriz < 0 ){
							LogWriter.stdErrLogln("unexpected low value in pw-alignt-cont");
							throw new GenericOpalException("");
					}
	
					
					modpair.subs[i][j] += weight * Math.min(badScore, bestSub);
					modpair.vLambdas[i][j] += weight * Math.min(badScore, bestExtVert);
					modpair.hLambdas[i][j] += weight * Math.min(badScore, bestExtHoriz);
					modpair.vGammaOpens[i][j] += weight * Math.min(badScore, bestOpenVert);
					modpair.hGammaOpens[i][j] += weight * Math.min(badScore, bestOpenHoriz);
					modpair.vGammaCloses[i][j] += weight * Math.min(badScore, bestCloseVert);				
					modpair.hGammaCloses[i][j] += weight * Math.min(badScore, bestCloseHoriz);
					
				}
				// turned off Aligner.switchParamID(0);
			}
		}
	}
	
	public int getInputSeqID (int i) {
		return root.leafOrderFromInput.get(i);
	}
	
	protected void cleanHashes (String pair) {
		int x = pairUseCount.get(pair).intValue() - 1;
		if (x == 0) {
			pairUseCount.remove(pair);
			pairSuboptMatrices.remove(pair); 
		} else {
			pairUseCount.put(pair, x);	
		}	
	}
	
	protected void calcPrework () {}

	protected abstract void setBlendParams (int neighborCnt);

	protected abstract int calcSub(int a, int b, int i, int j, ConsistencyModifiers_Pair modpair, int neighborCnt);
	
	protected abstract int calcVLambda(int i, int j, ConsistencyModifiers_Pair modpair, int neighborCnt);
	protected abstract int calcVGammaOpen(int i, int j, ConsistencyModifiers_Pair modpair, int neighborCnt);
	protected abstract int calcVGammaClose(int i, int j, ConsistencyModifiers_Pair modpair, int neighborCnt);

	protected abstract int calcHLambda(int i, int j, ConsistencyModifiers_Pair modpair, int neighborCnt);
	protected abstract int calcHGammaOpen(int i, int j, ConsistencyModifiers_Pair modpair, int neighborCnt);
	protected abstract int calcHGammaClose(int i, int j, ConsistencyModifiers_Pair modpair, int neighborCnt);
	
	
	public void calcAllPairModifiers (TreeNode node, ConsistencyModifiers_AllPairs mods) {
		TreeNode left = node.leftChild;
		TreeNode right = node.rightChild;
					
		calcPrework();
		
		for (int p = left.firstLeaf; p<=left.lastLeaf; p++) {
			int a = root.leafOrderFromInput.get(p);
			for (int q = right.firstLeaf; q<=right.lastLeaf; q++) {
				int b = root.leafOrderFromInput.get(q);
				
				M = origSeqs[a].length;
				N = origSeqs[b].length;

				calcAB_Subopts (a, p-left.firstLeaf, b, q-right.firstLeaf, mods);
				
				String pair = SequenceIdPair.makeString(a,b);
								
				LinkedList<Integer> neighbors = nearestNeighbors.get(pair);
				LinkedList<Float> weights = neighborDistances.get(pair);
				
				
				for (int x=0; x<neighbors.size(); x++) {
					int k = neighbors.get(x);
					float weight = weights.get(x) * neighbors.size(); // push value sup high enough to avoid funny rounding effects 
	
					calcAB_Modifiers(a,p-left.firstLeaf,b,q-right.firstLeaf,k,mods, weight);
	
					//clear alignments that we don't need any more
					pair = SequenceIdPair.makeString(a, k);
					cleanHashes(pair);
					pair = SequenceIdPair.makeString(b, k);
					cleanHashes(pair);
				}
				
				setBlendParams(neighbors.size());
//				float c_weight = neighbors.size()==0? 0 : consistency_weight/neighbors.size();
//				float normalizer = 1/(PairSuboptimalityMatrices.delta*(1+ neighbors.size()==0? 0 : consistency_weight)); // will make the worse-case cost just double the cost of an unmodified edge

				
				ConsistencyModifiers_Pair modpair = mods.modifiers[p-left.firstLeaf][q-right.firstLeaf];

				
				for (int i=0; i<=M; i++) {
					for (int j= (i==0?1:0); j<=N; j++) {//it's not meaningful to modify 0,0														
						
						if (i>0 && j>0) {																												
							modpair.subs[i][j] = calcSub(a,b,i,j, modpair, neighbors.size());
						} else {
							modpair.subs[i][j] = bigBadVal;

						}

						
						if (i>0) {
							modpair.vLambdas[i][j] = calcVLambda(i,j, modpair, neighbors.size());
							modpair.vGammaOpens[i][j] = calcVGammaOpen(i,j, modpair, neighbors.size());
							modpair.vGammaCloses[i][j] = calcVGammaClose(i,j, modpair, neighbors.size());
						} else {
							modpair.vLambdas[i][j] = modpair.vGammaOpens[i][j] = modpair.vGammaCloses[i][j] = bigBadVal;
						}
 
						if (j>0) {
							modpair.hLambdas[i][j] = calcHLambda(i,j, modpair, neighbors.size());
							modpair.hGammaOpens[i][j] = calcHGammaOpen(i,j, modpair, neighbors.size());
							modpair.hGammaCloses[i][j] = calcHGammaClose(i,j, modpair, neighbors.size());
						}	else {
							modpair.hLambdas[i][j] = modpair.hGammaOpens[i][j] = modpair.hGammaCloses[i][j] = bigBadVal;
						}			
						
					}
				}

				
				postProcessAB(a, b, modpair);
			}
		}				
	}


	protected void postProcessAB (int a, int b, ConsistencyModifiers_Pair modpair ) {
		// nothing in most cases
	}
	
}

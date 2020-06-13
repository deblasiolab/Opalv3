package opal.polish;

import java.text.NumberFormat;

import opal.makers.AlignmentMaker;
import opal.tree.TreeNode;
import com.traviswheeler.libs.LogWriter;
import opal.align.Aligner;
import opal.align.Alignment;
import opal.align.ExactCountAligner_Time;
import opal.align.Aligner.AlignmentType;

public class ExhaustiveTwoCutPolisher extends TreePolisher {

	public int passCnt = 0;
	
	public ExhaustiveTwoCutPolisher(TreeNode treeNode, Aligner al, int verbosity) {
		super(treeNode, al, verbosity);
	}
	

	final public void polish ( int reps ) {
		Alignment alignment = root.alignment;
		int K = alignment.K;
		passCnt = 0;
		
		root.assignLeafRanges();
		int low = root.firstLeaf;
		int high = root.lastLeaf;

		int numLeaves = high-low+1;
		int numNodes = 2*numLeaves-1;
	
		TreeNode[] nodeList = new TreeNode[numNodes];
		root.fillNodeList(0, nodeList);
				
		TreeNode node;
		PolishHelper helper;
				
		int[][] fullAlignment = alignment.seqs;
		long prevCost;
	
		if (verbosity>1) {
			LogWriter.stdErrLog("Initial alignment formed"); 
			if (AlignmentMaker.showCost) {
				prevCost = Aligner.calcCost(fullAlignment, alignment.seqIds, aligner.config, alignment.in);
				LogWriter.stdErrLog(", with a cost of " + NumberFormat.getInstance().format( prevCost ));
			}
			LogWriter.stdErrLogln("");
			LogWriter.stdErrLogln("Beginning polishing phase. ");
		}
		int testedCnt = 0;
		int index=0; // this will be a leaf in nodelist

		Aligner fastAligner = aligner;
		Aligner slowAligner = new ExactCountAligner_Time();	
		Aligner origAligner = null;
		boolean stillDoingFast = true;
		while (stillDoingFast) {
			
			if (aligner == slowAligner) { // this will only be true if default aligner is profile, and we've already exhaustively polished the tree with that aligner 
				stillDoingFast = false;
			}
			while (testedCnt < numNodes-1) {
				node = nodeList[index]; 
				 
				int[][] A, B;
				int j,x;
				int sizeA = node.lastLeaf-node.firstLeaf+1;
				int sizeB = K - sizeA;
				
				int[] ABreorder = new int[K];
	
				for (j=0; j<K; j++)  
					ABreorder[j] = -1;
				
	//			split out 2 alignments
					
				A = new int[sizeA][];
				B = new int[sizeB][];
	
				int[] idsA = new int[sizeA];
				int[] idsB = new int[sizeB];				
				
				boolean[] used = new boolean[K];
				for (j=0; j<K; j++) 
					used[j] = false;
				
				int pos;
				for (j=0; j<sizeA; j++){ 
					pos = j + node.firstLeaf;
					ABreorder[j] = pos;
					A[j] = fullAlignment[pos];
					idsA[j] = root.leafOrderFromInput.get(pos);
					used[pos] = true;
				}
				x=0;
				for (j=0; j<K; j++) {
					if (!used[j]) {
						ABreorder[x+sizeA] = j;
						B[x] = fullAlignment[j];
						idsB[x] = root.leafOrderFromInput.get(j);
						x++;
					}
				}
					
				int[][] AB = new int[sizeA + sizeB][];
				int jj = 0;
				for(j=0; j<sizeA; j++) 
					AB[jj++] = A[j];			
				for(j=0; j<sizeB; j++) 
					AB[jj++] = B[j];

				
				int[] idsAB = new int[idsA.length + idsB.length];
				System.arraycopy(idsA, 0, idsAB, 0, idsA.length);
				System.arraycopy(idsB, 0, idsAB, idsA.length, idsB.length);
				
				prevCost = Aligner.calcCost(AB, sizeA, sizeB, idsAB, aligner.config, alignment.in);
				
				
				Alignment alA = Alignment.buildNewAlignment(Alignment.getDegappedCopy(A), idsA, aligner.config, alignment.in);
				Alignment alB = Alignment.buildNewAlignment(Alignment.getDegappedCopy(B), idsB, aligner.config, alignment.in);
				
				boolean usedExactCost;
				if (Polisher.polishAligmentMethod == AlignmentType.mixed &&
						alA.K * alB.K > Aligner.mixedAlignmentCutoff) {
					helper = new PolishHelper(alA,alB,ABreorder,fastAligner);
					usedExactCost = false;
				} else {
					helper = new PolishHelper(alA,alB,ABreorder,aligner);
					usedExactCost = (aligner instanceof ExactCountAligner_Time);
				}
				long cost = helper.getCost();
	
	
				
				if (cost < prevCost &&  !usedExactCost && 
						Polisher.polishAligmentMethod != AlignmentType.profile) {
					//I just did a fast alignment that improved things, 
					//   and exact might be better, (but only if I'm planning on doing them at some point),  so ...
					//LogWriter.stdErrLog(":");
					if ( Polisher.polishAligmentMethod != AlignmentType.mixed || 					
								alA.K * alB.K <= Aligner.mixedAlignmentCutoff)  {					
							helper = new PolishHelper(alA,alB,ABreorder,slowAligner);
							long exactCost = helper.getCost();
							if (exactCost<cost) {
								cost = exactCost;
								usedExactCost = true;
							}
					}
				}

				

				if (cost < prevCost) {					
					root.alignment = helper.getAlignment();
					fullAlignment = root.alignment.seqs;
					prevCost = cost;
					testedCnt = 0;
					if (verbosity>1) {
						if ( usedExactCost )
							LogWriter.stdErrLog("o");
						else
							LogWriter.stdErrLog("*");
					}
					if (origAligner!=null && passCnt < reps-1 ) { // go back to fast aligner, to run through tree again
						stillDoingFast = true;
						aligner = origAligner;
						origAligner = null;
					}
					
				} else {
					if (verbosity>1)
						LogWriter.stdErrLog(".");
				}
				testedCnt++;
				
				index++;
				if (index == numNodes-1) { // that'll be root ... which doesn't have an edge above it to cut
					passCnt++;
					if (verbosity>1) {
						LogWriter.stdErrLog("\nPass " + passCnt + " complete.");
						if (AlignmentMaker.showCost) {
							long c = Aligner.calcCost(fullAlignment, root.alignment.seqIds, aligner.config, alignment.in); // I think those ids are right
							LogWriter.stdErrLog(" Current cost is " + NumberFormat.getInstance().format( c ));
						}
						LogWriter.stdErrLog("\n");
					}					
					if (passCnt == reps || 	(passCnt == reps-1 && origAligner == null)) {
						testedCnt = numNodes-1; //want to bounce out of heuristic alignment, and do one more round of exact. 	
						
					}

					if (passCnt != reps)
						index = 0;
				}
				
				
			}
			// made it exhaustively through the tree using heuristic alignment.
			// Now go through again with exact
			if (passCnt != reps && origAligner == null && Polisher.polishAligmentMethod != AlignmentType.profile ) {
				testedCnt = 0;
				//index = 0;
				origAligner = aligner;
				aligner = slowAligner;
			} else {
				stillDoingFast = false;
				if (passCnt!=reps) passCnt++;
				if (index < numNodes-1) {
					if (verbosity>1) {
						LogWriter.stdErrLogln("\nPass " + passCnt + " complete (other edges have already been tested).");
						if (AlignmentMaker.showCost) { 
							long c = Aligner.calcCost(fullAlignment, root.alignment.seqIds, aligner.config, alignment.in);
							LogWriter.stdErrLog(" Current cost is " + NumberFormat.getInstance().format( c ));
						}
						LogWriter.stdErrLog("\n");
					}
				}
			}
		}
		if (verbosity>1) LogWriter.stdErrLogln("");
	}	
		
}

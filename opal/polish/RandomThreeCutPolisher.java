package opal.polish;

import java.text.NumberFormat;
import java.util.Date;
import java.util.Random;

import opal.makers.AlignmentMaker;
import opal.tree.TreeNode;
import com.traviswheeler.libs.LogWriter;
import opal.align.Aligner;
import opal.align.Alignment;
import opal.align.ExactCountAligner_Time;
import opal.align.Aligner.AlignmentType;

public class RandomThreeCutPolisher extends TreePolisher {

	
	public RandomThreeCutPolisher(TreeNode treeNode, Aligner al, int verbosity) {
		super(treeNode, al, verbosity);
		
		if (!inSeedAssigned) {
			Random rand = new Random();
			seed = rand.nextLong();
		}else{
			seed = inSeed;
		}
		randomGenerator = new Random(seed);
	}
	
	final public void polish ( int reps ) {
		
		
		Alignment alignment = root.alignment;
		int K = alignment.K;
					
		Aligner fastAligner, slowAligner;
		slowAligner = new ExactCountAligner_Time();
		fastAligner = aligner;

		
		root.assignLeafRanges();
		int low = root.firstLeaf;
		int high = root.lastLeaf;

		int numLeaves = high-low+1;
		int numNodes = 2*numLeaves-1;
		TreeNode[] nodeList = new TreeNode[numNodes];
		root.fillNodeList(0, nodeList);
		
		int indexA, indexB;
		TreeNode nodeA, nodeB;
		PolishHelper helper;

		int[][] fullAlignment = alignment.seqs;
		long prevCost = 0;
		if (verbosity>1) {
			LogWriter.stdErrLog("Initial alignment formed"); 
			if (AlignmentMaker.showCost) {
				prevCost = Aligner.calcCost(fullAlignment, alignment.seqIds, aligner.config);
				LogWriter.stdErrLog(", with a cost of " + NumberFormat.getInstance().format( prevCost ));
			}
			LogWriter.stdErrLogln("");
			
			LogWriter.stdErrLogln("Beginning polishing phase, with " + reps + " random three-cut iterations");
			LogWriter.stdErrLogln("Random seed " + seed);
		}
		
		int exactReps = (reps == polishIterations) ? polishIterations_exact : 0;
		int startExhaustive = polishIterations - exactReps + 1;
		int noChangeReps = 0;
//		int startExhaustive = polishIterations - polishIterations_exact + 1;
		Date date1 = new Date();
		for (int i=1; i<=reps; i++) {

			/*
			if (noChangeReps == 100) {
				//bypass the rest of the fast ones
				if (exactReps>0 && i < startExhaustive) {
					i = startExhaustive;
					reps = i + exactReps; 
				} else {
					reps = i;
				}
			}
			*/
			
			if (i == startExhaustive ) 
				aligner = slowAligner;
			
			indexB = indexA = randomGenerator.nextInt( numNodes-1 ); // "-1" so I don't pick root

			while (indexB == indexA ||    
					// can't cut the two edges under the root ... then there's no "group C"
					(nodeList[indexA]==root.leftChild && nodeList[indexB]==root.rightChild) ||  
					(nodeList[indexB]==root.leftChild && nodeList[indexA]==root.rightChild)
					) {
				indexB = randomGenerator.nextInt( numNodes-1 );
			}
			
			nodeA = nodeList[indexA]; 
			nodeB = nodeList[indexB]; 
			
			// I want to have either distinct subtrees, or A inside B
			if (nodeB.firstLeaf >= nodeA.firstLeaf && nodeB.lastLeaf <= nodeA.lastLeaf) {
				nodeB = nodeList[indexA]; 
				nodeA = nodeList[indexB];
			}
			boolean distinctSubTrees = true;
			if (nodeA.firstLeaf >= nodeB.firstLeaf && nodeA.lastLeaf <= nodeB.lastLeaf) {
				distinctSubTrees = false;
			}			
			
			int[][] A, B, C;
			int j,x,y,z;
			int sizeA = nodeA.lastLeaf-nodeA.firstLeaf+1;
			int sizeB = nodeB.lastLeaf-nodeB.firstLeaf+1;
			int sizeC;
			
			int[] ABCreorder = new int[K];
			int[] ACBreorder = new int[K];
			int[] BCAreorder = new int[K];
			
			for (j=0; j<K; j++)  
				ABCreorder[j] = ACBreorder[j] = BCAreorder[j] = -1;

			int[] idsA = new int[sizeA];
			int[] idsB = new int[sizeB];				
			int[] idsC;
			
//			split out 3 alignments
			if ( distinctSubTrees ) { 
				sizeC = K-sizeA-sizeB;

				idsC = new int[sizeC];

				A = new int[sizeA][];
				B = new int[sizeB][];
				C = new int[sizeC][];

				boolean[] used = new boolean[K];
				for (j=0; j<K; j++) 
					used[j] = false;
				
				int pos;
				for (j=0; j<sizeA; j++){ 
					pos = j + nodeA.firstLeaf;
					ABCreorder[j] = ACBreorder[j] = BCAreorder[j+(sizeB+sizeC)] = pos;
					A[j] = fullAlignment[pos];
					idsA[j] = root.leafOrderFromInput.get(pos);
					used[pos] = true;
				}
				for (j=0; j<sizeB; j++) {
					pos = j + nodeB.firstLeaf;
					ABCreorder[j+sizeA] = ACBreorder[j+sizeA+sizeC] = BCAreorder[j] = pos;
					B[j] = fullAlignment[pos];
					idsB[j] = root.leafOrderFromInput.get(pos);
					used[pos] = true;
				}
				x=0;
				for (j=0; j<K; j++) {
					if (!used[j]) {
						ABCreorder[x+sizeA+sizeB] = ACBreorder[x+sizeA] = BCAreorder[x+sizeB] = j;
						C[x] = fullAlignment[j];
						idsC[x] = root.leafOrderFromInput.get(j);
						x++;
					}
				}
				
			} else {// A is in B's subtree
				sizeC = K-sizeB;
				idsC = new int[sizeC];
				
				A = new int[sizeA][];
				B = new int[sizeB-sizeA][];
				C = new int[sizeC][];

				boolean[] used = new boolean[K];
				for (j=0; j<K; j++) 
					used[j] = false;

				
				int pos;
				for (j=0; j<sizeA; j++){ 
					pos = j + nodeA.firstLeaf;
					ABCreorder[j] = ACBreorder[j] = BCAreorder[j+(K-sizeA)] = pos;
					A[j] = fullAlignment[pos];
					idsA[j] = root.leafOrderFromInput.get(pos);
					used[pos] = true;
				}
				x=0;
				for (j=0; j<sizeB; j++) {
					pos = j + nodeB.firstLeaf;
					if (!used[pos]) {
						ABCreorder[x+sizeA] = ACBreorder[x+sizeA+sizeC] = BCAreorder[x] = pos;
						B[x] = fullAlignment[pos];
						idsB[x] = root.leafOrderFromInput.get(pos);
						used[pos] = true;
						x++;
					}
				}
				x=0;
				for (j=0; j<K; j++) {
					if (!used[j]) {
						ABCreorder[x+sizeB] = ACBreorder[x+sizeA] = BCAreorder[x+sizeB-sizeA] = j;
						C[x] = fullAlignment[j];
						idsC[x] = root.leafOrderFromInput.get(j);
						x++;
					}
				}
				sizeB -= sizeA;
			}
			

			//also want copies of AB, AC and BC  from original alignment
			int[][] AB = new int[A.length + B.length][];
			int[][] AC = new int[A.length + C.length][];
			int[][] BC = new int[B.length + C.length][];
			int[] idsAB = new int[A.length + B.length];
			int[] idsAC = new int[A.length + C.length];
			int[] idsBC = new int[B.length + C.length];
			
			x=y=z=0;
			for(j=0; j<A.length; j++) {
				AB[x] = AC[y] = A[j];
				idsAB[x] = idsAC[y] = idsA[j];
				x++;
				y++;
			}
			for(j=0; j<B.length; j++) {
				AB[x] = BC[z] = B[j];
				idsAB[x] = idsBC[z] = idsB[j];
				x++;
				z++;
			}
			for(j=0; j<C.length; j++) {
				AC[y] = BC[z] = C[j];
				idsAC[y] = idsBC[z] = idsC[j];
				y++;
				z++;
			}			
			
			//now make copies (removing gap columns), so we don't break the original			
			Alignment alA = Alignment.buildNewAlignment(Alignment.getDegappedCopy(A), idsA, aligner.config);
			Alignment alB = Alignment.buildNewAlignment(Alignment.getDegappedCopy(B), idsB, aligner.config);
			Alignment alC = Alignment.buildNewAlignment(Alignment.getDegappedCopy(C), idsC, aligner.config);
			Alignment alAB = Alignment.buildNewAlignment(Alignment.getDegappedCopy(AB), idsAB, aligner.config);
			Alignment alAC = Alignment.buildNewAlignment(Alignment.getDegappedCopy(AC), idsAC, aligner.config);
			Alignment alBC = Alignment.buildNewAlignment(Alignment.getDegappedCopy(BC), idsBC, aligner.config);


			
			//boolean usedExactCost = (aligner instanceof ExactCountAligner_Time);
			boolean usedExactCost;		
			long cost, tmpCost;
			PolishHelper bestPolish;
			
			
			long AB_cost = Aligner.calcCost(AB, sizeA, sizeB, idsAB, aligner.config);  
			long AC_cost = Aligner.calcCost(AC, sizeA, sizeC, idsAC, aligner.config);
			long BC_cost = Aligner.calcCost(BC, sizeB, sizeC, idsBC, aligner.config);
			
			cost = prevCost = AB_cost + AC_cost + BC_cost;
 

			//1
			int bigX = alAB.K * alC.K;
			if (i >= startExhaustive && 
					Polisher.polishAligmentMethod == AlignmentType.mixed  &&
					bigX > Aligner.mixedAlignmentCutoff) {
				helper = new PolishHelper(alAB,alC,ABCreorder,fastAligner);
				usedExactCost = false;
			} else {
				// This is what will happen for all runs except the final runs
				//  for mixed method, where large alignments are handled by the fast aligner
				helper = new PolishHelper(alAB,alC,ABCreorder,aligner);
				usedExactCost = (aligner instanceof ExactCountAligner_Time);
			}
			bestPolish = helper; // this will be at least as good as the input alignment
			tmpCost = helper.getCost() + AB_cost;
			if (tmpCost < cost) {	
				cost = tmpCost;
				if ( !usedExactCost  && polishAligmentMethod != AlignmentType.profile ) { 
					//I just did a fast alignment that improved things, 
					//   and exact might be better, (but only if I'm planning on doing them at some point),  so ...
					if ( polishAligmentMethod != AlignmentType.mixed || 					
							bigX <= Aligner.mixedAlignmentCutoff)  {
						helper = new PolishHelper(alAB,alC,ABCreorder,slowAligner);
						long exactCost = helper.getCost();
						if (exactCost<cost) {
							cost = exactCost;
							bestPolish = helper;
							usedExactCost = true;
						}
					}
				} 
			}
			
			
			//2
			bigX = alAC.K * alB.K;
			if (i >= startExhaustive && 
					Polisher.polishAligmentMethod == AlignmentType.mixed  &&
					bigX > Aligner.mixedAlignmentCutoff) {
				helper = new PolishHelper(alAC,alB,ACBreorder,fastAligner);
				usedExactCost = false;
			} else {
				// This is what will happen for all runs except the final runs
				//  for mixed method, where large alignments are handled by the fast aligner
				helper = new PolishHelper(alAC,alB,ACBreorder,aligner);
				usedExactCost = (aligner instanceof ExactCountAligner_Time);
			}
			tmpCost = helper.getCost() + AC_cost;
				
			if (tmpCost < cost) {	
				cost = tmpCost;
				bestPolish = helper;
				if ( !usedExactCost  && polishAligmentMethod != AlignmentType.profile ) { 
					if ( polishAligmentMethod != AlignmentType.mixed || 					
							bigX <= Aligner.mixedAlignmentCutoff)  {
						helper = new PolishHelper(alAC,alB,ACBreorder,slowAligner);
						long exactCost = helper.getCost();
						if (exactCost<cost) {
							cost = exactCost;
							bestPolish = helper;
							usedExactCost = true;
						}
					}
				}
			}
			
			//3
			bigX = alBC.K * alA.K;
			if (i >= startExhaustive && 
					Polisher.polishAligmentMethod == AlignmentType.mixed  &&
					bigX > Aligner.mixedAlignmentCutoff) {
				helper = new PolishHelper(alBC,alA,BCAreorder,fastAligner);
				usedExactCost = false;
			} else {
				// This is what will happen for all runs except the final runs
				//  for mixed method, where large alignments are handled by the fast aligner
				helper = new PolishHelper(alBC,alA,BCAreorder,aligner);
				usedExactCost = (aligner instanceof ExactCountAligner_Time);
			}
			tmpCost = helper.getCost() + BC_cost;

			if (tmpCost < cost) {	
				cost = tmpCost;
				bestPolish = helper;
				if ( !usedExactCost  && polishAligmentMethod != AlignmentType.profile ) { 
					if ( polishAligmentMethod != AlignmentType.mixed || 					
							bigX <= Aligner.mixedAlignmentCutoff)  {
						helper = new PolishHelper(alBC,alA,BCAreorder,slowAligner);
						long exactCost = helper.getCost();
						if (exactCost<cost) {
							cost = exactCost;
							bestPolish = helper;
							usedExactCost = true;
						}
					}
				}
			}

			
			//4
			bigX = Math.max(alA.K * alC.K, AC.length * alB.K);			
			if (i >= startExhaustive && 
					Polisher.polishAligmentMethod == AlignmentType.mixed  &&
					bigX > Aligner.mixedAlignmentCutoff) {
				helper = new PolishHelper(alA,alC,alB,ACBreorder,fastAligner);
				usedExactCost = false;
			} else {
				helper = new PolishHelper(alA,alC,alB,ACBreorder,aligner);
				usedExactCost = (aligner instanceof ExactCountAligner_Time);
			}
			tmpCost = helper.getCost();			
			
			if (tmpCost < cost) {	
				cost = tmpCost;
				bestPolish = helper;
				if ( !usedExactCost  && polishAligmentMethod != AlignmentType.profile ) { 
					if ( polishAligmentMethod != AlignmentType.mixed || 					
							bigX <= Aligner.mixedAlignmentCutoff)  {
						helper = new PolishHelper(alA,alC,alB,ACBreorder,slowAligner);
						long exactCost = helper.getCost();
						if (exactCost<cost) {
							cost = exactCost;
							bestPolish = helper;
							usedExactCost = true;
						}
					}
				}
			}

			//5
			bigX = Math.max(alA.K * alB.K, AB.length * alC.K);			
			if (i >= startExhaustive && 
					Polisher.polishAligmentMethod == AlignmentType.mixed  &&
					bigX > Aligner.mixedAlignmentCutoff) {
				helper = new PolishHelper(alA,alB,alC,ABCreorder,fastAligner);
				usedExactCost = false;
			} else {
				helper = new PolishHelper(alA,alB,alC,ABCreorder,aligner);
				usedExactCost = (aligner instanceof ExactCountAligner_Time);
			}
			tmpCost = helper.getCost();			

			if (tmpCost < cost) {	
				cost = tmpCost;
				bestPolish = helper;
				if ( !usedExactCost  && polishAligmentMethod != AlignmentType.profile ) { 
					if ( polishAligmentMethod != AlignmentType.mixed || 					
							bigX <= Aligner.mixedAlignmentCutoff)  {
						helper = new PolishHelper(alA,alB,alC,ABCreorder,slowAligner);
						long exactCost = helper.getCost();
						if (exactCost<cost) {
							cost = exactCost;
							bestPolish = helper;
							usedExactCost = true;
						}
					}
				}
			}
			
			
			//6
			bigX = Math.max(alB.K * alC.K, BC.length * alA.K);			
			if (i >= startExhaustive && 
					Polisher.polishAligmentMethod == AlignmentType.mixed  &&
					bigX > Aligner.mixedAlignmentCutoff) {
				helper = new PolishHelper(alB,alC,alA,BCAreorder,fastAligner);
				usedExactCost = false;
			} else {
				helper = new PolishHelper(alB,alC,alA,BCAreorder,aligner);
				usedExactCost = (aligner instanceof ExactCountAligner_Time);
			}
			tmpCost = helper.getCost();			
			
			if (tmpCost < cost) {	
				cost = tmpCost;
				bestPolish = helper;
				if ( !usedExactCost  && polishAligmentMethod != AlignmentType.profile ) { 
					if ( polishAligmentMethod != AlignmentType.mixed || 					
							bigX <= Aligner.mixedAlignmentCutoff)  {
						helper = new PolishHelper(alB,alC,alA,BCAreorder,slowAligner);
						long exactCost = helper.getCost();
						if (exactCost<cost) {
							cost = exactCost;
							bestPolish = helper;
							usedExactCost = true;
						}
					}
				}
			}
			

			if (!keepGoing) {
				try {  // waste time while the printout does it's job
					Thread.sleep(100000);
				} catch (InterruptedException e) {	}  
			}

			if (cost < prevCost) {
				root.alignment = bestPolish.getAlignment();
				fullAlignment = root.alignment.seqs;
				prevCost = cost;
				if (verbosity>1) {
					if ( usedExactCost )
						LogWriter.stdErrLog("o");
					else
						LogWriter.stdErrLog("*");
				}
				noChangeReps = 0;
			} else {
				noChangeReps++;
				if (verbosity>1)
					LogWriter.stdErrLog(".");
			}

			if (i%20 == 0) { 
				Date date2 = new Date();
				long diff2 = date2.getTime() - date1.getTime();
				if (verbosity>0)
					LogWriter.stdErrLogln(" Round " + i + " complete. (" + diff2 + " ms)");
				date1 = new Date();
			}
			
		}
		if (verbosity>1) LogWriter.stdErrLogln("");

	}	
		
}

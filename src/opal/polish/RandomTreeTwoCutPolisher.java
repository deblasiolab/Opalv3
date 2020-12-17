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

public class RandomTreeTwoCutPolisher extends TreePolisher {

	public RandomTreeTwoCutPolisher(TreeNode treeNode, Aligner al, int verbosity) {
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
		slowAligner.config = fastAligner.config;
		
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
			
			LogWriter.stdErrLogln("Beginning polishing phase, with " + reps + " random two-cut iterations");
			LogWriter.stdErrLogln("Random seed " + seed);
		}
		
		int exactReps = (reps == polishIterations) ? polishIterations_exact : 0;
		int startExhaustive = polishIterations - exactReps + 1;
		Date date1 = new Date();
		int noChangeReps = 0;
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
			
			int index = randomGenerator.nextInt( numNodes-1 );
		  System.err.println(aligner.config.toString() + "Random Index: " + index);	
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

			
			if (i >= startExhaustive && 
					Polisher.polishAligmentMethod == AlignmentType.mixed  &&
					alA.K * alB.K > Aligner.mixedAlignmentCutoff) {
				helper = new PolishHelper(alA,alB,ABreorder,fastAligner);
				usedExactCost = false;
			} else {
				// This is what will happen for all runs except the final runs
				//  for mixed method, where large alignments are handled by the fast aligner
				helper = new PolishHelper(alA,alB,ABreorder,aligner);
				usedExactCost = (aligner instanceof ExactCountAligner_Time);
			}
			long cost = helper.getCost();
						

			if (cost < prevCost &&  !usedExactCost  &&
					polishAligmentMethod != AlignmentType.profile ) { 
					//exactReps > 0 ) {
				//I just did a fast alignment that improved things, 
				//   and exact might be better, (but only if I'm planning on doing them at some point),  so ...
				if ( polishAligmentMethod != AlignmentType.mixed || 					
						alA.K * alB.K <= Aligner.mixedAlignmentCutoff)  {
					//LogWriter.stdErrLog(":");
					helper = new PolishHelper(alA,alB,ABreorder,slowAligner);
					long exactCost = helper.getCost();
					if (exactCost<cost) {
						cost = exactCost;
						usedExactCost = true;
					}
				}
			}
			

			if (!keepGoing) {
				try {  // waste time while the printout does it's job
					Thread.sleep(100000);
				} catch (InterruptedException e) {	}  
			}

			if (cost < prevCost) {
				root.alignment = helper.getAlignment();
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
				if (verbosity>1) 
					LogWriter.stdErrLogln(" Round " + i + " complete. (" + diff2 + " ms)");
				date1 = new Date();
			}
						
		}
		if (verbosity>1) LogWriter.stdErrLogln("");

	}	
		
}

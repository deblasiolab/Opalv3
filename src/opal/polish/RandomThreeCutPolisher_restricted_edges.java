package opal.polish;

import java.util.Random;
import opal.tree.TreeNode;
import com.traviswheeler.libs.LogWriter;
import opal.align.Aligner;
import opal.align.Alignment;
import opal.align.ExactCountAligner_Time;
import opal.align.ProfileAligner;
import opal.align.Aligner.AlignmentType;

public class RandomThreeCutPolisher_restricted_edges extends TreePolisher {

	private Random randomGenerator;
	static long seed;
	static boolean seedAssigned = false;
	
	public RandomThreeCutPolisher_restricted_edges(TreeNode treeNode, Aligner al, int verbosity) {
		super(treeNode, al, verbosity);
		
		if (!seedAssigned) {
			Random rand = new Random();
			seed = rand.nextLong();
			seedAssigned = true;
		}
		randomGenerator = new Random(seed);
	}
	
	public static void setRandomSeed (long s) {
		seedAssigned = true;
		seed = s;
	}

	public long getRandomSeed () {
		return seed;
	}

	final public void polish ( int reps ) {
		Alignment alignment = root.alignment;
		int K = alignment.K;
				
		Aligner fastAligner = new ProfileAligner(aligner);
		fastAligner.setPessimistic(true);

		Aligner slowAligner;
		if (aligner instanceof ExactCountAligner_Time) 
			slowAligner = aligner;
		else 
			slowAligner = new ExactCountAligner_Time();

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
		long prevCost = Aligner.calcCost(fullAlignment, alignment.seqIds, aligner.config, alignment.in);
		int startExhaustive = polishIterations - polishIterations_exact + 1;

		if (verbosity>1) {
			LogWriter.stdErrLogln("Initial alignment formed, with a cost of " + prevCost);
			LogWriter.stdErrLogln("Beginning polishing phase, with " + reps + " random three-cut iterations");
			LogWriter.stdErrLogln("Random seed " + seed);
		}
		for (int i=1; i<=reps; i++) {

			if (i == startExhaustive && aligner instanceof ProfileAligner) {
				aligner = slowAligner;
			}
			
			while (true) { //keep going 'till I've got nodes that are part of separate subtrees, and not siblings
				indexA = randomGenerator.nextInt( numNodes-1 ); // "-1" so I don't pick root
				indexB = randomGenerator.nextInt( numNodes-1 ); // "-1" so I don't pick root
				
				nodeA = nodeList[indexA]; 
				nodeB = nodeList[indexB]; 

				if ( (nodeA.firstLeaf >= nodeB.firstLeaf && nodeA.lastLeaf <= nodeB.lastLeaf) || /*A is in B's subtree*/
						(nodeB.firstLeaf >= nodeA.firstLeaf && nodeB.lastLeaf <= nodeA.lastLeaf) || /* B is in As subtree */
						nodeA.parent == nodeB.parent
					)  {
					
				} else {
					break;
				}
			}			
			
			
			int[][] A, B, C;
			int j,x,y,z;
			int sizeA = nodeA.lastLeaf-nodeA.firstLeaf+1;
			int sizeB = nodeB.lastLeaf-nodeB.firstLeaf+1;
			int sizeC = K-sizeA-sizeB;
			
			int[] ABCreorder = new int[K];
			int[] ACBreorder = new int[K];
			int[] BCAreorder = new int[K];

			int[] idsA = new int[sizeA];
			int[] idsB = new int[sizeB];				
			int[] idsC = new int[sizeC];

			for (j=0; j<K; j++)  
				ABCreorder[j] = ACBreorder[j] = BCAreorder[j] = -1;

			
//			split out 3 alignments
			
			sizeC = K-sizeA-sizeB;
			
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
			Alignment alA = Alignment.buildNewAlignment(Alignment.getDegappedCopy(A), idsA, aligner.config, alignment.in);
			Alignment alB = Alignment.buildNewAlignment(Alignment.getDegappedCopy(B), idsB, aligner.config, alignment.in);
			Alignment alC = Alignment.buildNewAlignment(Alignment.getDegappedCopy(C), idsC, aligner.config, alignment.in);
			Alignment alAB = Alignment.buildNewAlignment(Alignment.getDegappedCopy(AB), idsAB, aligner.config, alignment.in);
			Alignment alAC = Alignment.buildNewAlignment(Alignment.getDegappedCopy(AC), idsAC, aligner.config, alignment.in);
			Alignment alBC = Alignment.buildNewAlignment(Alignment.getDegappedCopy(BC), idsBC, aligner.config, alignment.in);

			boolean usedExactCost = (aligner instanceof ExactCountAligner_Time);

			long cost, tmpCost;
			PolishHelper bestPolish;
			bestPolish = helper = new PolishHelper(alA,alB,alC,ABCreorder,aligner);
			cost = helper.getCost();
			if (cost < prevCost && Aligner.alignmentMethod == AlignmentType.exact  && aligner instanceof ProfileAligner) {
				//I just did a fast alignment that improved things, and exact might be better, so ...
				bestPolish = helper = new PolishHelper(alA,alB,alC,ABCreorder,slowAligner);
				long exactCost = helper.getCost();
				if (exactCost<cost) {
					cost = exactCost;
					usedExactCost = true;
				} else {
					usedExactCost = false;
				}
			}

			
			helper = new PolishHelper(alAB,alC,ABCreorder,aligner);
			tmpCost = helper.getCost();
			if (tmpCost < prevCost && Aligner.alignmentMethod == AlignmentType.exact && aligner instanceof ProfileAligner) {
				helper = new PolishHelper(alAB,alC,ABCreorder,slowAligner);
				long exactCost = helper.getCost();
				if (exactCost<tmpCost ) {
					tmpCost = exactCost;
					if (tmpCost<cost ) usedExactCost = true;
				} else {
					if (tmpCost<cost ) usedExactCost = false;
				}				
			}
			if (tmpCost < cost) {
				cost = tmpCost;
				bestPolish = helper; 
			}
				
			helper = new PolishHelper(alA,alC,alB,ACBreorder,aligner);
			tmpCost = helper.getCost();
			if (tmpCost < prevCost && Aligner.alignmentMethod == AlignmentType.exact && aligner instanceof ProfileAligner) {
				helper = new PolishHelper(alA,alC,alB,ACBreorder,slowAligner);
				long exactCost = helper.getCost();
				if (exactCost<tmpCost ) {
					tmpCost = exactCost;
					if (tmpCost<cost ) usedExactCost = true;
				} else {
					if (tmpCost<cost ) usedExactCost = false;
				}
			}
			if (tmpCost < cost) {
				cost = tmpCost;
				bestPolish = helper; 
			}
			
			helper = new PolishHelper(alAC,alB,ACBreorder,aligner);
			tmpCost = helper.getCost();			
			if (tmpCost < prevCost && Aligner.alignmentMethod == AlignmentType.exact && aligner instanceof ProfileAligner) {
				helper = new PolishHelper(alAC,alB,ACBreorder,slowAligner);
				long exactCost = helper.getCost();
				if (exactCost<tmpCost ) {
					tmpCost = exactCost;
					if (tmpCost<cost ) usedExactCost = true;
				} else {
					if (tmpCost<cost ) usedExactCost = false;
				}					
			}
			if (tmpCost < cost) {
				cost = tmpCost;
				bestPolish = helper; 
			}

		
			helper = new PolishHelper(alB,alC,alA,BCAreorder,aligner);
			tmpCost = helper.getCost();
			if (tmpCost < prevCost && Aligner.alignmentMethod == AlignmentType.exact && aligner instanceof ProfileAligner) {
				helper = new PolishHelper(alB,alC,alA,BCAreorder,slowAligner);
				long exactCost = helper.getCost();
				if (exactCost<tmpCost ) {
					tmpCost = exactCost;
					if (tmpCost<cost ) usedExactCost = true;
				} else {
					if (tmpCost<cost ) usedExactCost = false;
				}
			}
			if (tmpCost < cost) {
				cost = tmpCost;
				bestPolish = helper; 
			}

			helper = new PolishHelper(alBC,alA,BCAreorder,aligner);
			tmpCost = helper.getCost();			
			if (tmpCost < prevCost && Aligner.alignmentMethod == AlignmentType.exact && aligner instanceof ProfileAligner) {
				helper = new PolishHelper(alBC,alA,BCAreorder,slowAligner);
				long exactCost = helper.getCost();
				if (exactCost<tmpCost ) {
					tmpCost = exactCost;
					if (tmpCost<cost ) usedExactCost = true;
				} else {
					if (tmpCost<cost ) usedExactCost = false;
				}
			}
			if (tmpCost < cost) {
				cost = tmpCost;
				bestPolish = helper; 
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
			} else {
				if (verbosity>1)
					LogWriter.stdErrLog(".");
			}
			if (verbosity>1) {
				if (i%20 == 0)
					LogWriter.stdErrLogln(" Round " + i + " complete. Current cost is " + prevCost);
			}
		}
		if (verbosity>1) LogWriter.stdErrLogln("");

	}	
		
}

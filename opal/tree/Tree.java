package opal.tree;

import java.util.ArrayList;
import java.util.Date;

import opal.IO.Configuration;
import opal.makers.AlignmentMaker;
import opal.makers.AlignmentMaker_SingleSequences;
import com.traviswheeler.libs.LogWriter;
import opal.align.Aligner;
import opal.align.Alignment;
import opal.align.PairAligner;
import opal.align.ProfileAligner;
import opal.align.Aligner.AlignmentType;
import opal.exceptions.GenericOpalException;

public abstract class Tree {

	public static enum TreeType {mst, entered};
	public static enum DistanceType {normcost, pctid, kmer, kmer_normcost};
	public static int iterations = -1;
	
	public static TreeType treeType = TreeType.mst;
	public static DistanceType distanceType = null; //DistanceType.kmer_normcost; //DistanceType.normcost;
	public static boolean justPWDists = false;
	
	protected ArrayList<TreeNode> roots = new ArrayList<TreeNode>();
	protected TreeNode nextA, nextB;
	protected Aligner defaultAligner, pairAligner, fastAligner;
	private TreeNode latestNode;
	protected TreeNode[] nodesToMerge;
	protected TreeNode[] leaves; 
	private int mergeIndex;
	protected float distances[][];
	private Configuration config;

	protected Distance dist;
	
	Date start;
	public boolean keepGoing = true;
	
	int verbosity = 1;
	
//	public Tree (int[][] seqs, Distance dist, Aligner al) {
	public Tree (Alignment[] alignments, Aligner al, Distance distance, int currIteration, int verbosity, Configuration c) {	
		dist = distance;
		defaultAligner = al;
		config = c;
		
		pairAligner = new PairAligner(al);		
		this.verbosity = verbosity;
		
		if (Aligner.alignmentMethod == AlignmentType.mixed) {
			fastAligner = new ProfileAligner();
		}
		int K = alignments.length;
		if (currIteration < iterations || AlignmentMaker.initAlignmentProvided)
			K = alignments[0].K;
			
		nodesToMerge = new TreeNode[K-1];
		leaves = new TreeNode[K];

		initializeStructures(K, alignments, currIteration);
		defaultAligner.preprocessTree(this, distances);

		if (justPWDists){ 
			boolean self = (distanceType==DistanceType.normcost);
			for (int p = 0; p<distances.length-(self?0:1); p++) {
				for (int q = p+(self?0:1); q<distances.length; q++) {
					LogWriter.stdOutLogln(p+","+q+ ": " + distances[p][q]);
				}
			}
			return;
		}

		distances = null;
		start = new Date();
	}
	
	protected void initializeStructures (int K, Alignment[] alignments, int currIteration) {
		for (int i=0; i<K; i++) {
			TreeNode node = new MST_TreeNode(K,i);
			if (currIteration<iterations || AlignmentMaker.initAlignmentProvided) {
				int[][] tmp = new int[1][];
				tmp[0] = alignments[0].seqs[i];
				Alignment tmpAl;
				tmpAl = Alignment.buildNewAlignment(Alignment.getDegappedCopy(tmp)[0], i, config);
				node.setAlignment(tmpAl);
			} else { 
				node.setAlignment(alignments[i]);
			}
			roots.add(node);
		}
		initializePriorityDS();
		
//		Distance dist = getDistance(pairAligner);
		distances = new float[K][K];
		if (verbosity> 1) {
			LogWriter.stdErrLogln("While performing initial pairwise alignments, Opal will print a dot for every sequence.");
			LogWriter.stdErrLogln("When " + K + " dots are printed, Opal will move on to constructing the multiple alignment.");
		}
		for (int i=0; i<K; i++) {
			TreeNode node_i = roots.get(i);
			if (verbosity > 1) {	
				LogWriter.stdErrLog(".");
				if ( (i+1)%20 == 0) LogWriter.stdErrLogln("  " + (i+1) + " complete");
			}
			
			for (int j=i; j<K; j++) {
				TreeNode node_j = roots.get(j);
//				LogWriter.stdErrLog(i + "," + j + ": " );

				float d;
				if (currIteration == iterations && !AlignmentMaker.initAlignmentProvided) {
					d = dist.calcDistance( node_i.alignment, node_j.alignment);
				} else {
					int[][] tmp = new int[2][];
					tmp[0]= alignments[0].seqs[i];
					tmp[1]= alignments[0].seqs[j];
					int[] ids = {i,j};
					Alignment tmpAl; 
					tmpAl = Alignment.buildNewAlignment(tmp, ids, config);
					d = ((NormCostDistance)dist).calcDistance( tmpAl);
				}
				int a = node_i.leafOrderFromInput.get(0);
				int b = node_j.leafOrderFromInput.get(0);
				distances[a][b] = distances[b][a] = d;

				if (i!=j) {
					TreeNodePair pair = new TreeNodePair(node_i, node_j, d);
		 			addToPriorityDS(pair);
					
					node_i.setDistance(j,d);
					node_j.setDistance(i,d);
				}
				
			}
		}
		if (verbosity > 1) LogWriter.stdErrLogln("\n Pairwise alignments complete.");
		
		finishInitialization();
	}	
	
	protected void finishInitialization() { // should be overridden for method like DAD
		while (roots.size() > 1) {
			clusterNext();
		}
		killPriorityDS();

		TreeNode root = roots.get(0);
		root.fillInternalNodeList(0, nodesToMerge);
	}
	
	public TreeNode mergeNext () {
		
		if (!keepGoing) {
			try {  // waste time while the printout does it's job ... it's about to quit
				Thread.sleep(100000);
			} catch (InterruptedException e) {	}  
		}

		
		TreeNode newNode = getNext();
		if (newNode == null ) return null;
		
		TreeNode left = newNode.leftChild;
		TreeNode right = newNode.rightChild;
		Aligner aligner;
		
		//align nextA and nextB
		Alignment alA = left.alignment;
		Alignment alB = right.alignment;
		if ( alA.K==1 && alB.K==1 && !AlignmentMaker_SingleSequences.consistency  ) 
			aligner = pairAligner;		
		else if (!AlignmentMaker_SingleSequences.consistency && 
				Aligner.alignmentMethod == AlignmentType.mixed && 
				alA.K*alB.K > Aligner.mixedAlignmentCutoff)
			aligner = fastAligner;
		else 	
			aligner = defaultAligner;

	
		aligner.setAlignments(alA, alB);
		aligner.setNode(newNode);
		aligner.config = config;
		
		if (verbosity > 2) {
			LogWriter.stdErrLog("merging nodes: (" );
			int cnt=0;
			for (int id: alA.seqIds) { 
				LogWriter.stdErrLog(""+id);
				if (++cnt<alA.K) LogWriter.stdErrLog(",");
			}
			LogWriter.stdErrLog(  ") to (" );
			cnt=0;
			for (int id: alB.seqIds) { 
				LogWriter.stdErrLog(""+id);
				if (++cnt<alB.K) LogWriter.stdErrLog(",");
			}
			LogWriter.stdErrLog(  ") " );

			
	/*		
			
			LogWriter.stdErrLog("merging nodes: (" + left.firstLeaf);
			if (left.firstLeaf != left.lastLeaf) LogWriter.stdErrLog( "-" + left.lastLeaf);
			LogWriter.stdErrLog(  ") to " );
			
			LogWriter.stdErrLog("(" + right.firstLeaf);
			if (right.firstLeaf != right.lastLeaf) LogWriter.stdErrLog( "-" + right.lastLeaf);
			LogWriter.stdErrLog( ")  " );
			LogWriter.stdErrLog("(lengths: " + alA.M + " and " + alB.M + ")   " );
	*/
	 		Runtime runtime = Runtime.getRuntime();
			long maxMemory = runtime.maxMemory();
			long allocatedMemory = runtime.totalMemory();
			long freeMemory = runtime.freeMemory();
			LogWriter.stdErrLog("(freemem: " + (freeMemory + (maxMemory - allocatedMemory)) / (1024*1024) + " MB),    ");
		}
		
		//Date date1 = new Date();		
		aligner.align();
		//Date date2 = new Date();
		//long diff2 = date2.getTime() - date1.getTime();		
		aligner.getAlignment();

		if (!keepGoing) {
			try {  // waste time while the printout does it's job
				Thread.sleep(100000);
			} catch (InterruptedException e) {	}  
		}

		if (verbosity > 2) {
			Date now = new Date();
			long diff = now.getTime() - start.getTime();
			LogWriter.stdErrLogln("(" + diff + " ms)");
			start = now;
		}

		
		left.setAlignment(null);
		right.setAlignment(null);
		
		aligner.cleanup();
		
		return newNode;
	}

	
	public void prepareToMerge () {
		mergeIndex = 0;
	}

	public int mergesRemaining () {
		return nodesToMerge.length - mergeIndex ;
	}

	
	protected TreeNode getNext () {
		//	For most cases, there's already a tree, so just pick next node
		//	For method like DAD, override this
		
		if (mergeIndex >= nodesToMerge.length ) return null;
		TreeNode next = nodesToMerge[mergeIndex++];
		return next;
	}
	
	public TreeNode[] getNodesToMerge () {
		return nodesToMerge;
	}
	
	protected TreeNode clusterNext () {
		pickNextPair(); //implemented in subclasses, sets node pointers nextA, nextB
		
		//make new node (latestNode), set alignment to result
		latestNode = new MST_TreeNode(roots.size());
				
		//calculate distance from new node to all others, adding them to heap (and nodes nextChoiceDS_pointers list )
//		 note: this must happen before clearOldValues is called on nextA/B
		latestNode.setChildren(nextA, nextB); 
		roots.add(latestNode);
		updateDistances(latestNode);

		
		//remove alignments from nextA and nextB, and remove nextA and nextB from root list
		nextA.clearOldValues();
		nextB.clearOldValues();
		roots.remove(nextA);
		roots.remove(nextB);

		return latestNode;
	}

	public ArrayList<TreeNode> getRoots () {
		return roots;
	}
		
	public void buildTreeString (StringBuffer sb, String[] names) {
		//TreeNode.resetLeafCounter();
		roots.get(0).buildTreeString(sb, names, roots.get(0).leafOrderFromInput);
	}
	
	protected abstract void initializePriorityDS() ;
	protected abstract void killPriorityDS() ;
	protected abstract void addToPriorityDS (TreeNodePair pair); 
	protected abstract void pickNextPair();
	protected abstract void updateDistances(TreeNode node);

	
}

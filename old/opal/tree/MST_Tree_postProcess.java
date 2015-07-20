package opal.tree;

import opal.align.Aligner;
import opal.align.Alignment;
import opal.structures.PairingHeap;

public class MST_Tree_postProcess extends Tree {

	PairingHeap heap;
	
	public MST_Tree_postProcess(Alignment[] alignments, Aligner al, Distance distance, int currIteration, int verbosity) {
		super(alignments, al, distance, currIteration, verbosity);
	}

	final protected void initializePriorityDS() {
		heap = new PairingHeap();
	}
	
	final protected void finishInitialization () {	
	}
	
	final protected void killPriorityDS() {
		
	}
	
	final protected void addToPriorityDS(TreeNodePair pair) {
		heap.insert(pair);
	}
	
	final protected void pickNextPair() {
		TreeNodePair nodePair = (TreeNodePair)(heap.deleteMin());
		
		//In this implementation, I don't bother cleaning out old heap entries when the
		// nodes they relate to have been merged into some new node. Instead,
		// I just loop through repeated calls to deleteMin until I find a "live"
		// TreeNodePair. 
		while (nodePair.A.alignment == null || nodePair.B.alignment == null) {
			nodePair = (TreeNodePair)(heap.deleteMin());
		}
		nextA = nodePair.A;
		nextB = nodePair.B;
	}
	
	final protected void updateDistances(TreeNode newNode) {
		//sets distance from new node to all existing nodes ... including the two merged ones.
		int x = -1, y= -1;
		int K = roots.size();  // note this list already includes the new node.
		for (int i=0; i<K-1; i++) {
			TreeNode node = (TreeNode)roots.get(i);

			float d = nextA.getDistance(i);
			float d2 = nextB.getDistance(i);
			if (d2<d)
				d = d2;
			
			// put distance to i into new node's array
			newNode.setDistance(i, d);

			// put distance to newNode at end of i's array
			// this position corresponds to the position in "roots" where the new node was just placed
			node.setDistance(K-1, d);

			//add distance to the priority queue
			TreeNodePair pair = new TreeNodePair(node, newNode, d);
			addToPriorityDS(pair);

			//capture indexes for nodes to be deleted
			if (node == nextA) 
				x = i;
			else if (node == nextB) 
				y = i;
		}
		if (x<y) { // order them so the removal below doesn't cause a problem
			int z = y;
			y=x;
			x=z;
		}
		
		//For all nodes, including the new one, pull out the distance values for the nodes to be removed.  
		//This keeps the distance lists in line with roots (where the nodes will be removed in a moment.
		for (int i=0; i<roots.size(); i++) {
			if (i==x || i==y) continue;
			TreeNode node = (TreeNode)roots.get(i);
			node.removeDistances(x,y);
		}
	}

	
/*	final protected void cleanPriorityDS( TreeNode node) {
		ArrayList pointers = node.clearDSTPointers();
		for (int i=0; i<pointers.size(); i++ ) {
			HeapNode hn = (HeapNode)pointers.get(i);
			if (hn.prev != hn) { // this is used as a marker for having already been deleted from the heap 
				heap.delete(hn); 
				hn.prev = hn; 
			}
		}
	}
*/
}

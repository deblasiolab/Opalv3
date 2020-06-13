package opal.align;

import opal.tree.TreeNode;

public class ConsistencyModifiers_AllPairs {

	public ConsistencyModifiers_Pair[][] modifiers;
	
	public ConsistencyModifiers_AllPairs (ConsistencyModifiers_Pair modpair) {
		
		modifiers = new ConsistencyModifiers_Pair[1][1];	
		modifiers[0][0] = modpair;
	}

	
	public ConsistencyModifiers_AllPairs (TreeNode node, PairwiseAlignmentsContainer cont) {
		TreeNode left = node.leftChild;
		TreeNode right = node.rightChild;
		int K = left.lastLeaf - left.firstLeaf + 1;
		int L = right.lastLeaf - right.firstLeaf + 1;
		
		modifiers = new ConsistencyModifiers_Pair[K][L];
		
		
		for (int i=0; i<K; i++) {
			int a = cont.getInputSeqID(i+left.firstLeaf);
			int M = cont.origSeqs[a].length;
			for (int j=0; j<L; j++) {
				int b = cont.getInputSeqID(j+right.firstLeaf); 
				int N = cont.origSeqs[b].length;
				modifiers[i][j] = new ConsistencyModifiers_Pair(M,N);
			}							
		}

		
		cont.calcAllPairModifiers(node, this);

	}

	
}

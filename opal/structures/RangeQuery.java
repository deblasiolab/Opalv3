package opal.structures;

import opal.exceptions.GenericOpalException;

public class RangeQuery {

	int[] leafPointers;
	//double[] maxRangeTree;
	int[] maxRangeTree;
	
	//public RangeQuery(double[] values) {
	public RangeQuery(int[] values) {
		int leafCount = values.length;
		int nodeCount = 2*leafCount -1;
		leafPointers = new int[leafCount];
		//maxRangeTree = new double [nodeCount];
		maxRangeTree = new int[nodeCount];
		
		assignLeafValues ( values);

		pushUpMaxValues();
	}

//	private void assignLeafValues ( double[] values ) {
	private void assignLeafValues ( int[] values ) {
		int leafCount = values.length;
		int nodeCount = 2*values.length -1;

		int fullDepth = (new Double (Math.log(nodeCount)/Math.log(2))).intValue()-1;
		int fullDepthNodes = (new Double(Math.pow((double)2,(double)fullDepth))).intValue(); // how many leaves there would be in a full binary tree of this depth
		int fullDepthInternalNodes = leafCount - fullDepthNodes ; // if this != 0, then this many sequences will have leaves on the next (unfull) level
		int nextLevelLeaves = 2 * fullDepthInternalNodes; 

		int i;
		int startOfFullLevel = fullDepthNodes-1;
		int startOfNextLevel = startOfFullLevel+fullDepthNodes;
		int pos;
		for (i=0; i<nextLevelLeaves; i++) {
			pos = i+startOfNextLevel;
			maxRangeTree[pos] = values[i];
			leafPointers[i] = pos;
		}
		for (i=0; i<leafCount-nextLevelLeaves; i++) {
			pos = i + fullDepthInternalNodes + startOfFullLevel;
			maxRangeTree[pos] = values[i + nextLevelLeaves];
			leafPointers[i + nextLevelLeaves] = pos;			
		}
	}
	
	
	private void pushUpMaxValues( ){
		int emptyNodes = maxRangeTree.length - leafPointers.length;
		for (int i = emptyNodes - 1; i>=0; i--) {
			maxRangeTree[i] = Math.max(maxRangeTree[2*i+1], maxRangeTree[2*i+2]);  // left and right child 
		}
	}

	public int getMax (int leftIndex, int rightIndex) {
		/*if (rightIndex == -1 ) { // special case at the boundaries of the DP table
			return 0;
		}*/
		if (leftIndex > rightIndex || leftIndex < 0 || rightIndex >= leafPointers.length) { 
			throw new GenericOpalException("Invalid request for max value in range " + leftIndex + " to " + rightIndex);
		}
//		double max = -1;
		int max = -1;
		if (rightIndex - leftIndex < 8 ) {
			//small range, just scan through values
			for (int i=leftIndex; i<=rightIndex; i++ ) {
				if (maxRangeTree[leafPointers[i]] > max)
					max = maxRangeTree[leafPointers[i]];
			}
		} else { 
			//large range, so use log(n) tree traversal method
			int leftNode = leafPointers[leftIndex];
			int rightNode = leafPointers[rightIndex];
//			double leftMax =  maxRangeTree[leftNode];
//			double rightMax =  maxRangeTree[rightNode];
			int leftMax =  maxRangeTree[leftNode];
			int rightMax =  maxRangeTree[rightNode];
			
			if (leftNode > rightNode) { 
				//left node is on lower node of the tree. 
				//Push it up one before starting the next lockstep move.
				if (leftNode%2 == 1) {//odd index
					//Left child of it's parent; possibly set parent's max = other sibling's max
					if (maxRangeTree[leftNode+1] /*right sibling*/ > leftMax)
						leftMax = maxRangeTree[leftNode+1];
				} // else, right child, so max is stable when moving up.
				leftNode = (leftNode-1)/2; // move up to parent
			}
			
			//Pushing maxes up both sides in lockstep
			while (leftNode != rightNode-1) { // go up until they siblings
				if (leftNode%2 == 1) { //odd index
					//Left child of it's parent; possibly set parent's max = other sibling's max
					if (maxRangeTree[leftNode+1] /*right sibling*/ > leftMax)
						leftMax = maxRangeTree[leftNode+1];
				} // else, right child, so max is stable when moving up.
				
				if (rightNode%2 == 0) {
					//It's an even index, so it's right child of it's parent.
					//That means go up a left-aiming edge, and possibly set parent's max = other sibling's max
					if (maxRangeTree[rightNode-1] /*left sibling*/ > rightMax)
						rightMax = maxRangeTree[rightNode-1];
				} // else, left child, so max is stable when moving up.
				
				leftNode = (leftNode-1)/2; 
				rightNode = (rightNode-1)/2;
			}
			max = Math.max(leftMax, rightMax);
		}
		return max;
	}
}

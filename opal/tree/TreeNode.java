package opal.tree;

import java.util.ArrayList;

import opal.align.Alignment;

public abstract class TreeNode {

	//public to avoid time spent on getters and setters
	public TreeNode leftChild = null;
	public TreeNode rightChild = null;
	public TreeNode parent = null;
	public Alignment alignment;
	public ArrayList<Integer> leafOrderFromInput = new ArrayList<Integer>(1);
	public int firstLeaf = -1; // holds the index (in the sorted order of tree leaves) of the first leaf of the subtree under this node 
	public int lastLeaf = -1; // same as above, but for last leaf of subtree 
	protected boolean isRoot = true;
	
	private int nextLeafID = 0;
	
	public TreeNode(int rootsCount, int i) {
		leafOrderFromInput.add(new Integer(i));
	}
	
	public TreeNode(int rootsCount) {
	}
	
	final public void setAlignment (Alignment a) {
		alignment = a;
	}
		
	final public void clearOldValues () {
		isRoot = false;
		leafOrderFromInput = null;
		clearDistances();
	}

	final public void setChildren (TreeNode a, TreeNode b) {
		leftChild = a;
		rightChild = b;
		a.parent = b.parent = this;
		leafOrderFromInput.addAll(a.leafOrderFromInput);
		leafOrderFromInput.addAll(b.leafOrderFromInput);
	}

	
	/*final public void resetLeafCounter () {
		nextLeafID = 0;
	}*/

	//must be called before doing tree-based polishing 
	final public void assignLeafRanges () {
		//nextLeafID = 0;
		doLeafCount(0);
	}
	
	final private int doLeafCount (int startLeafID) {
		if (leftChild == null && rightChild == null) {
			firstLeaf = lastLeaf = startLeafID++;
		} else {
			startLeafID = leftChild.doLeafCount(startLeafID);
			startLeafID = rightChild.doLeafCount(startLeafID);
			firstLeaf = leftChild.firstLeaf;
			lastLeaf = rightChild.lastLeaf;
		}
		return startLeafID;
	}
	
	final public int fillNodeList (int index, TreeNode[] list) {
		if (leftChild != null) {
			index = leftChild.fillNodeList(index, list);
			index = rightChild.fillNodeList(index, list);
		}
		list[index++] = this;
		return index;
	}
	
	final public int fillInternalNodeList (int index, TreeNode[] list) {
		if (leftChild != null) {
			index = leftChild.fillInternalNodeList(index, list);
			index = rightChild.fillInternalNodeList(index, list);		
			list[index++] = this;
		}
		return index;
	}
	
	final public int fillLeafList (int index, TreeNode[] list) {
		if (leftChild != null) {
			index = leftChild.fillLeafList(index, list);
			index = rightChild.fillLeafList(index, list);
		} else {
			list[index++] = this;
		}
		return index;
	}
	
	final public void buildTreeString (StringBuffer sb, String[] names, ArrayList<Integer> leafOrderList) {
		if (null == leftChild) {
			sb.append(names[leafOrderList.get(nextLeafID)]);
			nextLeafID++;
		} else {
			sb.append("(");
			leftChild.buildTreeString(sb, names, leafOrderList);
			sb.append(",");
			rightChild.buildTreeString(sb, names, leafOrderList);
			sb.append(")");
		} 

		
	}	
	
	public abstract float getDistance (int i);
	public abstract void setDistance (int i, float d) ;
	public abstract void removeDistances (int i, int j) ;
	public abstract void clearDistances () ;
		
	//	public abstract void addDSTPointers(HeapNode heapNode) ;
	//	public abstract ArrayList clearDSTPointers();
	
}

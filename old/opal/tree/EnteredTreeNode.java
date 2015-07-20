package opal.tree;

public class EnteredTreeNode extends TreeNode {
	
	public EnteredTreeNode(int rootsCount, int i) {
		super(rootsCount,i);
	}
	
	public EnteredTreeNode(int rootsCount) {
		super(rootsCount);
	}

	final public void setDistance (int i, float d) {
	}

	final public float getDistance (int i) {
		return 0;
	}

	//this removes entries for two nodes that are about to be removed from "roots",
	// and adds a new entry to the end, which will be used when a new merge is performed 
	final public void removeDistances (int i, int j) {
	}	
	
	final public void clearDistances () {
	}
	
}

package opal.tree;

import java.util.ArrayList;

public class MST_TreeNode extends TreeNode {

	ArrayList<Float> distances ;
	
	public MST_TreeNode(int rootsCount, int i) {
		super(rootsCount,i);
		initList(rootsCount);
	}
	
	public MST_TreeNode(int rootsCount) {
		super(rootsCount);
		initList(rootsCount);
	}
	
	private void initList (int rootsCount) {
		int x = rootsCount+1; // one extra, to allow a new node to be added later
		distances = new ArrayList<Float>(x); 

		// for some reason, initializing the ArrayList(i) creates an empty list
		// ... and distances.ensureCapacity() isn't resetting the size (at least on Blackdown), 
		// so do it manually. 
		x -= distances.size();
		for (int j=0; j<x; j++) {
			distances.add(Float.valueOf(0));
		}
	}
	

	final public void setDistance (int i, float d) {
		distances.set(i, Float.valueOf(d));
	}

	final public float getDistance (int i) {
		return ((Float)distances.get(i)).floatValue();
	}

	//this removes entries for two nodes that are about to be removed from "roots",
	// and adds a new entry to the end, which will be used when a new merge is performed 
	final public void removeDistances (int i, int j) {
		distances.remove(i);
		distances.remove(j);		
		distances.add(Float.valueOf(0));
	}	
	
	final public void clearDistances () {
		distances = null;
	}
	
}

package opal.tree;

public class TreeNodePair implements Comparable<TreeNodePair> {

	public TreeNode A, B; // public to avoid time spent on getters and setters
	float distance;
	
	public TreeNodePair (TreeNode A, TreeNode B, float distance) {
		this.A = A;
		this.B = B;
		this.distance = distance;
	}
	
	public float getDistance () {
		return distance;
	}
		
	public int compareTo(TreeNodePair tp) {

		float otherDist = tp.getDistance();
		return distance<otherDist ? -1 : distance>otherDist ? 1 : 0 ;
	}

}

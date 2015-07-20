package opal.polish;


import opal.tree.TreeNode;
import opal.align.Aligner;

public abstract class TreePolisher extends Polisher {
	public TreePolisher(TreeNode treeNode, Aligner al, int verbosity) {
		super(treeNode, al, verbosity);
	}

	public abstract void polish ( int reps ) ;


}

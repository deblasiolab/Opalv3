package opal.polish;

import java.util.Random;

import opal.align.Aligner;
import opal.align.Aligner.AlignmentType;
import opal.tree.TreeNode;

public abstract class Polisher {

	public static enum PolishType {exhaust_twocut, exhaust_threecut, random_tree_twocut, random_twocut, random_threecut};
	public static PolishType polishMethod = PolishType.random_tree_twocut;
	public static int polishIterations = -1;
	public static int polishIterations_exact = -1;
	public static boolean justPolish = false;
	
	public static AlignmentType polishAligmentMethod = null;// AlignmentType.profile;
	
	TreeNode root;
	Aligner aligner;
	public boolean keepGoing = true;
	protected int verbosity;
	
	Random randomGenerator;
	static long inSeed;
	static boolean inSeedAssigned = false;
	long seed;
	
		
	public static void setRandomSeed (long s) {
		inSeedAssigned = true;
		inSeed = s;
	}
	
	public long getRandomSeed () {
		return seed;
	}

//	public Polisher(TreeNode treeNode, Tree tree, Aligner al, int verbosity) {
//		this(treeNode, al); // no need for tree
//	}
	
	public Polisher(TreeNode treeNode, Aligner al, int verbosity) {
		root = treeNode;
		aligner = al;
		this.verbosity = verbosity;
	}

	public abstract void polish ( int reps ) ;

}

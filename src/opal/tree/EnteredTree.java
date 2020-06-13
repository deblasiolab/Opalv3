package opal.tree;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import com.traviswheeler.libs.LogWriter;
import opal.align.Aligner;
import opal.align.Alignment;
import opal.align.ConsistencyAligner;
import opal.exceptions.GenericOpalException;
import opal.makers.AlignmentMaker;
import opal.IO.Configuration;

public class EnteredTree extends Tree {
	
	public static String treeString;
	public static String treeFile;
	public static String names[];
	
	public EnteredTree (Alignment[] alignments, Aligner al, Distance distance, int currIteration, int verbosity, Configuration c) {
		super(alignments, al, distance, currIteration, verbosity, c);
	}
	
	protected void initializeStructures (int K, Alignment[] alignments, int currIteration) {
		parseString (treeString, alignments);
		
		TreeNode root = roots.get(0);
		
		root.leafOrderFromInput.addAll(root.leftChild.leafOrderFromInput);
		root.leafOrderFromInput.addAll(root.rightChild.leafOrderFromInput);
		root.leftChild.clearOldValues();
		root.rightChild.clearOldValues();

		root.fillInternalNodeList(0, nodesToMerge);
		
		root.fillLeafList (0, leaves);
		
		if (defaultAligner instanceof ConsistencyAligner) {
			distances = new float[K][K];
			for (int i=0; i<K; i++) {
				TreeNode node_i = leaves[i];
				for (int j=i; j<K; j++) {
					TreeNode node_j = leaves[j];
					float d = dist.calcDistance( node_i.alignment, node_j.alignment);					
					int a = root.leafOrderFromInput.get(i);
					int b = root.leafOrderFromInput.get(j);
					distances[a][b] = distances[b][a] = d;
				}
			}
		}		
	}
	
	
	public static void setTreeFromFile (String filename) {
		treeFile = filename;
		InputStream is = null;
		try {
			is = new FileInputStream(filename);
		} catch (FileNotFoundException e) {
			LogWriter.stdErrLogln("The tree file '" + filename + "' cannot be found");
			throw new GenericOpalException(e.getMessage());
		}
	
		byte b[] = null;
		try {
			int x= is.available();
			b= new byte[x];
			is.read(b);
			is.close();
		} catch ( IOException e) {
			LogWriter.stdErrLogln("Error reading tree file '" + filename + "'");
			throw new GenericOpalException(e.getMessage());			
		}
		treeString = new String(b);
	}

	public static void setTreeString (String newick) {
		treeString = newick;
	}
	
	private void parseString (String st, Alignment[] alignments) {
	
		TreeNode newNode;
		TreeNode currentNode = new EnteredTreeNode(0);
		roots.add(currentNode);

		Pattern pattern = Pattern.compile("[(),]");
		Matcher matcher=pattern.matcher(st);

		int prev = -1;
		int curr;
		while (matcher.find() ) {
			curr = matcher.start();
			String ch = st.substring(curr, curr+1);
			String name = st.substring(prev+1,curr);
			prev = curr;

			if (ch.equals("(")) {
				newNode = new EnteredTreeNode(0);
				currentNode.leftChild = newNode;
				newNode.parent = currentNode;
				currentNode = newNode;
			} else if (ch.equals(",")) {
				assignPosition(currentNode, name, alignments);
				newNode = new EnteredTreeNode(0);
				currentNode.parent.rightChild = newNode;
				newNode.parent = currentNode.parent;
				currentNode = newNode;
			} else if (ch.equals(")")) {
				assignPosition(currentNode, name, alignments);
				currentNode = currentNode.parent;
			} 		
		}		
	}

	private void assignPosition (TreeNode node, String name, Alignment[] a) {

		if (!name.equals("") && name.indexOf(":")>-1) 
			name = name.substring(0,name.indexOf(":"));//trim off branch length annotation
		
		if (name.equals("")){ // internal node
			node.leafOrderFromInput.addAll(node.leftChild.leafOrderFromInput);
			node.leafOrderFromInput.addAll(node.rightChild.leafOrderFromInput);
			node.leftChild.clearOldValues();
			node.rightChild.clearOldValues();
			return; 
		}
		
		
		//could speed this up with sort/binary-search ... but this is certainly not a bottleneck
		boolean found = false;
		for (int i=0; !found && i<names.length; i++) {
			if (name.equals(names[i])) {
				found = true;
				node.leafOrderFromInput = new ArrayList<Integer>(1);
				node.leafOrderFromInput.add(Integer.valueOf(i));
				if ( AlignmentMaker.initAlignmentProvided) {
					node.setAlignment(null);
				} else {
					node.setAlignment(a[i]);
				}
			}
		}
		if (!found) {
			LogWriter.stdErrLogln("Name from tree file doesn't match a name in the fasta file. Quitting.");
			throw new GenericOpalException("Name from tree file doesn't match a name in the fasta file. Quitting.");
		}	
	}
	
	final protected void initializePriorityDS() {}
	final protected void killPriorityDS() {}
	final protected void addToPriorityDS(TreeNodePair pair) {}
	
	final protected void pickNextPair() {}
	
	final protected void updateDistances(TreeNode newNode) { }
	
}

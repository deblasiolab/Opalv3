package opal.makers;

import opal.align.Aligner;

import opal.IO.Configuration;
import opal.IO.Inputs;


public abstract class AlignmentMaker {
	public static int outputWidth = -1;
	Aligner al = null;

	public static enum OutputOrderType {tree, input};
	public static OutputOrderType outputOrderMethod = OutputOrderType.input;
	
	//public static String outputName = null;
	public static boolean showCost = false;
	public static boolean initAlignmentProvided = false;
	
/*	public AlignmentMaker (int outputWidth) {
		this.outputWidth = outputWidth;
	}
	*/
	public static void setPolishIterations (int i) {

	}

	/*public static void setOutputName(String name) {
		outputName = name;
	}*/


	//single input ... mulitiple alignment
	abstract public void initialize (Configuration c, Inputs i); 
	
	//two inputs ... aligning alignments
	abstract public void initialize (String fileA, String structFileA, 
					String fileB, String structFileB); 

	abstract public int[][] buildAlignment ();

	public boolean printOutput(int[][] alignmentInstance, String fname){
		return printOutput(alignmentInstance, fname, true);
	}
	abstract public boolean printOutput(int[][] alignmentInstance, String fname, boolean printRealignmentLines);

}

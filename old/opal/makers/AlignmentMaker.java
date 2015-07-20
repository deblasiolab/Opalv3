package opal.makers;

import opal.align.Aligner;


public abstract class AlignmentMaker {
	public static int outputWidth = -1;
	Aligner al = null;

	public static enum OutputOrderType {tree, input};
	public static OutputOrderType outputOrderMethod = OutputOrderType.input;
	
	public static String outputName = null;
	public static boolean showCost = false;
	public static boolean initAlignmentProvided = false;
	
/*	public AlignmentMaker (int outputWidth) {
		this.outputWidth = outputWidth;
	}
	*/
	public static void setPolishIterations (int i) {

	}

	public static void setOutputName(String name) {
		outputName = name;
	}


	//single input ... mulitiple alignment
	abstract public void initialize (String file, String structFile); 
	
	//two inputs ... aligning alignments
	abstract public void initialize (String fileA, String structFileA, 
					String fileB, String structFileB); 

	abstract public int[][] buildAlignment ( String costName, int verbosity, boolean toUpper);

}

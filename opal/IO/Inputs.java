/**
 * 
 */
package opal.IO;

/**
 * @author deblasio
 *
 */
public class Inputs {

	/**
	 * 
	 */
	public int verbosity;
	public boolean toUpper; 
	public boolean justDoConvert;
	public boolean justDoSubOpt;
	public boolean justTree;
	public String fileA;
	public String fileB;
	public String structFileA;
	public String structFileB;
	public String configOutputFile;
	public String bestOutputFile;
	public String bestOutputFileIncludePreRealignment;
	public String featureOutputFile;	
	public String preRealignmentOutputFile;
	public String bestPreRealignmentOutputFile;
	public String bestPreRealignmentsRealignmentOutputFile;
	public String bestPreRealignmentsRealignmentOutputFileIncludePreRealignment;
	
	public StructureFileReader structure;
    
	public Inputs() {
		// TODO Auto-generated constructor stub
	}
	
	public Inputs(Inputs in){
		verbosity = in.verbosity;
		toUpper = in.toUpper; 
		justDoConvert = in.justDoConvert;
		justDoSubOpt = in.justDoSubOpt;
		justTree = in.justTree;
		fileA = in.fileA;
		fileB = in.fileB;
		structFileA = in.structFileA;
		structFileB = in.structFileB;
		configOutputFile = in.configOutputFile;
		bestOutputFile = in.bestOutputFile;
		featureOutputFile = in.featureOutputFile;
		preRealignmentOutputFile = in.preRealignmentOutputFile;
		bestPreRealignmentOutputFile = in.bestPreRealignmentOutputFile;
		bestPreRealignmentsRealignmentOutputFile = in.bestPreRealignmentsRealignmentOutputFile;

		bestOutputFileIncludePreRealignment = in.bestOutputFileIncludePreRealignment;
		bestPreRealignmentsRealignmentOutputFileIncludePreRealignment = in.bestPreRealignmentsRealignmentOutputFileIncludePreRealignment;
		
		structure = new StructureFileReader(in.structure);
		
	}

}

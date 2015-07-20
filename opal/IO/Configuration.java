/**
 * 
 */
package opal.IO;

/**
 * @author deblasio
 * generated on 15 March 2015
 * 
 * The intended use of a configuration object is to store the parameters 
 * for a single alignment.  An array of these would be used when aligning 
 * the same input under different parameters for advising.
 */

import java.util.Scanner;
import java.util.StringTokenizer;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import com.traviswheeler.libs.LogWriter;

import opal.align.Aligner;
import opal.align.Alignment;
import opal.align.StructureAlignment;
import opal.exceptions.GenericOpalException;
public class Configuration {
	public int gamma = -1;
	public int lambda = -1;
	public int gammaTerm = -1;
	public int lambdaTerm = -1;
	public CostMatrix cost;
	public SequenceConverter sc;
	public int repetition = -1;
	public boolean useStructure = false;
	
	public opal.align.StructureAlignment.ParamModel modelType;
	public int gapLevelCnt = 1;
 	public int[] gapOpenMods;
	public int[] gapExtMods;

	public String toString(){
		if(useStructure){
			//if(modelType != null) return modelType.toString();
			/*String rtn = gapLevelCnt + "." + gamma + "." + gammaTerm + ".";
			for(int i=0;i<gapLevelCnt;i++){
				if(i!=0) rtn += "_";
				rtn += gapOpenMods[i];
			}
			
			rtn += "." + lambda + "." + lambdaTerm + ".";

			for(int i=0;i<gapLevelCnt;i++){
				if(i!=0) rtn += "_";
				rtn += gapExtMods[i];
			}
			
			rtn += "." + subHelixHelix + "." + subHelixSheet + "." + subSheetSheet + "." + subHelixLoop + "." + subSheetLoop + "." + subLoopLoop;
			return rtn;*/
			return modelType.toString() + "." + cost.costName + "." + gamma + "." + gammaTerm + "." + lambda + "." + lambdaTerm;
		}
		return cost.costName + "." + gamma + "." + gammaTerm + "." + lambda + "." + lambdaTerm;
	}
	
	public Configuration(){
		initialize(null, -1, -1, -1, -1);
	}
	
	public Configuration(String input){
		initializeFromString(input);
	}
	
	public Configuration(Configuration c){
		initialize(c.cost.costName, c.gamma, c.gammaTerm, c.lambda, c.lambdaTerm);
	}

	public Configuration(String matrix, int ga, int ga_t, int la, int la_t, String fileA){
		if(matrix == null || matrix.equals("")){
			if(useStructure) matrix="VTML200";
			else matrix = getCostName(fileA);
		}
		initialize(matrix, ga, ga_t, la, la_t);
	}
	
	public Configuration(String matrix, int ga, int ga_t, int la, int la_t){
		initialize(matrix, ga, ga_t, la, la_t);
	}
	
	private void initialize(String matrix, int ga, int ga_t, int la, int la_t){
		cost = new CostMatrix();
		if(matrix == null) matrix = "VTML200";
		if(matrix.equals("BLOSUM62_new")) matrix = "BLOSUM62";
		
		cost.initialize(matrix, this);
		
		if (useStructure) {
			cost.isDNA = false;
		} else {
			if (matrix.equals("DNA")) {
				cost.isDNA = true;
				if (ga == -1 ) gamma = CostMatrix.dnaDefaultGamma;
				else gamma = ga;
				gammaTerm = gamma;
				
				if (la == -1 )  lambda = CostMatrix.dnaDefaultLambda;   
				else lambda = la;
				lambdaTerm = lambda;
			} else { //blosum62
				cost.isDNA = false;
				
				if (ga == -1 ) gamma = CostMatrix.protDefaultGamma;
				else gamma = ga;
				
				if (la == -1 )  lambda = CostMatrix.protDefaultLambda;   
				else lambda = la;
				
				if (ga_t == -1 ) gammaTerm = CostMatrix.protDefaultGammaTerm;
				else gammaTerm = ga_t;
				
				if (la_t == -1 )  lambdaTerm = CostMatrix.protDefaultLambdaTerm;   
				else lambdaTerm = la_t;
			}		
		}
		
		sc = new SequenceConverter(cost.getChars());
	}
	
	public void initializeFromString(String input){
		//System.err.println("Configuration: " + input);
		int numTokens =(new StringTokenizer(input,"\\.")).countTokens();
		if(numTokens==1){
			try{
				useStructure = true;
				StructureAlignment.setParams(StructureAlignment.ParamModel.valueOf(input), this);
				
				cost = new CostMatrix();
				cost.initialize("VTML200", this);
				sc = new SequenceConverter(cost.getChars());
				
				//System.out.println("New Configuration for Structure: " + input );
				//System.out.println("Set the subHelixHelix: " + subHelixHelix);
			}catch(java.lang.IllegalArgumentException e){
				LogWriter.stdErrLogln("Configuration is formatted as structure parameter model, but model name "+input+" is not found.");
				throw new GenericOpalException("Configuration is formatted as structure parameter model, but model name "+input+" is not found.");
			}
		}else if(numTokens==5){
			cost = new CostMatrix();
			Scanner s = new Scanner(input).useDelimiter("\\.");
			String temp_cost_name = s.next();
			if(temp_cost_name.equals("BLOSUM62_new")) temp_cost_name = "BLOSUM62";
			cost.initialize(temp_cost_name, this);
			sc = new SequenceConverter(cost.getChars());
			gamma = s.nextInt();
			gammaTerm = s.nextInt();
			lambda = s.nextInt();
			lambdaTerm = s.nextInt();
			s.close();
		}else if(numTokens==6){
			cost = new CostMatrix();
			Scanner s = new Scanner(input).useDelimiter("\\.");
			
			String model = s.next();
			useStructure = true;
			StructureAlignment.setParams(StructureAlignment.ParamModel.valueOf(model), this);
			
			String temp_cost_name = s.next();
			if(temp_cost_name.equals("BLOSUM62_new")) temp_cost_name = "BLOSUM62";
			cost.initialize(temp_cost_name,this);
			sc = new SequenceConverter(cost.getChars());
			gamma = s.nextInt();
			gammaTerm = s.nextInt();
			lambda = s.nextInt();
			lambdaTerm = s.nextInt();
			s.close();
		}else{
			LogWriter.stdErrLogln("Configuration must be a '.' delimited 5-tuple or structure parameter model, confirguration malformed ("+input+").");
			throw new GenericOpalException("Configuration must be a '.' delimited 5-tuple or structure parameter model, confirguration malformed ("+input+").");
		}
	}
	
	private static String getCostName (String file) {
		SequenceFileReader seqReader = new SequenceFileReader(file, true);
		char[][] seqs = seqReader.getSeqs();

		//this is a list of characters that are amino acids, but not in the DNA
		// ambiguity code list
		Pattern forceAA_pattern = Pattern.compile("[QEILFPZ]",Pattern.CASE_INSENSITIVE);	
        Matcher matcher;       
		for (int i=0; i<seqs.length; i++) {
			String str = String.valueOf(seqs[i]);
			matcher = forceAA_pattern.matcher(str);
			if (matcher.find() ) {
				return "VTML200"; 			   
			}
			//if I wanted to be really anal, I'd account for the exceedingly 
			//rare case where the sequences are proteins, but contain none of those
			//protein-only characters, and maybe count the number of ACGTs, but
			//that could still be wrong.
		}
		return "DNA";
	}
	
	
	public void increaseGapCosts(int increase) {
		lambda += increase;
		lambdaTerm += increase;
		/*
		if (lambdas != null){
			for(int i=0; i<lambdas.length; i++) {
				lambdas[i] += increase;
				lambdaTerms[i] += increase;
			}
		}
		*/
	}

	public void multiplyGapCosts(float x) {
		lambda = Math.round(lambda * x);
		lambdaTerm = Math.round(lambdaTerm * x);
		gamma = Math.round(gamma * x);
		gammaTerm = Math.round(gammaTerm * x);
		
	}	
	
	public int getStructureLevelFromProbability(float probability){
		int lvl = (int)Math.floor( probability * gapLevelCnt);
		if (lvl == gapLevelCnt) lvl--;
		return lvl;
	}
}
package opal.IO;

import opal.align.StructureAlignment.ParamModel;
import opal.exceptions.GenericOpalException;

import com.traviswheeler.libs.LogWriter;

import java.util.*;
import java.io.File;
import java.io.FileNotFoundException;

public class CostMatrix {
	
	public char[] chars;
	public int[][] costs;
	public String costName;
	
	public boolean isDNA = false;

	public static int dnaSubAG = 56;
	public static int dnaSubCT = 53;
	
	public static final int dnaDefaultGamma = 260;
	public static final int dnaDefaultGammaTerm = 100;
	public static final int dnaDefaultLambda = 69;
	public static final int dnaDefaultLambdaTerm = 66;
	
	//VTML200.45.11.42.40
	public static final int protDefaultGamma = 45;
	public static final int protDefaultGammaTerm = 11;
	public static final int protDefaultLambda = 42;
	public static final int protDefaultLambdaTerm = 40;
	

	/*public static final int protDefaultGamma = 60;
	public static final int protDefaultGammaTerm = 16;
	public static final int protDefaultLambda = 38;
	public static final int protDefaultLambdaTerm = 36;*/
		
	public int subHelixHelix;
	public int subHelixSheet;
	public int subSheetSheet;
	public int subHelixLoop;
	public int subSheetLoop;
	public int subLoopLoop;
	
	/*public void initialize (String name, int gamma, int gammaTerm, int lambda, int lambdaTerm) {
	       initialize(name);

	       Aligner.gamma = gamma;
	       Aligner.gammaTerm = gammaTerm;
	       Aligner.lambda = lambda;
	       Aligner.lambdaTerm = lambdaTerm;
	}
	
	public static void initialize (String name, String gammas, String gammaTerms, String lambda, String lambdaTerms logger) {
		initialize(name);
		
		Aligner.gamma = Integer.parseInt(gammas);
		Aligner.gammaTerm = Integer.parseInt(gammaTerms);
		Aligner.lambda = Integer.parseInt(lambda);
		Aligner.lambdaTerm = Integer.parseInt(lambdaTerms);
	}*/

	public void initialize (String name, Configuration conf) {
		costName = name;
		int [][] tmp_cost = null;
		if (name.equals("BLOSUM62")) {
			tmp_cost = setBLOSUM62();
		} else if (name.equals("BLOSUM62_modified")){
			tmp_cost = setBLOSUM62_modified();
		} else if (name.equals("BLOSUM50")){
			tmp_cost = setBLOSUM50();
		} 
        
        // ADDED BY D.DeBlasio Feb 2012
        else if (name.equals("BLOSUM60")){
			tmp_cost = setBLOSUM60();
		} else if (name.equals("BLOSUM70")){
			tmp_cost = setBLOSUM70();
		} else if (name.equals("BLOSUM40")){
			tmp_cost = setBLOSUM40();
		} else if (name.equals("BLOSUM80")){
			tmp_cost = setBLOSUM80();
		} else if (name.equals("BLOSUM45")){
			tmp_cost = setBLOSUM45();
		} else if (name.equals("PAM120")){
			tmp_cost = setPAM120();
		} else if (name.equals("PAM40")){
			tmp_cost = setPAM40();
		} else if (name.equals("PAM250")){
			tmp_cost = setPAM250();
		} else if (name.equals("PAM80")){
			tmp_cost = setPAM80();
		} else if (name.equals("PAM320")){
			tmp_cost = setPAM320();
		} else if (name.equals("VTML20")){
			tmp_cost = setVTML20();
		} else if (name.equals("VTML40")){
			tmp_cost = setVTML40();
		} else if (name.equals("VTML80")){
			tmp_cost = setVTML80();
		} else if (name.equals("VTML120")){
			tmp_cost = setVTML120();
		} else if (name.equals("VTML200")){
			tmp_cost = setVTML200();
		} else if (name.equals("B62_P160")){
			tmp_cost = setB62_P160();
		} 
        
        else if (name.equals("DNA")){
			tmp_cost = setDNA();
		} else {
            // Custom Matrix Added 2015
            // try to open a custom file
            try{
                File f = new File(name);
                if(f.exists() && !f.isDirectory()) {
                    tmp_cost = setCustomMatrix(f);
                }else{
                    throw new FileNotFoundException("File is Directory or Not Found");
                }
            }catch(FileNotFoundException e){
                LogWriter.stdErrLogln("You asked me to open the cost matrix " + name + " ... but the code isn't written or the file does not exist");
                throw new GenericOpalException("You asked me to open the cost matrix " + name + " ... but the code isn't written or the file does not exist");
            }
		}
		int stats[] = setCosts(tmp_cost);
		double multiplier = ((float)stats[1] - stats[0])/88.0;
		if (conf.useStructure && conf.modelType == ParamModel.G1) {
			subHelixHelix = (int)(	-20	* multiplier);
			subHelixSheet = 0;
			subSheetSheet = (int)(	-45	* multiplier);
			subHelixLoop = (int)(	9	* multiplier);
			subSheetLoop = (int)(	5	* multiplier);
			subLoopLoop = 0;			
		} else if (conf.useStructure && conf.modelType == ParamModel.G4) {
			subHelixHelix = (int)(	-15	* multiplier);
			subHelixSheet = 0;
			subSheetSheet = (int)(	-41	* multiplier);
			subHelixLoop = (int)(	14	* multiplier);
			subSheetLoop = (int)(	10	* multiplier);
			subLoopLoop = 0;
		} else if (conf.useStructure && conf.modelType == ParamModel.G6) {
			subHelixHelix = (int)(	-11	* multiplier);
			subHelixSheet = 0;
			subSheetSheet = (int)(	-48	* multiplier);
			subHelixLoop = (int)(	10	* multiplier);
			subSheetLoop = (int)(	9	* multiplier);
			subLoopLoop = 0;			
		} else if (conf.useStructure && conf.modelType == ParamModel.G8) {
			subHelixHelix = (int)(	-12	* multiplier);
			subHelixSheet = 0;
			subSheetSheet = (int)(	-46	* multiplier);
			subHelixLoop = (int)(	15	* multiplier);
			subSheetLoop = (int)(	14	* multiplier);
			subLoopLoop = 0;			
		}
	}
	
    public int[][] setCustomMatrix(File f) throws FileNotFoundException{
        Scanner sc = new Scanner(f);
        sc.useLocale(Locale.US);
        chars = sc.next().toCharArray();
        int[][] tmp = new int[chars.length][];
        for(int i=0; i<chars.length; i++){
            tmp[i] = new int[i+1];
            for(int j=0; j<i+1; j++){
                tmp[i][j] = sc.nextInt();
            }
        }
        sc.close();
        return tmp;
    }
    
    
	public int[][] getCosts () {
		return costs;
	}
	
	public char[] getChars () {
		return chars;
	}

	private int[][] setBLOSUM62_modified() {
		//60,38,15,36  = gamma, lambda, gamma_term, lambda_term
		chars = "ARNDCQEGHILKMFPSTWYVBZX".toCharArray();
		return getBLOSUM62_modified_matrix();
	}
	
	public int[][] getBLOSUM62_modified_matrix() {
		int[][] tmp  = {
            {41},
            {73,32},
            {73,66,31},
            {75,73,57,29},
            {67,81,79,84,11},
            {69,59,64,65,80,33},
            {69,66,66,55,85,55,35},
            {65,77,67,74,80,76,77,31},
            {74,65,61,70,78,61,66,76,19},
            {73,80,82,85,72,79,83,87,83,39},
            {74,77,82,84,72,76,80,87,79,55,41},
            {70,53,65,67,82,57,60,74,68,80,78,37},
            {69,73,74,81,71,67,75,83,71,58,53,71,30},
            {77,81,79,85,75,82,85,84,73,65,61,84,63,27},
            {69,76,74,73,80,72,71,76,77,80,80,71,77,83,21},
            {59,70,62,65,70,66,66,67,69,78,79,66,73,78,69,40},
            {65,71,65,71,72,67,69,75,74,69,73,67,69,77,70,56,36},
            {80,78,81,88,80,77,83,80,73,76,72,80,75,56,84,81,80,0},
            {75,73,75,80,75,72,76,82,56,72,70,75,70,47,81,74,73,54,23},
            {66,79,80,83,70,76,79,83,82,49,60,77,63,69,78,74,66,79,72,41},
            {74,70,45,42,82,65,60,71,66,83,83,66,78,82,73,64,68,85,78,82,43},
            {69,63,65,59,83,46,43,77,64,82,78,59,72,83,71,66,68,81,74,78,64,44},
            {68,70,69,72,74,68,70,74,71,71,70,69,68,73,73,68,68,76,71,69,71,69,70}
				//*A*/		{40},
				//*R*/		{72,32},
				//*N*/		{72,68,32},
				//*D*/		{76,72,56,28},
				//*C*/		{68,84,80,84,12},
				//*Q*/		{68,60,64,64,80,32},
				//*E*/		{68,64,64,56,84,52,36},
				//*G*/		{64,76,68,72,80,76,76,32},
				//*H*/		{72,64,60,72,80,60,64,76,20},
				//*I*/		{72,80,84,84,72,80,84,88,84,40},
				//*L*/		{72,76,84,84,72,76,80,84,80,56,40},
				//*K*/		{68,52,64,68,84,56,60,72,68,80,80,36},
				//*M*/		{68,72,76,84,72,68,76,80,72,56,52,72,32},
				//*F*/		{76,80,80,84,80,84,84,84,72,64,60,84,64,28},
				//*P*/		{68,76,76,72,80,72,72,76,76,80,80,72,80,84,20},
				//*S*/		{56,68,60,64,68,64,64,64,68,80,80,64,72,80,68,40},
				//*T*/		{64,72,64,72,68,68,68,72,76,68,72,68,68,76,72,56,36},
				//*W*/		{80,80,88,88,76,76,80,80,80,80,72,80,72,60,84,80,80,0},
				//*Y*/		{76,76,76,84,80,72,76,84,52,72,72,76,68,48,80,76,72,52,24},
				//*V*/		{64,80,80,84,68,76,80,84,84,48,60,76,60,68,80,72,64,80,72,40},
				//*B*/		{72,72,44,40,84,64,60,68,68,84,84,68,80,84,76,64,68,88,80,84,44},
				//*Z*/		{68,64,64,60,84,44,40,76,64,84,80,60,72,84,72,64,68,80,76,80,64,44},
				//*X*/		{68,72,72,72,76,68,68,72,72,72,72,68,68,72,72,68,68,76,72,68,72,68,72}
				/*            A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X*/
				};
			  
				return tmp;
	}
	

	private int[][] setBLOSUM62() {
		//68,43,17,41  = gamma, lambda, gamma_term, lambda_term
		chars = "ARNDCQEGHILKMFPSTWYVBZX".toCharArray();
		int[][] tmp  = {
			{46},
			{83,36},
			{83,75,36},
			{85,83,64,33},
			{76,92,90,95,12},
			{78,67,73,74,91,38},
			{79,75,75,62,96,62,40},
			{73,88,77,84,90,86,87,35},
			{85,74,69,79,88,69,75,86,22},
			{82,91,93,96,81,90,95,98,94,44},
			{84,88,93,95,82,86,91,98,89,63,46},
			{80,60,73,77,93,65,68,84,78,91,89,42},
			{78,83,84,92,80,76,85,94,81,66,60,81,35},
			{88,92,90,97,86,93,96,95,82,74,69,95,72,31},
			{79,86,84,83,91,82,81,87,88,91,91,81,87,94,24},
			{67,79,70,74,79,75,75,77,78,89,89,75,83,89,78,46},
			{74,80,74,81,81,76,79,85,84,79,83,77,79,88,79,63,41},
			{91,89,92,100,91,88,94,91,83,87,81,91,85,64,95,93,91,0},
			{85,83,85,91,86,82,86,93,64,82,80,85,80,54,92,85,83,61,27},
			{75,89,90,95,79,87,90,95,93,56,68,87,71,78,88,85,75,90,82,46},
			{84,80,51,48,93,74,68,80,75,95,94,75,88,94,83,72,77,96,88,93,49},
			{78,72,74,67,94,53,49,87,73,93,89,67,81,95,81,75,78,92,84,89,72,50},
			{77,80,79,82,84,78,79,84,81,80,80,78,78,83,83,77,77,87,81,79,80,79,80}
		};
	  
		return tmp;
	}

    
    //ADDED BY D. DeBlasio Feb 2012
    
    private int[][] setBLOSUM60() {
		chars = "ARNDCQEGHILKMFPSTWYVBZX".toCharArray();
		int[][] tmp  = {
            {47},
            {84,36},
            {83,76,36},
            {85,84,64,34},
            {76,92,88,96,11},
            {78,68,73,74,91,38},
            {80,75,76,62,97,63,40},
            {74,89,77,84,91,87,88,35},
            {84,73,68,78,88,68,76,86,22},
            {83,92,93,97,83,90,95,99,94,45},
            {84,88,93,96,83,86,92,100,89,63,47},
            {81,60,73,76,94,66,68,85,77,92,90,43},
            {79,84,84,92,81,77,85,94,82,66,60,82,35},
            {89,93,90,99,85,93,97,96,82,74,69,95,72,31},
            {80,88,84,83,91,82,81,87,88,92,92,82,88,95,24},
            {67,80,70,75,79,76,76,77,79,89,90,76,84,90,79,46},
            {74,81,74,81,81,77,80,86,83,79,83,77,80,88,80,64,42},
            {92,89,91,99,92,87,94,91,82,87,82,92,84,64,95,93,92,0},
            {86,83,85,90,85,83,87,94,64,82,79,85,79,54,92,86,84,61,27},
            {75,90,91,96,79,88,91,95,93,56,69,88,71,79,89,85,75,89,82,47},
            {84,80,51,47,92,74,68,81,74,95,94,75,88,95,83,73,78,95,88,94,49},
            {79,72,75,67,94,54,49,88,73,93,90,67,82,95,82,76,79,91,85,90,73,51},
            {77,80,79,82,85,78,80,84,81,81,81,79,78,83,84,77,77,87,81,79,81,79,80}
		};
		return tmp;
	}
    
    
    private int[][] setBLOSUM70() {
		chars = "ARNDCQEGHILKMFPSTWYVBZX".toCharArray();
		int[][] tmp  = {
            {44},
            {83,34},
            {83,76,33},
            {85,83,63,32},
            {77,93,92,95,13},
            {78,66,72,74,91,35},
            {78,75,75,62,99,61,37},
            {72,88,77,83,92,86,87,34},
            {84,73,69,80,92,68,75,88,21},
            {82,91,94,97,82,91,94,99,94,43},
            {84,89,94,97,82,87,92,98,90,62,45},
            {79,59,73,77,94,64,67,84,77,91,88,40},
            {78,83,85,93,80,73,85,93,82,65,59,81,33},
            {88,92,92,96,86,93,96,94,82,74,68,95,72,30},
            {78,86,85,83,92,82,81,88,88,92,91,80,89,95,23},
            {65,78,69,74,80,74,74,76,79,88,89,75,82,88,77,44},
            {72,80,73,80,81,76,78,85,83,78,83,76,78,87,80,62,39},
            {91,88,93,100,94,87,94,90,83,87,83,91,85,65,97,91,91,0},
            {86,83,86,91,86,83,86,94,62,82,79,85,78,53,92,84,83,60,25},
            {74,90,90,94,79,87,89,95,93,55,67,88,69,77,88,84,74,90,82,45},
            {84,80,49,46,94,73,68,80,75,96,96,75,89,94,84,72,77,97,89,93,48},
            {78,71,74,67,96,51,47,87,72,93,90,66,80,95,81,74,77,92,85,88,72,48},
            {76,80,79,82,85,77,79,84,81,80,80,78,77,82,83,76,76,87,81,78,81,78,80}
		};
		return tmp;
	}
    
    private int[][] setBLOSUM40() {
		chars = "ARNDCQEGHILKMFPSTWYVBZX".toCharArray();
		int[][] tmp  = {
            {52},
            {85,38},
            {81,79,44},
            {86,82,70,36},
            {81,94,86,91,4},
            {81,70,77,79,88,41},
            {83,79,79,66,90,69,46},
            {77,89,79,86,92,89,88,38},
            {87,78,71,78,89,73,81,90,19},
            {83,91,91,96,88,91,92,96,93,48},
            {86,88,90,94,87,84,90,96,89,65,50},
            {83,64,74,79,94,71,70,85,79,91,90,49},
            {79,86,83,93,83,82,84,93,78,69,64,85,39},
            {91,92,87,93,84,95,93,95,83,74,70,93,73,34},
            {84,91,85,86,93,83,81,84,89,90,94,85,89,91,27},
            {73,81,75,80,82,79,79,79,85,89,87,79,86,85,84,52},
            {79,84,75,84,85,82,83,85,85,81,82,78,83,85,83,67,49},
            {91,82,85,95,100,82,89,88,91,82,77,84,89,67,85,93,87,0},
            {88,83,85,87,91,83,88,95,74,79,80,87,78,58,89,86,84,60,30},
            {77,86,90,94,81,89,91,93,93,58,72,90,75,77,90,86,77,86,83,50},
            {84,81,58,51,89,78,72,83,75,94,92,77,88,91,86,78,80,90,86,92,54},
            {82,76,78,70,90,58,54,89,78,92,88,70,83,93,82,79,82,86,86,90,77,56},
            {80,81,80,82,86,80,81,84,83,82,81,80,80,83,84,79,79,84,82,81,81,80,81}
		};
		return tmp;
	}
    
    private int[][] setBLOSUM80() {
		chars = "ARNDCQEGHILKMFPSTWYVBZX".toCharArray();
		int[][] tmp  = {
            {41},
            {80,32},
            {80,73,30},
            {82,81,61,30},
            {75,91,90,95,14},
            {75,63,69,72,89,32},
            {75,72,72,60,98,58,35},
            {70,85,74,80,92,83,84,32},
            {81,70,67,78,93,65,71,85,20},
            {79,89,92,95,79,88,91,97,92,40},
            {80,86,92,96,80,84,90,95,87,59,42},
            {74,56,69,75,92,62,65,81,74,88,86,37}, 
            {76,80,82,91,80,70,82,90,81,62,56,78,30},
            {85,90,90,93,85,90,93,92,79,72,66,92,69,28},
            {73,82,83,81,90,80,77,86,84,89,88,78,86,92,21},
            {62,75,66,72,79,71,71,73,76,85,86,72,78,85,75,40},
            {69,77,71,78,79,73,75,82,80,75,80,73,74,83,77,59,36},
            {90,86,92,100,95,87,93,90,82,86,82,90,82,65,97,90,88,0},
            {83,82,84,90,88,81,85,93,58,80,77,83,77,50,90,82,81,56,24},
            {71,87,88,92,77,84,86,92,89,52,65,85,66,75,85,81,71,88,81,42},
            {81,77,47,44,93,71,65,77,73,94,94,72,87,92,82,69,75,96,87,90,45},
            {75,68,71,65,95,48,44,84,69,90,87,64,77,92,78,71,75,91,83,85,69,45},
            {73,77,77,80,84,75,76,81,78,78,77,75,74,80,81,73,73,86,79,76,78,76,77}
		};
		return tmp;
	}
    
    private int[][] setBLOSUM45() {
		chars = "ARNDCQEGHILKMFPSTWYVBZX".toCharArray();
		int[][] tmp  = {
            {51},
            {87,39},
            {84,79,43},
            {89,84,69,36},
            {80,95,86,96,7},
            {82,71,77,76,90,41},
            {84,78,80,65,94,68,45},
            {77,90,79,87,93,90,90,38},
            {86,78,72,79,89,74,81,89,23},
            {85,95,92,98,89,91,94,100,95,48},
            {87,91,93,96,88,86,92,100,88,66,50},
            {84,63,75,80,94,70,71,87,80,95,92,46},
            {82,87,85,93,85,82,86,95,83,70,64,86,39},
            {91,93,89,97,86,95,95,96,83,76,70,95,74,34},
            {85,90,86,86,92,85,81,87,88,91,95,86,91,93,26},
            {72,83,74,80,82,80,81,79,84,90,91,80,88,89,84,51},
            {79,86,76,84,83,83,82,87,85,82,85,81,83,87,83,67,48},
            {91,84,90,96,97,81,90,89,85,83,81,87,88,66,87,96,91,0},
            {90,83,88,89,90,85,89,96,71,82,82,88,79,57,91,87,85,60,30},
            {77,92,90,98,82,90,92,96,94,58,72,92,74,79,93,87,78,88,83,51},
            {86,81,57,51,91,77,72,83,76,95,95,77,89,93,86,77,80,93,88,94,54},
            {83,75,79,69,93,58,54,90,78,93,90,71,85,95,82,81,83,86,87,92,77,55},
            {80,83,81,84,87,81,82,85,83,83,83,82,81,84,85,80,80,85,83,82,82,81,82}
		};
		return tmp;
	}
    
    private int[][] setPAM120() {
		chars = "ARNDCQEGHILKMFPSTWYVBZX".toCharArray();
		int[][] tmp  = {
            {44},
            {73,28},
            {62,66,40},
            {62,75,48,36},
            {73,79,82,91,17},
            {65,56,61,56,92,33},
            {60,73,56,43,92,48,36},
            {56,79,61,61,81,72,64,36},
            {73,56,52,61,79,47,63,77,26},
            {66,72,71,75,73,75,72,79,78,32},
            {73,79,77,86,95,70,81,85,73,53,33},
            {71,48,56,64,92,60,65,74,67,72,78,36},
            {68,65,74,80,91,66,74,78,78,53,47,58,21},
            {80,84,80,93,88,88,92,85,72,58,58,91,64,23},
            {56,64,68,72,78,61,68,69,65,76,76,71,75,85,30},
            {54,63,55,62,62,67,64,57,68,71,78,64,70,76,57,44},
            {54,70,59,65,75,69,67,66,72,61,73,62,65,79,63,52,39},
            {92,55,81,97,100,88,100,97,76,91,72,83,87,65,92,71,89,0},
            {79,85,69,85,63,85,83,91,62,70,72,84,80,39,90,76,75,67,19},
            {60,76,74,76,72,73,73,70,74,43,56,78,55,72,70,69,61,96,76,36},
            {62,71,44,41,87,58,49,61,57,73,82,60,77,87,70,59,62,90,78,75,59},
            {63,66,58,49,92,41,41,67,56,73,76,63,71,91,65,66,68,95,84,73,46,59},
            {63,68,63,67,79,66,67,68,67,66,71,68,67,76,67,63,63,84,75,66,67,68,68}
		};
		return tmp;
	}
    
    private int[][] setPAM40() {
		chars = "ARNDCQEGHILKMFPSTWYVBZX".toCharArray();
		int[][] tmp  = {
            {26},
            {68,17},
            {57,63,20},
            {56,78,39,20},
            {67,72,80,91,12},
            {59,50,56,52,91,18},
            {53,74,52,36,91,40,21},
            {50,76,55,56,76,68,59,24},
            {69,50,44,57,70,41,61,74,15},
            {60,63,62,70,66,71,64,80,75,17},
            {66,74,69,86,94,62,76,80,66,49,23},
            {68,43,48,59,91,54,59,69,65,65,72,23},
            {61,59,73,80,90,58,68,73,78,47,42,50,8},
            {73,77,74,93,87,88,91,76,65,52,53,90,58,15},
            {50,58,64,70,72,54,63,65,58,73,69,67,72,78,18},
            {46,55,45,57,54,62,59,50,64,67,73,58,63,67,50,25},
            {47,66,51,60,71,63,64,64,69,53,67,56,58,74,58,43,22},
            {89,51,72,95,97,87,100,95,70,90,65,83,86,60,90,62,87,0},
            {71,79,60,82,58,84,74,90,56,66,68,77,81,38,89,68,67,62,12},
            {52,71,71,72,65,67,67,64,66,38,52,75,49,71,64,65,54,96,70,22},
            {56,71,30,29,86,54,43,55,51,66,78,54,77,84,67,51,56,85,72,71,46},
            {55,64,54,43,91,30,29,63,52,67,70,57,64,90,59,61,64,94,78,67,46,46},
            {56,63,57,63,75,60,61,63,62,61,65,62,62,71,62,55,57,81,70,60,62,63,62}
		};
		return tmp;
	}
    
    private int[][] setPAM250() {
		chars = "ARNDCQEGHILKMFPSTWYVBZX".toCharArray();
		int[][] tmp  = {
            {62},
            {75,45},
            {68,69,61},
            {68,74,61,54},
            {77,84,84,90,21},
            {71,64,66,62,91,53},
            {68,73,63,55,90,59,54},
            {64,79,68,67,83,74,68,50},
            {74,63,63,66,83,57,66,77,43},
            {71,77,76,79,78,77,77,79,79,51},
            {77,81,81,85,93,76,82,85,77,59,45},
            {74,55,65,69,91,66,69,76,69,77,80,50},
            {74,71,76,79,90,73,78,80,78,60,54,67,43},
            {83,87,83,92,86,88,91,88,76,65,62,90,68,33},
            {65,70,71,73,80,68,71,71,70,77,79,74,77,87,46},
            {65,70,66,68,69,71,69,65,72,75,80,70,75,82,65,63},
            {64,72,67,70,78,72,71,69,74,69,76,69,71,81,68,64,59},
            {92,60,85,96,100,88,97,97,80,90,77,83,86,68,91,79,90,0},
            {83,86,77,86,68,85,86,90,69,73,73,87,79,41,89,80,80,70,28},
            {68,79,76,78,77,77,76,74,78,54,62,79,62,74,74,73,68,94,79,52},
            {68,72,61,57,87,64,59,67,65,77,83,67,78,88,72,67,69,91,82,77,69},
            {69,69,64,58,90,56,56,71,62,77,80,68,76,89,70,70,71,93,86,76,54,69},
            {70,72,70,72,81,71,72,72,71,72,74,72,72,78,72,70,70,85,78,71,72,73,72}
		};
		return tmp;
	}
    
    private int[][] setPAM80() {
		chars = "ARNDCQEGHILKMFPSTWYVBZX".toCharArray();
		int[][] tmp  = {
            {36},
            {71,23},
            {60,64,31},
            {59,76,43,28},
            {70,76,81,91,15},
            {62,53,58,54,91,25},
            {57,73,53,39,92,44,29},
            {53,78,58,59,79,70,62,30},
            {71,53,47,59,75,44,62,76,20},
            {63,68,67,73,70,73,69,79,76,25},
            {70,77,74,86,95,66,79,83,70,51,29},
            {70,45,52,62,92,57,62,72,66,69,76,30},
            {65,62,73,80,90,62,72,76,78,50,44,54,14},
            {77,81,78,93,87,88,92,81,69,55,55,91,61,19},
            {53,61,66,71,75,58,66,67,62,75,73,69,74,82,24},
            {50,59,50,59,58,65,62,54,66,69,76,61,67,72,54,35},
            {50,68,55,63,73,67,66,65,71,57,70,59,62,77,61,47,31},
            {90,53,77,96,99,87,100,96,73,90,69,83,86,63,91,67,88,0},
            {76,83,65,83,61,84,79,90,59,68,70,81,80,38,89,72,71,65,16},
            {56,74,72,74,69,71,70,67,71,40,54,77,52,72,67,67,57,96,74,29},
            {59,71,37,35,86,56,46,58,54,70,80,57,77,86,69,55,59,88,75,73,53},
            {59,64,56,45,91,36,35,65,54,71,73,60,68,90,62,63,66,95,81,70,49,53},
            {60,66,60,65,77,63,64,66,65,64,68,65,65,74,65,59,60,83,73,64,65,65,65}
		};
		return tmp;
	}
    
    private int[][] setPAM320() {
		chars = "ARNDCQEGHILKMFPSTWYVBZX".toCharArray();
		int[][] tmp  = {
            {68},
            {77,52},
            {71,71,67},
            {71,75,66,60},
            {79,85,85,89,23},
            {73,68,69,66,90,61},
            {71,74,67,61,90,64,61},
            {68,80,71,69,83,75,71,56},
            {76,67,67,69,84,63,69,78,51},
            {74,79,78,80,80,79,79,80,80,58},
            {78,82,82,85,93,79,83,85,79,63,51},
            {75,60,69,71,91,69,72,77,71,79,82,57},
            {76,74,78,80,90,76,79,81,79,64,59,72,53},
            {84,88,84,91,86,88,90,89,78,68,64,90,71,37},
            {68,72,73,74,81,71,73,73,73,78,81,75,79,88,53},
            {69,73,70,71,72,73,72,68,74,77,81,72,77,84,69,68},
            {69,74,71,72,80,74,73,71,76,72,77,72,74,83,71,69,66},
            {93,63,87,96,100,88,96,97,81,90,78,83,86,69,92,82,90,0},
            {84,87,80,87,70,86,87,90,73,74,73,88,79,43,89,82,82,71,33},
            {72,80,78,79,79,78,78,76,79,60,65,80,66,75,76,75,71,94,80,59},
            {71,73,66,63,87,68,64,70,68,79,84,70,79,88,74,71,72,91,84,78,72},
            {72,71,68,63,90,63,62,73,66,79,81,70,78,89,72,72,73,93,87,78,62,72},
            {73,75,73,74,82,74,74,74,74,74,76,74,74,80,74,73,73,86,79,74,75,75,75}
		};
		return tmp;
	}
    
    private int[][] setVTML20() {
		chars = "ARNDCQEGHILKMFPSTWYVBZX".toCharArray();
		int[][] tmp  = {
            {19},
            {66,15},
            {64,60,13},
            {63,85,47,15},
            {53,68,69,92,2},
            {60,50,55,58,88,12},
            {60,78,60,45,91,46,17},
            {56,68,59,62,68,69,65,17},
            {66,54,52,58,65,48,62,68,7},
            {67,72,73,86,61,76,77,97,74,17},
            {66,70,75,95,86,64,72,82,67,48,20},
            {62,42,54,61,88,49,51,68,60,73,70,17},
            {61,64,67,75,58,57,72,77,84,48,44,61,7},
            {70,78,78,100,88,70,93,80,59,59,53,91,52,12},
            {58,66,70,64,74,60,64,70,65,79,68,62,78,73,13},
            {49,62,50,59,54,57,59,58,59,75,71,60,71,66,57,17},
            {53,64,55,62,61,59,62,70,62,59,67,58,57,70,63,46,16},
            {76,70,79,79,97,95,98,74,64,64,64,73,87,52,74,70,94,0},
            {71,67,63,90,58,85,67,79,47,69,64,70,82,41,94,64,70,49,10},
            {53,73,75,75,55,66,69,77,70,40,54,69,54,63,69,70,56,89,69,18},
            {63,73,30,31,81,57,53,61,55,80,85,57,71,89,67,54,59,79,76,75,30},
            {60,64,58,52,89,29,31,67,55,77,68,50,65,82,62,58,61,96,76,68,55,30},
            {60,64,61,69,69,61,65,68,60,66,65,62,63,67,66,59,60,72,65,63,65,63,64},
            {100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,40}
		};
		return tmp;
	}
    
    private int[][] setVTML40() {
		chars = "ARNDCQEGHILKMFPSTWYVBZX".toCharArray();
		int[][] tmp  = {
            {24},
            {68,18},
            {66,62,17},
            {65,83,47,18},
            {55,71,72,91,3},
            {62,51,57,59,87,16},
            {62,76,61,46,90,47,21},
            {58,71,61,64,71,72,68,20},
            {69,55,53,60,68,49,64,70,9},
            {68,76,77,89,63,78,80,96,77,21},
            {69,73,79,94,85,67,76,85,69,49,24},
            {65,42,56,62,87,50,52,70,61,76,73,21},
            {63,67,70,78,60,60,75,81,83,49,44,64,11},
            {73,81,81,100,87,74,92,84,61,61,54,91,53,15},
            {60,68,72,66,76,62,66,73,67,82,71,65,80,77,15},
            {50,65,51,61,55,59,61,60,61,77,74,62,73,69,59,22},
            {55,67,57,64,62,61,64,72,64,61,70,60,59,73,66,47,20},
            {79,73,82,84,97,95,98,78,66,67,66,76,86,53,78,73,94,0},
            {74,70,65,90,60,83,70,83,48,71,66,73,80,41,94,66,72,49,12},
            {55,76,77,78,56,69,72,80,73,41,55,72,55,65,72,72,58,88,72,23},
            {65,73,32,33,81,58,53,63,57,83,86,59,74,91,69,56,60,83,77,78,32},
            {62,64,59,53,88,31,34,70,57,79,71,51,67,83,64,60,63,96,77,71,56,32},
            {62,66,63,70,70,63,67,71,61,68,67,64,64,69,68,61,62,74,67,66,67,65,66},
            {100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,46}
		};
		return tmp;
	}
    
    private int[][] setVTML80() {
		chars = "ARNDCQEGHILKMFPSTWYVBZX".toCharArray();
		int[][] tmp  = {
            {31},
            {71,24},
            {68,65,24},
            {68,82,50,24},
            {57,75,75,91,5},
            {66,54,60,61,86,24},
            {65,74,63,48,89,49,28},
            {61,74,64,67,74,75,71,24},
            {72,58,56,63,72,52,66,74,15},
            {70,80,81,91,65,80,84,97,80,29},
            {72,77,82,94,84,71,80,90,74,51,31},
            {68,44,59,65,87,53,55,73,64,80,77,29},
            {67,71,75,82,64,64,78,85,82,51,46,69,18},
            {76,84,84,100,86,78,93,89,64,63,56,91,56,20},
            {63,71,74,69,79,66,69,76,71,84,76,68,82,81,19},
            {53,68,54,64,58,62,64,63,64,79,78,65,75,74,62,31},
            {58,69,60,67,65,64,67,75,67,65,73,64,63,76,69,50,28},
            {84,78,86,90,97,95,99,83,69,71,69,81,85,54,82,78,94,0},
            {77,74,69,90,64,83,75,87,50,74,68,77,78,41,93,70,76,50,17},
            {59,79,80,82,59,74,76,84,77,43,57,76,57,68,76,74,62,88,75,31},
            {68,73,37,37,83,61,55,65,60,86,88,62,78,92,72,59,63,88,79,81,37},
            {66,64,62,54,88,37,39,73,59,82,76,54,71,85,67,63,66,97,79,75,58,38},
            {65,69,66,72,72,66,70,74,64,71,70,67,67,72,72,64,66,77,69,69,69,68,69},
            {100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,54}
		};
		return tmp;
	}
    
    private int[][] setVTML120() {
		chars = "ARNDCQEGHILKMFPSTWYVBZX".toCharArray();
		int[][] tmp  = {
            {39},
            {74,29},
            {70,67,31},
            {71,81,52,30},
            {60,78,77,90,8},
            {69,56,63,63,86,32},
            {68,74,65,50,89,53,34},
            {64,77,66,69,76,77,73,28},
            {74,61,59,66,75,56,68,76,20},
            {72,83,84,93,68,81,86,97,82,36},
            {75,80,85,95,83,75,83,93,77,53,37},
            {71,47,62,67,87,56,58,76,66,82,80,36},
            {70,74,78,85,67,69,80,88,82,54,49,72,25},
            {79,86,86,100,86,82,93,92,66,65,58,91,58,24},
            {66,73,76,72,81,68,71,78,73,86,79,71,84,85,22},
            {57,70,58,66,61,65,66,65,67,81,80,67,77,77,65,39},
            {61,72,63,70,67,67,69,76,70,68,75,67,67,79,71,53,35},
            {87,81,89,94,97,95,99,87,71,74,71,85,84,55,86,82,94,0},
            {80,77,73,90,67,83,79,89,53,76,70,80,78,42,94,74,79,51,21},
            {62,82,82,85,62,77,79,86,80,46,59,79,59,70,79,76,65,88,77,38},
            {70,74,42,41,84,63,57,68,62,89,90,64,82,93,74,62,66,91,81,84,41},
            {69,65,64,57,87,42,44,75,62,83,79,57,74,87,70,66,68,97,81,78,60,43},
            {68,71,69,74,73,69,72,77,67,73,73,70,70,74,74,67,68,79,72,71,72,70,72},
            {100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,59}
		};
		return tmp;
	}
    
    private int[][] setVTML200() {
		chars = "ARNDCQEGHILKMFPSTWYVBZX".toCharArray();
		int[][] tmp  = {
            {52},
            {77,40},
            {74,71,44},
            {75,80,58,40},
            {65,82,80,90,13},
            {73,62,68,67,86,47},
            {73,74,68,55,89,60,46},
            {68,80,70,73,79,79,76,35},
            {78,66,65,70,79,63,72,80,32},
            {75,86,88,94,71,84,88,97,85,47},
            {78,84,88,95,83,80,87,95,81,57,46},
            {75,53,67,70,87,62,64,79,70,85,84,49},
            {75,79,83,88,73,76,84,91,83,59,54,78,40},
            {83,89,88,100,85,86,94,95,70,69,62,91,63,33},
            {71,77,78,75,83,73,75,81,77,88,83,74,86,89,29},
            {64,74,65,70,67,70,71,70,72,83,84,72,80,82,70,53},
            {66,75,68,73,71,72,73,79,74,73,78,72,73,82,74,61,50},
            {91,86,92,98,97,96,100,92,75,78,75,89,84,56,91,87,95,0},
            {83,81,78,91,72,84,84,92,58,78,73,84,78,45,94,79,82,52,28},
            {69,84,85,88,66,81,83,89,83,52,62,83,63,73,82,79,71,88,80,49},
            {74,76,51,49,85,67,61,71,68,91,91,69,86,94,77,67,71,95,85,86,50},
            {73,68,68,61,87,53,53,78,67,86,83,63,80,90,74,71,73,98,84,82,64,53},
            {73,75,74,77,76,73,76,80,72,77,76,74,75,77,77,73,73,81,75,76,76,75,75},
            {100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,66}
		};
		return tmp;
	}
    
    private int[][] setB62_P160() {
		chars = "ARNDCQEGHILKMFPSTWYVBZX".toCharArray();
		int[][] tmp  = {
            
		};
		return tmp;
	}
    //END ADDITION D.DeBlasio Feb 2012

	
	private int[][] setBLOSUM50() {
		chars = "ARNDCQEGHILKMFPSTWYVBZX".toCharArray();
		return getBLOSUM50matrix();
	}
	public int[][] getBLOSUM50matrix() {
		int[][] tmp  = {
				{49},
				{85,37},
				{83,77,40},
				{87,84,66,34},
				{76,93,85,94,9},
				{79,69,75,74,90,40},
				{81,76,77,63,93,66,42},
				{75,89,78,85,90,87,89,36},
				{85,75,70,77,85,71,78,86,23},
				{84,93,93,98,86,90,94,98,96,46},
				{86,90,92,95,85,86,92,100,87,64,48},
				{82,61,74,77,94,67,69,86,78,93,91,44},
				{80,86,83,91,81,80,85,94,83,68,62,84,37},
				{90,94,90,98,84,94,95,97,83,76,69,95,72,32},
				{82,89,84,84,90,83,81,87,89,91,94,84,90,94,24},
				{70,82,72,77,79,78,80,77,80,90,90,78,85,90,82,48},
				{76,83,75,82,80,79,81,86,84,81,84,79,81,88,81,65,45},
				{91,86,88,96,91,84,91,89,82,83,81,89,86,63,90,95,90,0},
				{88,83,86,90,86,84,86,95,67,81,80,85,79,56,91,87,83,60,29},
				{75,91,91,97,80,89,92,95,93,57,70,90,72,79,91,85,77,88,83,49},
				{85,81,54,49,90,75,69,82,74,95,94,76,88,94,84,75,79,92,88,94,51},
				{81,73,76,68,92,56,51,88,75,93,90,69,83,95,81,79,80,88,85,91,74,53},
				{79,81,79,83,85,79,81,85,81,82,81,80,79,83,84,79,78,85,82,80,81,80,81}		
		};
	  
		return tmp;		
	}
	private int[][] setDNA() {
		//260, 69  = gamma, lambda  
		
		chars = "ACGTUKMRYSWBVHDNX".toCharArray();
		int minSub = Math.min(dnaSubAG, dnaSubCT);
		int[][] tmp  = {
		/*A*/		{   0  },
		/*C*/		{  100,       0},
		/*G*/		{dnaSubAG,   100,      0},
		/*T*/		{  100,   dnaSubCT,   100,       0},
		/*U*/		{  100,   dnaSubCT,   100,       0,        0},
		/*K  G,T*/  {dnaSubAG,dnaSubCT,    0,        0,        0,       0},
		/*M  A,C*/  {   0,        0,    dnaSubAG, dnaSubCT, dnaSubCT, minSub, 0},
		/*R  A,G*/  {   0,       100,      0,       100,      100,      0,    0,   0},
		/*Y  C,T*/  {  100,       0,      100,       0,        0,       0,    0,  100,  0},
		/*S  C,G*/  {dnaSubAG,    0,       0,     dnaSubCT, dnaSubCT,   0,    0,   0,   0,   0},
		/*W  A,T*/ 	{   0,    dnaSubCT, dnaSubAG,    0,        0,       0,    0,   0,   0, minSub, 0},
		/*B  C,G,T*/{dnaSubAG,    0, 	   0,        0,        0,       0,    0,   0,   0,   0,    0,   0},
		/*V  A,C,G*/{ 	0,        0,       0,     dnaSubCT, dnaSubCT,   0,    0,   0,   0,   0,    0,   0,  0},
		/*H  A,C,T*/{ 	0,        0,    dnaSubAG,    0,        0,       0,    0,   0,   0,   0,    0,   0,  0,  0},
		/*D  A,G,T*/{ 	0,    dnaSubCT,    0,        0,        0,       0,    0,   0,   0,   0,    0,   0,  0,  0,  0},
		/*N  any*/  {   0,        0,       0,        0,        0,       0,    0,   0,   0,   0,    0,   0,  0,  0,  0,  0},
		/*X  any*/  {   0,        0,       0,        0,        0,       0,    0,   0,   0,   0,    0,   0,  0,  0,  0,  0,  0}
		
		/*              A         C        G         T         U        K     M    R    Y    S     W    B   V   H   D   N   X*/
		/* these ambiguity costs are optimistic: they give a cost that is the lowest it could be*/
		};
	  
		return tmp;
	}



	private int[] setCosts(int[][] tmp){
		costs = new int[tmp.length][tmp.length];
		int stats[] = new int[3];
		stats[0] = tmp[0][0];
		stats[1] = tmp[0][0];
		int ct = 0;
		
		for (int i=0; i<tmp.length; i++){
			for (int j=i; j<tmp.length; j++){
				costs[i][j] = costs[j][i] = tmp[j][i];	
				if(tmp[j][i]<stats[0]) stats[0] = tmp[j][i];
				if(tmp[j][i]>stats[1]) stats[1] = tmp[j][i];
				stats[2] += tmp[j][i];
				ct++;
			}
		}
		stats[2] /= ct;
		return stats;
	}

	public void increaseCosts(int x) {
		for (int i=0; i<costs.length; i++) 
			for (int j=i; j<costs.length; j++) 
				costs[i][j] = costs[j][i] = costs[j][i] + x; 
	}

	public void multiplyCosts(float x) {
		for (int i=0; i<costs.length; i++) 
			for (int j=i; j<costs.length; j++) 
				costs[i][j] = costs[j][i] = Math.round(costs[j][i] * x); 
	}
	
}
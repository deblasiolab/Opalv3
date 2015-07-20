package opal.IO;

import opal.align.Aligner;
import opal.exceptions.GenericOpalException;
import com.traviswheeler.libs.LogWriter;

public class CostMatrix {
	
	public static char[] chars;
	static int[][] costs;	
	
	public static boolean isDNA = false;

	public static int dnaSubAG = 56;
	public static int dnaSubCT = 53;
	
	public static final int dnaDefaultGamma = 260;
	public static final int dnaDefaultGammaTerm = 100;
	public static final int dnaDefaultLambda = 69;
	public static final int dnaDefaultLambdaTerm = 66;
	
	
	public static final int protDefaultGamma = 60;
	public static final int protDefaultGammaTerm = 16;
	public static final int protDefaultLambda = 38;
	public static final int protDefaultLambdaTerm = 36;
		
	
	public static void initialize (String name, int gamma, int gammaTerm, int lambda, int lambdaTerm) {
	       initialize(name);

	       Aligner.gamma = gamma;
	       Aligner.gammaTerm = gammaTerm;
	       Aligner.lambda = lambda;
	       Aligner.lambdaTerm = lambdaTerm;
	}
	/*
	public static void initialize (String name, String gammas, String gammaTerms, String lambda, String lambdaTerms logger) {
		initialize(name);
		
		Aligner.gamma = Integer.parseInt(gammas);
		Aligner.gammaTerm = Integer.parseInt(gammaTerms);
		Aligner.lambda = Integer.parseInt(lambda);
		Aligner.lambdaTerm = Integer.parseInt(lambdaTerms);
	}*/

	public static void initialize (String name) {
		int [][] tmp_cost = null;
		if (name.equals("BLOSUM62_new")) {
			tmp_cost = setBLOSUM62_new();
		} else if (name.equals("BLOSUM62")){
			tmp_cost = setBLOSUM62();
		} else if (name.equals("BLOSUM50")){
			tmp_cost = setBLOSUM50();
		} else if (name.equals("DNA")){
			tmp_cost = setDNA();
		} else {
			LogWriter.stdErrLogln("You asked me to open the cost file " + name + " ... but the code isn't written for that, yet");
			throw new GenericOpalException("You asked me to open the cost file " + name + " ... but the code isn't written for that, yet");
		}
		setCosts(tmp_cost);
	}
	
	public static int[][] getCosts () {
		return costs;
	}
	
	public static char[] getChars () {
		return chars;
	}

	private static int[][] setBLOSUM62() {
		//60,38,15,36  = gamma, lambda, gamma_term, lambda_term
		chars = "ARNDCQEGHILKMFPSTWYVBZX".toCharArray();
		return getBLOSUM62matrix();
	}
	
	public static int[][] getBLOSUM62matrix() {
		int[][] tmp  = {
				/*A*/		{40},
				/*R*/		{72,32},
				/*N*/		{72,68,32},
				/*D*/		{76,72,56,28},
				/*C*/		{68,84,80,84,12},
				/*Q*/		{68,60,64,64,80,32},
				/*E*/		{68,64,64,56,84,52,36},
				/*G*/		{64,76,68,72,80,76,76,32},
				/*H*/		{72,64,60,72,80,60,64,76,20},
				/*I*/		{72,80,84,84,72,80,84,88,84,40},
				/*L*/		{72,76,84,84,72,76,80,84,80,56,40},
				/*K*/		{68,52,64,68,84,56,60,72,68,80,80,36},
				/*M*/		{68,72,76,84,72,68,76,80,72,56,52,72,32},
				/*F*/		{76,80,80,84,80,84,84,84,72,64,60,84,64,28},
				/*P*/		{68,76,76,72,80,72,72,76,76,80,80,72,80,84,20},
				/*S*/		{56,68,60,64,68,64,64,64,68,80,80,64,72,80,68,40},
				/*T*/		{64,72,64,72,68,68,68,72,76,68,72,68,68,76,72,56,36},
				/*W*/		{80,80,88,88,76,76,80,80,80,80,72,80,72,60,84,80,80,0},
				/*Y*/		{76,76,76,84,80,72,76,84,52,72,72,76,68,48,80,76,72,52,24},
				/*V*/		{64,80,80,84,68,76,80,84,84,48,60,76,60,68,80,72,64,80,72,40},
				/*B*/		{72,72,44,40,84,64,60,68,68,84,84,68,80,84,76,64,68,88,80,84,44},
				/*Z*/		{68,64,64,60,84,44,40,76,64,84,80,60,72,84,72,64,68,80,76,80,64,44},
				/*X*/		{68,72,72,72,76,68,68,72,72,72,72,68,68,72,72,68,68,76,72,68,72,68,72}
				/*            A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X*/
				};
			  
				return tmp;
	}
	

	private static int[][] setBLOSUM62_new() {
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

	
	private static int[][] setBLOSUM50() {
		chars = "ARNDCQEGHILKMFPSTWYVBZX".toCharArray();
		return getBLOSUM50matrix();
	}
	public static int[][] getBLOSUM50matrix() {
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
	private static int[][] setDNA() {
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



	private static void setCosts(int[][] tmp){
		costs = new int[tmp.length][tmp.length];
		for (int i=0; i<tmp.length; i++) 
			for (int j=i; j<tmp.length; j++) 
				costs[i][j] = costs[j][i] = tmp[j][i];	
	}

	public static void increaseCosts(int x) {
		for (int i=0; i<costs.length; i++) 
			for (int j=i; j<costs.length; j++) 
				costs[i][j] = costs[j][i] = costs[j][i] + x; 
	}

	public static void multiplyCosts(float x) {
		for (int i=0; i<costs.length; i++) 
			for (int j=i; j<costs.length; j++) 
				costs[i][j] = costs[j][i] = Math.round(costs[j][i] * x); 
	}
	
}
package opal.align;

import opal.IO.SequenceConverter;
import opal.IO.StructureFileReader;
import opal.IO.Configuration;

public class StructureAlignment extends Alignment {


	public static enum ParamModel { G1, G4, G6, G8};

	
	
	StructureFileReader strucReader;
	public int[][] origSeqIndices ; 
	
	public double[] helixProbSums;
	public double[] sheetProbSums;
	public double[] loopProbSums;
	
	public double[] gapExtSums;

	public float[] gapOpenSums_1; //  sum of the gap-open values for B[j] 
	public float[] gapOpenSums_10; // sum of the gap-open values for B[j], where col j-1 has a char, and col j doesn't
	public float[] gapOpenSums_0; // sum of the gap-open values for the most recent column (<j) with a char, with form f0[j]
	public float[] gapOpenSums_00; // sum of the gap-open values for the most recent column (<j) with a char, with form f00[j]

	
	public StructureAlignment(int[] A, int id, Configuration c ) {
		this(A, id, c, 0);
	}

	public StructureAlignment(int[][] A, int[] ids, boolean isReverse, Configuration c) {
		this(A, ids, isReverse, c, 0);
	}

	public StructureAlignment(int[][] A, int[] ids, Configuration c) {
		this(A, ids, c, 0);
	}

	
	public StructureAlignment(int[] A, int id, Configuration c, int startIndex ) {
		super(A, id, c, startIndex);
	}

	public StructureAlignment(int[][] A, int[] ids, boolean isReverse, Configuration c, int startIndex) {
		super(A, ids, isReverse, c, startIndex);
	}

	public StructureAlignment(int[][] A, int[] ids, Configuration c, int startIndex) {
		super(A, ids, c, startIndex);
	}


	protected final void makeProfile () {
		super.makeProfile();
		
		origSeqIndices = new int[K][M];
		
		helixProbSums = new double[M+1];
		sheetProbSums = new double[M+1];
		loopProbSums = new double[M+1];

		gapExtSums = new double[M+1];
		
		gapOpenSums_1 = new float[M+1];
		gapOpenSums_10 = new float[M+1];
		gapOpenSums_0 = new float[M+1];
		gapOpenSums_00 = new float[M+1];
		
		int[] prevOpenLvl = new int[K];		
		int[] ms = new int[K];
		for (int k=0; k<K; k++) 
			ms[k] = 0;
		int x;
		
		int extGapVal = SequenceConverter.GAP_VAL * 2;
		for (int m=0; m<M; m++) {

			helixProbSums[m+1] = 0;
			sheetProbSums[m+1] = 0;
			loopProbSums[m+1] = 0;

			for (int k=0; k<K; k++) {
				int seqLen = StructureFileReader.helices[seqIds[k]].length ;
				
				if (seqs[k][m] == SequenceConverter.GAP_VAL) {
					if (ms[k]==0 /*external gap*/ || ms[k] == seqLen) {
						origSeqIndices[k][m] = x = extGapVal;
					} else {
						x = SequenceConverter.GAP_VAL;
						origSeqIndices[k][m] = origSeqIndices[k][m-1];
					}
				} else {
					origSeqIndices[k][m] = x = ms[k]++;
				}

				if ( x == SequenceConverter.GAP_VAL ) {
					gapOpenSums_0[m+1] += conf.gapOpenMods[prevOpenLvl[k]];
					if ( seqs[k][m-1] == SequenceConverter.GAP_VAL)  /* x==gapval condition, m-1 is a valid index*/			
						gapOpenSums_00[m+1] += conf.gapOpenMods[prevOpenLvl[k]];
					
				} else if ( x != extGapVal ) {

					int struct_x = x;
					if (isReverse) {
						struct_x = StructureFileReader.helices[seqIds[k]].length - x - 1;
					}

					helixProbSums[m+1] += StructureFileReader.helices[seqIds[k]][struct_x];
					sheetProbSums[m+1] += StructureFileReader.sheets[seqIds[k]][struct_x];
					loopProbSums[m+1] += StructureFileReader.loops[seqIds[k]][struct_x];							
					
					gapExtSums[m+1] += conf.gapExtMods[conf.getStructureLevelFromProbability(StructureFileReader.structureLevels[seqIds[k]][struct_x])];
			
					if (x < seqLen-1) { // at least one more character follows, so opening a gap after this won't give a terminal gap
						if (isReverse) {
							struct_x--;
						}
						prevOpenLvl[k] = conf.getStructureLevelFromProbability(StructureFileReader.structureNeighborLevels[seqIds[k]][struct_x]);
						gapOpenSums_1[m+1] += conf.gapOpenMods[prevOpenLvl[k]];
						if ( seqs[k][m+1] == SequenceConverter.GAP_VAL) 
							gapOpenSums_10[m+2] += conf.gapOpenMods[prevOpenLvl[k]];
					}
				}
			}
		}
	}


	public static void setParams (ParamModel type, Configuration conf) {
		conf.modelType = type;
		
		if (type == ParamModel.G1) {
			conf.gapLevelCnt = 1;
			conf.gamma = 62;
			conf.setGammaTerm(27);
			int[] open = {0};
			conf.gapOpenMods = open;
			
			conf.lambda = 37;
			conf.setLambdaTerm(34);
			int[] ext = {0};
			conf.gapExtMods = ext;
			
		} else if (type == ParamModel.G4) {
			conf.gapLevelCnt = 4;
			conf.gamma = 51;
			conf.setGammaTerm(23);
			int[] open = {0, 9, 32, 32};
			conf.gapOpenMods = open;
			
			conf.lambda = 36;
			conf.setLambdaTerm(36);
			int[] ext = {0, 2, 4, 4};
			conf.gapExtMods = ext;
			
		} else if (type == ParamModel.G6) {
			conf.gapLevelCnt = 6;
			conf.gamma = 45;
			conf.setGammaTerm(11);
			int[] open = {0, 6, 19, 30, 49, 49};
			conf.gapOpenMods = open;
			
			conf.lambda = 36;
			conf.setLambdaTerm(35);
			int[] ext = {0, 1, 2, 2, 2, 2 };
			conf.gapExtMods = ext;	

		} else if (type == ParamModel.G8) {
			conf.gapLevelCnt = 8;
			conf.gamma = 47;
			conf.setGammaTerm(19);
			int[] open = {0, 0, 12, 19, 33, 36, 57, 57};
			conf.gapOpenMods = open;
			
			conf.lambda = 38;
			conf.setLambdaTerm(37);
			int[] ext = {0, 0, 1, 1, 1, 2, 4, 4};
			conf.gapExtMods = ext;		
		}
		
	}
	
}

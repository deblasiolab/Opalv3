package opal.align;

import com.traviswheeler.libs.LogWriter;
import opal.IO.SequenceConverter;
import opal.IO.StructureFileReader;

public class StructureAlignment extends Alignment {

	public static int gapLevelCnt = 1;

	public static enum ParamModel {firstGuess, UN1, UN4, UN6, UN8};

	
	public static ParamModel modelType;
	
	public static int subHelixHelix;
	public static int subHelixSheet;
	public static int subSheetSheet;
	public static int subHelixLoop;
	public static int subSheetLoop;
	public static int subLoopLoop;
 	public static int[] gapOpenMods;
	public static int[] gapExtMods;
	
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

	
	public StructureAlignment(int[] A, int id ) {
		super(A, id);
	}

	public StructureAlignment(int[][] A, int[] ids, boolean isReverse) {
		super(A, ids, isReverse);
	}

	public StructureAlignment(int[][] A, int[] ids) {
		super(A, ids);
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
					gapOpenSums_0[m+1] += gapOpenMods[prevOpenLvl[k]];
					if ( seqs[k][m-1] == SequenceConverter.GAP_VAL)  /* x==gapval condition, m-1 is a valid index*/			
						gapOpenSums_00[m+1] += gapOpenMods[prevOpenLvl[k]];
					
				} else if ( x != extGapVal ) {

					int struct_x = x;
					if (isReverse) {
						struct_x = StructureFileReader.helices[seqIds[k]].length - x - 1;
					}

					helixProbSums[m+1] += StructureFileReader.helices[seqIds[k]][struct_x];
					sheetProbSums[m+1] += StructureFileReader.sheets[seqIds[k]][struct_x];
					loopProbSums[m+1] += StructureFileReader.loops[seqIds[k]][struct_x];							
					
					gapExtSums[m+1] += gapExtMods[StructureFileReader.structureLevels[seqIds[k]][struct_x]];
			
					if (x < seqLen-1) { // at least one more character follows, so opening a gap after this won't give a terminal gap
						if (isReverse) {
							struct_x--;
						}
						prevOpenLvl[k] = StructureFileReader.structureNeighborLevels[seqIds[k]][struct_x];
						gapOpenSums_1[m+1] += gapOpenMods[prevOpenLvl[k]];
						if ( seqs[k][m+1] == SequenceConverter.GAP_VAL) 
							gapOpenSums_10[m+2] += gapOpenMods[prevOpenLvl[k]];
					}
				}
			}
		}
	}


	public static void setParams (ParamModel type) {
		modelType = type;
		
		if (type == ParamModel.firstGuess) {
			gapLevelCnt = 1;
			Aligner.gamma = 60;
			Aligner.gammaTerm = 16;
			int[] open = {0};
			gapOpenMods = open;
			
			Aligner.lambda = 38;
			Aligner.lambdaTerm = 36;
			int[] ext = {0};
			gapExtMods = ext;

			subHelixHelix = -10;
			subHelixSheet = 10;
			subSheetSheet = -20;
			subHelixLoop = 5;
			subSheetLoop = 5 ;
			subLoopLoop = 0;
		} else if (type == ParamModel.UN1) {
			gapLevelCnt = 1;
			Aligner.gamma = 62;
			Aligner.gammaTerm = 27;
			int[] open = {0};
			gapOpenMods = open;
			
			Aligner.lambda = 37;
			Aligner.lambdaTerm = 34;
			int[] ext = {0};
			gapExtMods = ext;

			subHelixHelix = -20;
			subHelixSheet = 0;
			subSheetSheet = -45;
			subHelixLoop = 9;
			subSheetLoop = 5 ;
			subLoopLoop = 0;			
		} else if (type == ParamModel.UN4) {
			gapLevelCnt = 4;
			Aligner.gamma = 51;
			Aligner.gammaTerm = 23;
			int[] open = {0, 9, 32, 32};
			gapOpenMods = open;
			
			Aligner.lambda = 36;
			Aligner.lambdaTerm = 36;
			int[] ext = {0, 2, 4, 4};
			gapExtMods = ext;

			subHelixHelix = -15;
			subHelixSheet = 0;
			subSheetSheet = -41;
			subHelixLoop = 14;
			subSheetLoop = 10;
			subLoopLoop = 0;
		} else if (type == ParamModel.UN6) {
			gapLevelCnt = 6;
			Aligner.gamma = 45;
			Aligner.gammaTerm = 11;
			int[] open = {0, 6, 19, 30, 49, 49};
			gapOpenMods = open;
			
			Aligner.lambda = 36;
			Aligner.lambdaTerm = 35;
			int[] ext = {0, 1, 2, 2, 2, 2 };
			gapExtMods = ext;

			subHelixHelix = -11;
			subHelixSheet = 0;
			subSheetSheet = -48;
			subHelixLoop = 10;
			subSheetLoop = 9;
			subLoopLoop = 0;			

		} else if (type == ParamModel.UN8) {
			gapLevelCnt = 8;
			Aligner.gamma = 47;
			Aligner.gammaTerm = 19;
			int[] open = {0, 0, 12, 19, 33, 36, 57, 57};
			gapOpenMods = open;
			
			Aligner.lambda = 38;
			Aligner.lambdaTerm = 37;
			int[] ext = {0, 0, 1, 1, 1, 2, 4, 4};
			gapExtMods = ext;

			subHelixHelix = -12;
			subHelixSheet = 0;
			subSheetSheet = -46;
			subHelixLoop = 15;
			subSheetLoop = 14;
			subLoopLoop = 0;			
		}
		
	}
	
}

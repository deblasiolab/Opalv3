package opal.align;

import java.util.ArrayList;

import opal.exceptions.GenericOpalException;
import opal.tree.Tree;
import opal.tree.TreeNode;

import opal.IO.CostMatrix;
import com.traviswheeler.libs.LogWriter; 
import opal.IO.SequenceConverter;
import opal.IO.StructureFileReader;

public abstract class Aligner  {

	public static enum AlignmentType {exact, profile, mixed};
	public static enum Direction {vert, horiz, diag };

	public static AlignmentType alignmentMethod = AlignmentType.mixed;
	
	public static int mixedAlignmentCutoff = 40;
	public static int linearCutoff = 100;
	public static int gamma;
	public static int gammaTerm;
	public static int lambda;
	public static int lambdaTerm;
	public static SequenceConverter seqConv;
	public static int[][] costs;

	public static int[] gammas = null;
	public static int[] lambdas = null;
	public static int[] gammaTerms = null;
	public static int[] lambdaTerms = null;
	
	
	public static boolean useStructure = false;
	
	//vars in the standard class
	int[][] C;
//	char[] alphabet;
	int K,L,M,N; //num seqs in A,B ... #cols in A,B
	long[][] H, V, D;
	ArrayList<Direction> path;
	public Alignment A = null, B = null;
	StructureAlignment structA = null, structB = null;

	protected Alignment resultAlignment = null;
	protected TreeNode node = null;
	
	//heuristic vars
	boolean isPessimistic = true;
	long estimatedCost = -1;
	long trueCost = -1;

	
	public static void setParams () {
		costs = CostMatrix.getCosts();
	}

	
	public static void setParams (String gammaVals, String gammaTermVals, String lambdaVals, String lambdaTermVals) {
//not used
		String[] tmp = gammaVals.split(",");
		gammas = new int[tmp.length]; //first one is already the default value
		for(int i=0; i<tmp.length; i++) {
			gammas[i] = Integer.parseInt(tmp[i]);
		}
		tmp = gammaTermVals.split(",");
		gammaTerms = new int[tmp.length]; 
		for(int i=0; i<tmp.length; i++) {
			gammaTerms[i] = Integer.parseInt(tmp[i]);
		}	
		tmp = lambdaVals.split(",");
		lambdas = new int[tmp.length]; 
		for(int i=0; i<tmp.length; i++) {
			lambdas[i] = Integer.parseInt(tmp[i]);
		}	
		tmp = lambdaTermVals.split(",");
		lambdaTerms = new int[tmp.length];
		for(int i=0; i<tmp.length; i++) {
			lambdaTerms[i] = Integer.parseInt(tmp[i]);
		}
		
		switchParamID(0);
		costs = CostMatrix.getCosts();
	}
	
	public static void switchParamID (int i) {
	//not used
		gamma = gammas[i];
		gammaTerm = gammaTerms[i];
		lambda = lambdas[i];		
		lambdaTerm = lambdaTerms[i];		
	}
	
	public static void increaseGapCosts(int increase) {
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

	public static void multiplyGapCosts(float x) {
		lambda = Math.round(lambda * x);
		lambdaTerm = Math.round(lambdaTerm * x);
		gamma = Math.round(gamma * x);
		gammaTerm = Math.round(gammaTerm * x);
		
		if (lambdas != null){
			for(int i=0; i<lambdas.length; i++) {
				lambdas[i] = Math.round(lambdas[i] * x);
				lambdaTerms[i] = Math.round(lambdaTerms[i] * x);
				gammas[i] = Math.round(gammas[i] * x);
				gammaTerms[i] = Math.round(gammaTerms[i] * x);
			}
		}
		
	}	
	public static void setSequenceConverter (SequenceConverter sc) {
		seqConv = sc;	
	}

	public Aligner( Aligner al) {
		this(al.A, al.B, al.isPessimistic);	
		this.node = al.node;
	}
	
	public Aligner( Alignment A, Alignment B) {
		this(A,B,true);
	}

	public Aligner() {
		this(true);
	}
	
	public Aligner( Alignment A, Alignment B, boolean pess) {
		this (pess);
		setAlignments(A, B);		
	}
	
	public Aligner(  boolean pess) {
		isPessimistic = pess;
	}

	
	public void setNode (TreeNode n) {
		node = n;
	}
	
	public void setAlignments (Alignment A, Alignment B) {
		this.A = A;
		this.B = B;
		if (A != null) {
			K = A.K;
			M = A.seqs[0].length;
		}
		if (B != null) {
			L = B.K;
			N = B.seqs[0].length;
		}

		if (A instanceof StructureAlignment && B instanceof StructureAlignment) {
			structA = (StructureAlignment)A;
			structB = (StructureAlignment)B;			
		}
		
	}
	
	abstract protected void initialize ();
	abstract protected void fillBoundary ();
	abstract protected void fillTable ();
	abstract protected boolean recover ();
	abstract public void cleanup ();
	
	public void align () {
		resultAlignment = null;
		trueCost = -1;
		estimatedCost = -1;

		initialize();
		fillBoundary();
		fillTable();
		if (!recover()) {
			LogWriter.stdErrLogln("problem generating alignment. Quitting.");
			throw new GenericOpalException("Opal: problem generating alignment. Quitting.");
		}
	}
	
	public ArrayList<Direction> getPath () {
		return path;
	}
	
	public long getEstimatedCost () {
		return estimatedCost;
	}
	public long getTrueCost () {
		if (trueCost > -1) {
			return trueCost;
		} else {
			Alignment al = getAlignment();
			trueCost = Aligner.calcCost(al.seqs, al.seqIds);
			return trueCost;
		}
	}
	
	public void setPessimistic (boolean p){
		isPessimistic = p;
	}
	

	public static long calcCost(char[][] alignment, int[] ids ) {
		return calcCost(alignment, alignment.length, 0, ids, true);
	}
	public static long calcCost(char[][] alignment, int K, int L, int[] ids ) {
		return calcCost(alignment, K, L, ids, false);
	}
	private static long calcCost(char[][] alignment, int K, int L, int[] ids, boolean allRows) {
		int[][]C = SequenceConverter.convertSeqsToInts(alignment);
		return calcCost(C, K, L, ids, allRows );
	}	
	
	public static long calcCost(int[][] alignment, int[] ids ) {
		return calcCost(alignment, alignment.length, 0, ids, true);
	}

	public static long calcCost(int[][] alignment, int K, int L, int[] ids ) {
		return calcCost(alignment, K, L, ids, false);
	}

		
	protected static long calcCost(int[][] C, int K, int L, int[] ids, boolean allRows) {
		int len = C[0].length;	
		double[] colCosts = new double[len];
		
		int a,b;
		Direction inGap = null;
		
		int xPos, yPos;
		
		for (int x=0; x < (allRows ? K+L-1 : K); x++) {
			for (int y=(allRows ? x+1 : K) ; y<K+L; y++) {
				//long cc = 0, cc2 = 0;
				inGap = null;
				xPos = yPos = -1;
				
				int i=0;
				//find terminal gap on left - skip the initial spaces
				while (SequenceConverter.GAP_VAL == C[x][i] && SequenceConverter.GAP_VAL == C[y][i]) i++;
				if (SequenceConverter.GAP_VAL == C[x][i] && SequenceConverter.GAP_VAL != C[y][i]) {
					//cc2 += gammaTerm;
					colCosts[i] += gammaTerm;
					while (SequenceConverter.GAP_VAL == C[x][i]) {
						if (SequenceConverter.GAP_VAL != C[y][i]) {
							//cc2 += lambdaTerm;
							colCosts[i] += lambdaTerm;
							yPos++;
						}
						i++;
					}
				} else if (SequenceConverter.GAP_VAL != C[x][i] && SequenceConverter.GAP_VAL == C[y][i]) {
					//cc2 += gammaTerm;
					colCosts[i] += gammaTerm;
					while ( SequenceConverter.GAP_VAL == C[y][i]) {
						if (SequenceConverter.GAP_VAL != C[x][i]){
							//cc2 += lambdaTerm;
							colCosts[i] += lambdaTerm;
							xPos++;
						}
						i++;
					}					
				}
			
				
				int beg = i;
				i = len-1;
				//find terminal gap on right - skip the initial spaces
				while (SequenceConverter.GAP_VAL == C[x][i] && SequenceConverter.GAP_VAL == C[y][i]) i--;
				if (SequenceConverter.GAP_VAL == C[x][i] && SequenceConverter.GAP_VAL != C[y][i]) {
					//cc2 += gammaTerm;
					colCosts[i] += gammaTerm;
					while (SequenceConverter.GAP_VAL == C[x][i]) {
						if (SequenceConverter.GAP_VAL != C[y][i]){
							//cc2 += lambdaTerm;
							colCosts[i] += lambdaTerm;
						}
						i--;
					}
				} else if (SequenceConverter.GAP_VAL != C[x][i] && SequenceConverter.GAP_VAL == C[y][i]) {
					//cc2 += gammaTerm;
					colCosts[i] += gammaTerm;
					while ( SequenceConverter.GAP_VAL == C[y][i]) {
						if (SequenceConverter.GAP_VAL != C[x][i]){
							//cc2 += lambdaTerm;
							colCosts[i] += lambdaTerm;
						}
						i--;
					}					
				}
				int end = i;

				
				for (i=beg; i<=end; i++) {
					a = C[x][i];
					b = C[y][i];
									
					
					float cc = 0;					
					if (a != SequenceConverter.GAP_VAL && b != SequenceConverter.GAP_VAL) { // substitution\
						cc += costs[a][b];
						xPos++;
						yPos++;
						if (useStructure) cc += getStructSubModifierPair(ids[x],ids[y],xPos,yPos);
						
						inGap = null;
					} else if (a != SequenceConverter.GAP_VAL) { //gap in B, i.e. VERT
						xPos++;
						cc += lambda;
						if (useStructure) cc += getStructGapExtModiferPair (ids[x],xPos);
						if (inGap != Direction.vert) { 
							cc += gamma; 
							inGap = Direction.vert;	
							if (useStructure) cc += getStructGapOpenModiferPair (ids[y],yPos);
						}						
					} else if (b != SequenceConverter.GAP_VAL) { //gap in A, i.e. HORIZ
						yPos++;
						cc += lambda;
						if (useStructure) cc += getStructGapExtModiferPair (ids[y],yPos);
						
						if (inGap != Direction.horiz) { 
							cc += gamma; 
							inGap = Direction.horiz;
							if (useStructure) cc += getStructGapOpenModiferPair (ids[x],xPos);
						}
					} // otherwise, it's a double-dash column 

					colCosts[i] += cc;
					//Logger.stdErrLogln("col " + i + " : " + cc);
				}
			}
		
		}	
		
		long cost = 0;
//		int x=1;
		for (double c: colCosts) {
			cost += Math.round( c + .00000005);
//			Logger.stdErrLogln("col " + (x++) + " : " + Math.round(c + .00000005));
		}		
		return cost;
	}	

	public Alignment getAlignment () {

		if (null == resultAlignment) {
			int[][] result = SequenceConverter.convertPathToIntAlignment(path, A, B);
			
			int[] ids = new int[A.seqIds.length + B.seqIds.length];
			for (int i=0; i<A.seqIds.length; i++) ids[i] = A.seqIds[i];
			for (int i=0; i<B.seqIds.length; i++) ids[A.seqIds.length+i] = B.seqIds[i];
			resultAlignment = Alignment.buildNewAlignment(result, ids);
		}

		if (null != node)
			node.setAlignment(resultAlignment);	

		return resultAlignment;
	}
	
	
	public void preprocessTree (Tree tree, float distances[][]) {
		TreeNode root = tree.getRoots().get(0);
		root.assignLeafRanges(); // this tells me first leaf and range of leaves, which gets me a translation to the input seqs
	}
	
	final long min (long x, long y, long z) {
		return x<y ? x<z?x:z  :  y<z?y:z ;
	}	
		
	public static long getStructSubModifier (StructureAlignment structA, StructureAlignment structB, int i, int j) {
		
		double x=0;
		
		double hA = structA.helixProbSums[i];
		double hB = structB.helixProbSums[j];
		double sA = structA.sheetProbSums[i];
		double sB = structB.sheetProbSums[j];
		double lA = structA.loopProbSums[i];
		double lB = structB.loopProbSums[j];

		//limit the significant digits
		double hh = Math.round(hA * hB * 1000000) / 1000000.0;   
		double hs = Math.round(hA * sB * 1000000) / 1000000.0;   
		double hl = Math.round(hA * lB * 1000000) / 1000000.0;   
		
		double sh = Math.round(sA * hB * 1000000) / 1000000.0;   
		double ss = Math.round(sA * sB * 1000000) / 1000000.0;   
		double sl = Math.round(sA * lB * 1000000) / 1000000.0;   

		double lh = Math.round(lA * hB * 1000000) / 1000000.0;   
		double ls = Math.round(lA * sB * 1000000) / 1000000.0;   
		double ll = Math.round(lA * lB * 1000000) / 1000000.0;   

			
		x= StructureAlignment.subHelixHelix * hh; 	
		x += StructureAlignment.subSheetSheet * ss;
		x += StructureAlignment.subLoopLoop * ll;

		x += StructureAlignment.subHelixLoop * (hl + lh);
		x += StructureAlignment.subHelixSheet * (hs + sh);
		x += StructureAlignment.subSheetLoop * (sl + ls);
		return Math.round(x + .00000005);
	}
	
/*
	protected static long getStructGapModifer (StructureAlignment structA, StructureAlignment structB, int i, int j, boolean isSub, boolean isPessimistic) {
		
		return getStructGapOpenModifer(structA, structB, i, j, isSub, isPessimistic) +
				getStructGapExtModifer(structA, structB, i, j, isSub) ;
	}
	*/
	
	protected static long getStructGapOpenModifer (StructureAlignment structA, StructureAlignment structB, int i, int j, Direction curDir, Direction prevDir, boolean isPessimistic ) {
		
		double mod = 0;
		
		if (curDir == Direction.vert) {
			if (prevDir == Direction.vert) {
				if (isPessimistic) {
					mod += structA.f01[i] * ( structB.gapOpenSums_0[j] + structB.gapOpenSums_1[j]) ;
				}
			} else if (prevDir == Direction.horiz) {
				mod += structA.f1[i] * structB.gapOpenSums_1[j] ;
				if (isPessimistic) {
					mod += structA.f1[i] * structB.gapOpenSums_0[j] ;
				}
			} else if (prevDir == Direction.diag) {
				mod += structA.f1[i] * structB.gapOpenSums_1[j] ;
				if (isPessimistic) {
					mod += structA.f01[i] * structB.gapOpenSums_0[j] ;
				}
			}
		} else if (curDir == Direction.horiz) {
			if (prevDir == Direction.vert) {
				mod += structB.f1[j] * structA.gapOpenSums_1[i] ;
				if (isPessimistic) {
					mod += structB.f1[j] * structA.gapOpenSums_0[i] ;
				}
			} else if (prevDir == Direction.horiz) {
				if (isPessimistic) {
					mod += structB.f01[j] * ( structA.gapOpenSums_0[i] + structA.gapOpenSums_1[i]) ;
				}
			} else if (prevDir == Direction.diag) {
				mod += structB.f1[j] * structA.gapOpenSums_1[i] ;
				if (isPessimistic) {
					mod += structB.f01[j] * structA.gapOpenSums_0[i] ;
				}
			}
		} else if (curDir == Direction.diag /*substitution*/) {
			if (prevDir == Direction.vert) {
				mod += structB.f1[j] * structA.gapOpenSums_10[i] ;
				if (isPessimistic) {
					mod += structB.f1[j] * structA.gapOpenSums_00[i] ;
					mod += structA.f01[i] * structB.gapOpenSums_0[j] ;
				}	
			} else if (prevDir == Direction.horiz) {
				mod += structA.f1[i] * structB.gapOpenSums_10[j] ;
				if (isPessimistic) {
					mod += structA.f1[i] * structB.gapOpenSums_00[j] ;
					mod += structB.f01[j] * structA.gapOpenSums_0[i] ;
				}
			} else if (prevDir == Direction.diag) {
				mod += structA.f1[i] * structB.gapOpenSums_10[j] ;
				mod += structB.f1[j] * structA.gapOpenSums_10[i] ;
				if (isPessimistic) {
					mod += structA.f01[i] * structB.gapOpenSums_00[j] ;
					mod += structB.f01[j] * structA.gapOpenSums_00[i] ;
				}	
			}
		}
		
		return Math.round(mod + .00000005);
	}
	
	
	public static long getStructGapExtModifer (StructureAlignment structA, StructureAlignment structB, int i, int j, Direction curDir) {
		/*e.g.
		 * A   x
		 * B   -
		 */
		
		double mod = 0;
		
		if (curDir == Direction.vert) {
			if (j>0)
				mod = structA.gapExtSums[i] * (structB.K - structB.gapsBeforeFirst[j] - structB.gapsAfterLast[j] - structB.lastLetterCount[j]); 
		} else if (curDir == Direction.horiz) {
			if (i>0)
				mod = structB.gapExtSums[j] * (structA.K - structA.gapsBeforeFirst[i] - structA.gapsAfterLast[i] - structA.lastLetterCount[i]); 
		} else {
			mod = structA.gapExtSums[i] * (structB.f0[j] - structB.gapsBeforeFirst[j] - structB.gapsAfterLast[j] );
			mod += structB.gapExtSums[j] * (structA.f0[i] - structA.gapsBeforeFirst[i] - structA.gapsAfterLast[i] );
		}
		return Math.round(mod + .00000005);
	}

	
	protected static double getStructSubModifierPair (int k, int l, int i, int j) {
		
		double hA = StructureFileReader.helices[k][i];
		double hB = StructureFileReader.helices[l][j];
		double sA = StructureFileReader.sheets[k][i];
		double sB = StructureFileReader.sheets[l][j];
		double lA = StructureFileReader.loops[k][i];
		double lB = StructureFileReader.loops[l][j];

		
		//limit the significant digits
		double hh = Math.round(hA * hB * 1000000) / 1000000.0;   
		double hs = Math.round(hA * sB * 1000000) / 1000000.0;   
		double hl = Math.round(hA * lB * 1000000) / 1000000.0;   
		
		double sh = Math.round(sA * hB * 1000000) / 1000000.0;   
		double ss = Math.round(sA * sB * 1000000) / 1000000.0;   
		double sl = Math.round(sA * lB * 1000000) / 1000000.0;   

		double lh = Math.round(lA * hB * 1000000) / 1000000.0;   
		double ls = Math.round(lA * sB * 1000000) / 1000000.0;   
		double ll = Math.round(lA * lB * 1000000) / 1000000.0;   

		
		double x = StructureAlignment.subHelixHelix * hh; 	
		x += StructureAlignment.subSheetSheet * ss;
		x += StructureAlignment.subLoopLoop * ll;

		x += StructureAlignment.subHelixLoop * (hl + lh);
		x += StructureAlignment.subHelixSheet * (hs + sh);
		x += StructureAlignment.subSheetLoop * (sl + ls);
		return x;
	}

	protected static double getStructGapOpenModiferPair (int seqid, int i) {
		/*
		 * A    x
		 * B  x -      seqid
 		 *    i
		 */
			
		return StructureAlignment.gapOpenMods[ StructureFileReader.structureNeighborLevels[seqid][i] ] ;
	}

	protected static double getStructGapExtModiferPair (int seqid, int i) {
		/*     
		 *     i
		 * A   x     seqid
		 * B   -
		 */
		
		return StructureAlignment.gapExtMods[ StructureFileReader.structureLevels[seqid][i] ] ;
	}

}

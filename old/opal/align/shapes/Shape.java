package opal.align.shapes;

import opal.align.Aligner;
import opal.align.Alignment;
import opal.align.StructureAlignment;
import com.traviswheeler.libs.LogWriter;
import opal.IO.SequenceConverter;

abstract public class Shape {
	protected int[] seqBlocks; //for each seq, the block it belongs to
	Shape parent = null;
	long shapeCost = 0;
	int maxBlock = 0;
	int aPos=0;
	int bPos=0;

	static protected Alignment A, B;
	static int K, L, M, N ;
	static Alignment alignment;
	static int[][] costs;
	
	
	// create a new flat shape
	public Shape () {
		seqBlocks = new int[K+L];
		for (int i=0; i<K+L; i++) seqBlocks[i] = 0;
	}
	
	//make a copy of the input shape, and set the original as the parent of the copy
	public Shape (Shape s) {
		seqBlocks = (int[])s.seqBlocks.clone();
		maxBlock = s.maxBlock;
		shapeCost = s.shapeCost;
		aPos = s.aPos;
		bPos = s.bPos;
		
		parent = s;
	}
	
	/*
	 * Just a way of reducing clutter in the Shape-instantiating code ... 
	 * let the static method figure out whether we're dealing with consistency shapes, 
	 * or whether the input is big enough to 
	 * bother using the linear methods (with accompanying overhead) 
	 */
	final public static Shape makeNewShape(Shape s) {
		if (s instanceof ConsistencyShape) 
			return new ConsistencyShape(s);
		else if (Shape.K * Shape.L > Aligner.linearCutoff) 
			return new ShapeLinear(s);
		else 
			return new ShapeQuadratic(s);
	}
	
	public static void setAlignments (Alignment alA, Alignment alB) {
		A = alA;
		B = alB;		
		K = alA.K;
		L = alB.K;
		M = alA.M;
		N = alB.M;
	}
		
	final public static void setParams(int[][] cost_arr){
		costs = cost_arr;
	}
		
			
	/* modify shape based on adding column(s) from A or B (or both) to shape s
	 * input:
	 * 		Shape s that this shape is based on
	 * 		a and b (column indices into A and B respectively) for columns to be added
	 * if a or b are <0, that means no column is appended for that alignment (gap)
	 */
	public void appendColumns (int a, int b) {
		int i=0;
		if (a<0 && b<0) return;

		Aligner.Direction dir = Aligner.Direction.horiz;
		if (a>0 )  // if not, then keep the default
			if (b>0) 
				dir = Aligner.Direction.diag;
			else
				dir = Aligner.Direction.vert;

		long boundaryCost = gapBoundaryCost(a,b,dir); // this is calculated before updating the shape, because it relies on the old shape with the new column(s)
		
		//determine upper portion of new shape 
		if (a>0) {
			aPos = a;
			for (i = 0; i < K; i++) {
				if (A.seqs[i][a-1] != SequenceConverter.GAP_VAL)
					seqBlocks[i] = maxBlock+1;
			}
		}
		//determine lower portion of new shape 
		if (b > 0) {
			bPos = b;
			for ( i = 0; i < L; i++) {
				if (B.seqs[i][b-1] != SequenceConverter.GAP_VAL)        
					seqBlocks[i+K] = maxBlock+1;  
			}
		} 
		
		maxBlock++;
		//compress block range   set flags (cnt) for blocks present in shape
		int[] cnt = new int[maxBlock+1];
		for (i = 0; i < K+L; i++) 
			cnt[seqBlocks[i]]++;

//		modifyOtherBlockLists (a, b, cnt);
		
		int[] cnv = new int[maxBlock+1];
		int j = -1;
		for (i = 0; i <= maxBlock; i++) {
		    if (cnt[i] > 0) { j++; }  //fill in the conversion array with a unique  
		    cnv[i] = j;               //block number for each block in shape        
		}	
		maxBlock = cnv[maxBlock];
		for (i = 0; i < K+L; i++) { 
			seqBlocks[i] = cnv[seqBlocks[i]]; //  set new blocks
		}
		
				
		//modify cost here  using gap/subs calcs.
		shapeCost += subCost(a,b) + boundaryCost; 		
		
	}

	public void freeUnusedStructures () {
		seqBlocks = null;
	}
		
//	protected void modifyOtherBlockLists (int a, int b, int[] cnt) {
		// nothing done here. This allows BurialShape to modify 
		// info that tracks columns for each block		
//	}

	/* For a column with p different characters in column a of A, 
	 * and q different characters in col b of B, this method takes
	 * O(pq) time for the diagonal calculation, and constant for the other two.
	 */
	protected long subCost (int a, int b) {
		long cost = 0;
	
		StructureAlignment structA = null;
		StructureAlignment structB = null;		
  
		if (Aligner.useStructure) {
			structA = (StructureAlignment)A;
			structB = (StructureAlignment)B;		
		}
		
		if (a<0) {//horiz	
			int termA = (aPos == M || aPos == 0) ? K : A.gapsBeforeFirst[aPos] + A.gapsAfterLast[aPos] + A.lastLetterCount[aPos];
			cost = B.f1[b] * (Aligner.lambda * (K-termA)  +  Aligner.lambdaTerm * termA );
			if (Aligner.useStructure) {
				cost += Aligner.getStructGapExtModifer(structA, structB, aPos, bPos, Aligner.Direction.horiz);
			}

		} else if (b<0) { //vert
			int termB = (bPos == N || bPos == 0) ? L : B.gapsBeforeFirst[bPos] + B.gapsAfterLast[bPos] + B.lastLetterCount[bPos];
			cost = A.f1[a] * (Aligner.lambda * (L-termB) +  Aligner.lambdaTerm * termB );
			if (Aligner.useStructure) {
				cost += Aligner.getStructGapExtModifer(structA, structB, aPos, bPos, Aligner.Direction.vert);
			}

		} else { //dir == DIAG
			int termA = A.gapsBeforeFirst[a]+A.gapsAfterLast[a];
			int termB = B.gapsBeforeFirst[b]+B.gapsAfterLast[b];
			// all substitution letter combos, including gap extensions
			int internalExtensionCount = (A.f0[a]-termA) * B.f1[b] + A.f1[a] * (B.f0[b]-termB) ;
			cost = Aligner.lambda * ( internalExtensionCount ) + Aligner.lambdaTerm * (termA * B.f1[b] + termB * A.f1[a]);
			
			for (int x=0; x<A.chars[a].length; x++){
				for (int y=0; y<B.chars[b].length; y++){
					int s = costs[A.chars[a][x]][B.chars[b][y]]; 
					cost += A.freqs[a][x]*B.freqs[b][y]*s;
				}	
			}	
			
			if (Aligner.useStructure) {
				cost += Aligner.getStructSubModifier(structA, structB, a,b);
				cost += Aligner.getStructGapExtModifer(structA, structB, aPos, bPos, Aligner.Direction.diag);
			}
			
		}
		return cost;
	}


	abstract protected long gapBoundaryCost (int a, int b, Aligner.Direction dir) ;
	
//	abstract public int countGapOpenOvercounts(Shape s);
	
	
	public long getCost () {
		return shapeCost;
	}
	public int getMaxBlock () {
		return maxBlock;
	}	
	public void setParent (Shape s) {
		parent = s;
	}
	public Shape getParent() {
		return parent;
	}
	
	public void printInfo () {
		LogWriter.stdOutLogln(aPos + ", " + bPos + " : " + shapeCost);
	}
	
	final public Aligner.Direction getParentDirection() {
		if (parent == null) 
			return null;
		else if (parent.aPos == aPos && parent.bPos == bPos-1) 
			return Aligner.Direction.horiz;
		else if (parent.aPos == aPos-1 && parent.bPos == bPos)
			return Aligner.Direction.vert;
		else if (parent.aPos == aPos-1 && parent.bPos == bPos-1)
			return Aligner.Direction.diag;
		else
			return null;
	}

}


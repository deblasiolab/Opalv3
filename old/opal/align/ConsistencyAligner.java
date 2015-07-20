package opal.align;

import java.util.ArrayList;
import java.util.Iterator;

import com.traviswheeler.libs.LogWriter;
import opal.IO.SequenceConverter;
import opal.align.shapes.ConsistencyShape;
import opal.align.shapes.ConsistencyShapeTester;
import opal.align.shapes.Shape;
import opal.exceptions.GenericOpalException;
import opal.tree.Tree;
import opal.tree.TreeNode;


public class ConsistencyAligner extends ExactCountAligner_Time {
	
	
	public static enum Direction {vert, horiz, diag };

	public static Aligner.AlignmentType alignmentMethod = AlignmentType.exact;

	boolean testing = false;
//	boolean testing = true;
	
	public Aligner simpleAligner;
	PairwiseAlignmentsContainer alignmentContainer;
	ConsistencyModifiers_AllPairs mods;
	
	
	public ConsistencyAligner(Aligner al) {
		this(true);
		simpleAligner = al;
	}

	public ConsistencyAligner() {
		this(true);
	}
	
	public ConsistencyAligner( Alignment A, Alignment B) {
		super(A,B);
		LogWriter.stdErrLogln("Not expecting to be used that way...");
		throw new GenericOpalException("");
	}
	
	public ConsistencyAligner( boolean pess) {
		super(pess);		
	}
	
	public ConsistencyAligner( Alignment A, Alignment B, boolean pess) {
		super(A,B,pess);
		LogWriter.stdErrLogln("Not expecting to be used that way...");
		throw new GenericOpalException("");
	}
	
	public void setAlignments (Alignment A, Alignment B) {
		super.setAlignments(A,B);
		simpleAligner.setAlignments(A, B);
	}

	public void setNode (TreeNode n) {
		super.setNode(n);
		simpleAligner.setNode(n);
	}

	
	public void align () {
		
		//int nodeCount = node.lastLeaf - node.firstLeaf + 1;
		int nodeCount = A.K + B.K;
		if ( nodeCount > PairwiseAlignmentsContainer.maxSubtreeSize) {
			simpleAligner.align();
	//		set vars to match those results
		//} else if (nodeCount == 2 || simpleAligner instanceof ProfileAligner) {
		} else if ( (A.K==1 && B.K==1) || alignmentMethod == AlignmentType.profile) {
			mods = new ConsistencyModifiers_AllPairs(node, alignmentContainer);
			ConsistencyHeuristicAligner al = new ConsistencyHeuristicAligner(this);
			al.align();
			path = al.getPath();
			estimatedCost = al.getEstimatedCost();
//			char[][] result = SequenceConverter.convertPathToCharAlignment(path, SequenceConverter.convertIntsToSeqs(A.seqs), SequenceConverter.convertIntsToSeqs(B.seqs));
//			int[][] result = SequenceConverter.convertPathToCharAlignment(al.getPath(), A, B);
		} else {
			super.align();
		}
	}

	public Alignment getAlignment () {
		int nodeCount = node.lastLeaf - node.firstLeaf + 1;
		if ( nodeCount > PairwiseAlignmentsContainer.maxSubtreeSize) {
			return simpleAligner.getAlignment();
		} else {
			return super.getAlignment();
		}
	}
		
	public void preprocessTree (Tree tree, float distances[][]) {
		super.preprocessTree(tree, distances);
		
		if (PairwiseAlignmentsContainer.blendMethod == PairwiseAlignmentsContainer.BlendType.symmetric ) {
			alignmentContainer = new PairwiseAlignmentContainer_symmetricBlend(tree, distances);
		} else if (PairwiseAlignmentsContainer.blendMethod == PairwiseAlignmentsContainer.BlendType.simple) {
			alignmentContainer = new PairwiseAlignmentContainer_simpleBlend(tree, distances);
		} else if (PairwiseAlignmentsContainer.blendMethod == PairwiseAlignmentsContainer.BlendType.asym_oneparam ) {
			alignmentContainer = new PairwiseAlignmentContainer_asymBlendOneParam(tree, distances);
		} else if (PairwiseAlignmentsContainer.blendMethod == PairwiseAlignmentsContainer.BlendType.asym_twoparam ) {
			alignmentContainer = new PairwiseAlignmentContainer_asymBlendTwoParam(tree, distances);
		}
			
	}
	

	protected void initialize() {

		//pass pairs to PairwiseAlignmentsContainer, do/get all p/w alignts.
		mods = new ConsistencyModifiers_AllPairs(node, alignmentContainer);

		firstColWithShape_curr = 0;
		lastColWithShape_curr = 0;
		firstColWithShape_next = -1;
		lastColWithShape_next = -1;
		
		
		//	pessimistic alignment.  Just used to establish an upper bound.
		ConsistencyHeuristicAligner al = new ConsistencyHeuristicAligner( this );
		al.setPessimistic(true);
		al.align();
		int[][] result = SequenceConverter.convertPathToIntAlignment(al.getPath(), A, B);
//    	char[][] result3 = SequenceConverter.convertPathToCharAlignment(al.getPath(), SequenceConverter.convertIntsToSeqs(A.seqs), SequenceConverter.convertIntsToSeqs(B.seqs));

		long estCost = al.getEstimatedCost();
		int pessimisticUpperBound = calcConsistencyCost(result); 
		if (estCost < 0) {
			LogWriter.stdErrLogln("Surprise: cost of alignment is negative");
			throw new GenericOpalException("Surprise: cost of alignment is negative");
		}
		if (pessimisticUpperBound > estCost) {
			LogWriter.stdErrLogln("Surprise: actual cost of alignment is worse than pessimistic estimate");
			throw new GenericOpalException("Surprise: actual cost of alignment is worse than pessimistic estimate");
		}
		
		//optimistic alignment.  Need to save the tables for bound pruning.
		//A and B are reversed because the tables are computed from end to beginning
		int[][] revA = SequenceConverter.buildReverseAlignment(A.seqs);
		int[][] revB = SequenceConverter.buildReverseAlignment(B.seqs);
		Alignment alA = Alignment.buildNewAlignment(revA, A.seqIds);
		Alignment alB = Alignment.buildNewAlignment(revB, B.seqIds);
		al = new ConsistencyHeuristicAligner(this);
		al.setPessimistic(false);
		al.setReverseMode(true);
		al.setAlignments(alA, alB);
		al.align();

//		char[][] result2 = SequenceConverter.convertPathToCharAlignment(al.getPath(), SequenceConverter.convertIntsToSeqs(alA.seqs), SequenceConverter.convertIntsToSeqs(alB.seqs));

		result = SequenceConverter.convertPathToIntAlignment (al.getPath(), alA, alB);
		estCost = al.getEstimatedCost();
		int costUpperBound = calcConsistencyCost(result,true);
		if (estCost < 0) {
			LogWriter.stdErrLogln("Surprise: cost of alignment is negative");
			throw new GenericOpalException("Surprise: cost of alignment is negative");
		}
		if (costUpperBound < estCost) {
			LogWriter.stdErrLogln("Surprise: actual cost of alignment is better than optimistic estimate");
			throw new GenericOpalException("Surprise: actual cost of alignment is better than optimistic estimate");
		}

		costUpperBound = Math.min(pessimisticUpperBound, costUpperBound); 
		shapeTester = new ConsistencyShapeTester (costUpperBound, al.getD(), al.getH(), al.getV());
		al.setAlignments(A, B); // back to the non-reverse seqs

		// let the GC clean these up
		alA = alB = null;
		revA = revB = result = null;		
					
		//static methods, so I don't need to keep passing these in.
		Shape.setAlignments (A, B);
		Shape.setParams(costs);
		ConsistencyShape.setMods(mods);
		
		currRow = 0;
		nextRow = 1;
			  
		//allocate dynamic programming table 
		dpRows = new ArrayList [3][N+1]; // 3 rows to cycle fill in the table. Each row of length N+1, each cell a list of shapes
		for (int i=0; i<3; i++){
			for(int j=0; j<= N; j++){
				dpRows[i][j] = new ArrayList<ConsistencyShape>();
			}
		}
		dpRows[currRow][0].add(new ConsistencyShape()); // initialize entry (0,0) with the flush shape and cost 0	  
	
	}

	
	protected void fillTable () {		
		int j;
		
		//fill in the dynamic programming table (row major order) 
		for (int i = 0; i <= M; i++) {
			
			j = firstColWithShape_curr;
			while (j <= lastColWithShape_curr){
				Iterator<ConsistencyShape> shapeIterator = dpRows[currRow][j].iterator();

				while (shapeIterator.hasNext()) {
					
					ConsistencyShape s = shapeIterator.next();
					ConsistencyShape s_new;
					if ( i < M ) {
						//*** propagate from D[i][j] to D[i+1][j] ***//*
						// determine new shape
						s_new = (ConsistencyShape)Shape.makeNewShape(s);
						s_new.appendColumns(i+1, -1);
						if ( !shapeTester.boundPrune(s_new) &&  !shapeTester.dominancePrune(s_new,dpRows[nextRow][j])) {
							dpRows[nextRow][j].add(s_new);
							if (firstColWithShape_next==-1 || firstColWithShape_next>j )
								firstColWithShape_next = j;
							if (lastColWithShape_next<j)
								lastColWithShape_next = j;
						}
					}
					if ( j < N ) {  
						//*** propagate from D[i][j] to D[i][j+1] ***//*
						s_new = (ConsistencyShape)Shape.makeNewShape(s);
						s_new.appendColumns(-1, j+1);
						if ( !shapeTester.boundPrune(s_new) &&  !shapeTester.dominancePrune(s_new,dpRows[currRow][j+1])) {
							dpRows[currRow][j+1].add(s_new);
							if (lastColWithShape_curr == j)
								lastColWithShape_curr = j+1;
						}
					}  	
					if ( i < M && j < N ) {  
						//*** propagate from D[i][j] to D[i+1][j+1] ***//*
						s_new = (ConsistencyShape)Shape.makeNewShape(s);
						s_new.appendColumns(i+1, j+1);
						if ( !shapeTester.boundPrune(s_new) &&  !shapeTester.dominancePrune(s_new,dpRows[nextRow][j+1])) {
							dpRows[nextRow][j+1].add(s_new);
							if (firstColWithShape_next==-1)
								firstColWithShape_next = j+1;
							if (lastColWithShape_next<j+1)
								lastColWithShape_next = j+1;

						}
					}
					if (i==M && j==N) {
						//for each shape in the last cell, 
						//need to add closing costs for all gaps still open at completion 
						s.closeGaps();
					}
					s.freeUnusedStructures();
				}
				j++;

			}			
			incrementRow();

		} 


	}  
	

	protected int calcConsistencyCost(int[][] C) {
		return calcConsistencyCost(C, false);
	}	
	
	protected int calcConsistencyCost(int[][] C, boolean reverseMode) {
		int len = C[0].length;	
		
		int cost = 0;
		int a,b;
		Direction inGap = null;
		for (int p=0; p<K; p++) {
			int mm = A.posInUgappedString[p][M];
			for (int q=0 ; q<L; q++) {
				int nn = B.posInUgappedString[q][N];
				
				int cc = 0, cc2 = 0;
				inGap = null;

				int beg=0;
				//skip the initial and terminal all-space columns
				while (SequenceConverter.GAP_VAL == C[p][beg] && SequenceConverter.GAP_VAL == C[q+K][beg]) beg++;
				int end = len-1;
				while (SequenceConverter.GAP_VAL == C[p][end] && SequenceConverter.GAP_VAL == C[q+K][end]) end--;
				
				ConsistencyModifiers_Pair modpair =  mods.modifiers[p][q];						
				
				int posA=0, posB=0;
				
				for (int i=beg; i<=end; i++) {
					a = C[p][i];
					if (a!=SequenceConverter.GAP_VAL) posA++;
					b = C[q+K][i];					
					if (b!=SequenceConverter.GAP_VAL) posB++;
							
					if (a != SequenceConverter.GAP_VAL && b != SequenceConverter.GAP_VAL) { // substitution\
						if (reverseMode) cc += modpair.subs[mm-(posA-1)][nn-(posB-1)]; 
						else cc += modpair.subs[posA][posB];
						
						if (inGap == Direction.vert) {
							if (reverseMode) cc += modpair.vGammaOpens[mm-(posA-1-1)][nn-(posB-1)];
							else cc += modpair.vGammaCloses[posA-1][posB-1];
						} else if (inGap == Direction.horiz) {
							if (reverseMode) cc += modpair.hGammaOpens[mm-(posA-1)][nn-(posB-1-1)];
							else cc += modpair.hGammaCloses[posA-1][posB-1];
						}
						inGap = null;
					} else if (a != SequenceConverter.GAP_VAL) { //gap in B, i.e. VERT
						if (reverseMode) cc += modpair.vLambdas[mm-(posA-1)][nn-posB];
						else  cc += modpair.vLambdas[posA][posB];
						
						if (inGap != Direction.vert) { 
							if (reverseMode) cc += modpair.vGammaCloses[mm-(posA-1)][nn-posB];
							else  cc += modpair.vGammaOpens[posA][posB];
							
							if (inGap == Direction.horiz)  
								if (reverseMode) cc += modpair.hGammaOpens[mm-(posA-1)][nn-(posB-1)];
								else cc += modpair.hGammaCloses[posA-1][posB];
							
							inGap = Direction.vert;	
						}
					} else if (b != SequenceConverter.GAP_VAL) { //gap in A, i.e. HORIZ
						if (reverseMode) cc += modpair.hLambdas[mm-posA][nn-(posB-1)];
						else  cc += modpair.hLambdas[posA][posB];
						
						if (inGap != Direction.horiz) { 
							if (reverseMode) cc += modpair.hGammaCloses[mm-posA][nn-(posB-1)];
							else  cc += modpair.hGammaOpens[posA][posB];
							
							if (inGap == Direction.vert)  
								if (reverseMode) cc += modpair.vGammaOpens[mm-(posA-1)][nn-(posB-1)];
								else  cc += modpair.vGammaCloses[posA][posB-1];
							
							inGap = Direction.horiz;	
						}
					} // otherwise, it's a double-dash column 

					
					if (testing) LogWriter.stdOutLogln("col " + i + " : " + cc);
					cc2 += cc;
					cc = 0;
				}

				//all done; did I just close a gap?
				int xx = 0;
				if (inGap == Direction.vert) {
					if (reverseMode) xx = modpair.vGammaOpens[1][0];
					else xx = modpair.vGammaCloses[mm][nn];
				} else if (inGap == Direction.horiz) {
					if (reverseMode) cc2 += modpair.hGammaOpens[0][1];
					else cc2 += modpair.hGammaCloses[mm][nn];
				}
				if (xx>0) {
					if (testing) LogWriter.stdOutLogln("a little more: " + xx);
					cc2 += xx;
				}

				cost += cc2;
			}
		}	
		
		return cost;
	}	

	

	public PairwiseAlignmentsContainer getAlignmentContainer () {
		return alignmentContainer;
	}

	public ConsistencyModifiers_AllPairs getConsistencyModifiers () {
		return mods;
	}
	
}

package opal.align;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import opal.align.shapes.*;
import opal.exceptions.GenericOpalException;
import com.traviswheeler.libs.LogWriter;
import opal.IO.SequenceConverter;


public class ExactCountAligner_Time extends ExactCountAligner {
	
	/* I don't parameterize this ArrayList because it's a "bad"  
	   thing to have an array of collections (should be a collection of collections)
	   ... but I need that for the speed to be acceptable */
	ArrayList[][] dpRows;
	
	
	int prevRow = 0;
	int currRow = 0;
	int nextRow = 1;
	long costUpperBound;	
	ShapeTester shapeTester;

	int firstColWithShape_curr;
	int lastColWithShape_curr;
	int firstColWithShape_next;
	int lastColWithShape_next;
	
	
	public ExactCountAligner_Time() {
		this(true);
	}
	public ExactCountAligner_Time( Alignment A, Alignment B ) {
		this(A,B,true);
	}
	
	public ExactCountAligner_Time( boolean pess ) {
		super(pess);		
	}
	
	public ExactCountAligner_Time( Alignment A, Alignment B, boolean pess ) {
		super(A,B,pess);		
	}
	
	protected void fillBoundary (){
	}
	
	protected void initialize () {		

		firstColWithShape_curr = 0;
		lastColWithShape_curr = 0;
		firstColWithShape_next = -1;
		lastColWithShape_next = -1;

		//optimistic alignment.  Need to save the tables for bound pruning.
		//A and B are reversed because the tables are computed from end to beginning
		int[][] revA = config.sc.buildReverseAlignment(A.seqs);
		int[][] revB = config.sc.buildReverseAlignment(B.seqs);
		
		Alignment alA = Alignment.buildNewAlignment(revA, A.seqIds, true, config);
		Alignment alB = Alignment.buildNewAlignment(revB, B.seqIds, true, config);
		
		
		ProfileAligner al = new ProfileAligner( alA, alB);
		al.config = config;
		al.setPessimistic(false);
		al.align();
		
		int[] idsAB = new int[A.seqIds.length + B.seqIds.length];
		System.arraycopy(A.seqIds, 0, idsAB, 0, A.seqIds.length);
		System.arraycopy(B.seqIds, 0, idsAB, A.seqIds.length, B.seqIds.length);

		
		int[][] alignedSeqs = null;
		if (config.useStructure) { // structure is stored in forward orientation. This is easiest way to get correct score calculated
			alignedSeqs = config.sc.buildReverseAlignment(al.getAlignment().seqs);
		} else {
			alignedSeqs = al.getAlignment().seqs;
		}
				
		costUpperBound = Aligner.calcCost(alignedSeqs,A.K,B.K,idsAB, config);

		if (costUpperBound < al.getEstimatedCost()) {
			LogWriter.stdErrLogln("Surprise: actual cost of alignment is better than optimistic estimate");
			throw new GenericOpalException("Surprise: actual cost of alignment is better than optimistic estimate");
		}
		
		if (A.seqs.length * B.seqs.length > Aligner.linearCutoff)
			shapeTester = new ShapeTesterLinear (costUpperBound, al.getD(), al.getH(), al.getV(), config);
		else
			shapeTester = new ShapeTesterQuadratic (costUpperBound, al.getD(), al.getH(), al.getV(), config);
		

		// let the GC clean these up
		alignedSeqs = null;
		alA = alB = null;
		revA = revB = null;		
			
		//	pessimistic alignment.  Just used to establish an upper bound.
		al = new ProfileAligner(this);
		al.setPessimistic(true);
		al.align();
		
		long pessimisticUpperBound = Aligner.calcCost(al.getAlignment().seqs,A.K,B.K, idsAB, config);

		if (pessimisticUpperBound > al.getEstimatedCost()) {
			LogWriter.stdErrLogln("Surprise: actual cost of alignment is worse than pessimistic estimate");
			throw new GenericOpalException("Surprise: actual cost of alignment is worse than pessimistic estimate");
		}
		
		if (pessimisticUpperBound < costUpperBound) 
			costUpperBound = pessimisticUpperBound;
		
		//static methods, so I don't need to keep passing these in.
		//Shape.setAlignments (A, B);
		//Shape.setParams(config.cost.costs);
		
		currRow = 0;
		nextRow = 1;
			  
		//allocate dynamic programming table 
		dpRows = new ArrayList [3][N+1]; // 3 rows to cycle fill in the table. Each row of length N+1, each cell a list of shapes
		for (int i=0; i<3; i++){
			for(int j=0; j<= N; j++){
				dpRows[i][j] = new ArrayList<Shape>();
			}
		}
		dpRows[currRow][0].add(new ShapeQuadratic(config, A, B)); // initialize entry (0,0) with the flush shape and cost 0	  
	}	
	
	
	public void cleanup () {
		dpRows = null;
	}

	
	protected void incrementRow() {
		prevRow = currRow;
		currRow = nextRow;
		nextRow = (nextRow+1)%3; // if I'm on index 2, recycle 0.
		for (int i= 0; i<dpRows[nextRow].length; i++) {
			dpRows[nextRow][i].clear(); // we'll be reusing this row, so remove existing shapes,
		}

		firstColWithShape_curr = firstColWithShape_next;
		lastColWithShape_curr = lastColWithShape_next;
		firstColWithShape_next = -1;
		lastColWithShape_next = -1;
	}

	
	protected void fillTable () {
	
		int j;
		//fill in the dynamic programming table (row major order) 
		for (int i = 0; i <= M; i++) {
			//for (int j = 0; j <= N; j++) {
			j = firstColWithShape_curr;

			while (j <= lastColWithShape_curr){

				if(j<0) System.err.println("J: " + j + "\trow:" + currRow + "\t" + config);
				Iterator shapeIterator = dpRows[currRow][j].iterator();
				
				//Logger.stdOutLogln(i + "," + j + " : " + dpRows[currRow][j].size() + " shapes");
				
				while (shapeIterator.hasNext()) {
					Shape s = (Shape)(shapeIterator.next());
					Shape s_new;
					if ( i < M ) {
						//*** propagate from D[i][j] to D[i+1][j] ***//*
						// determine new shape
						s_new = Shape.makeNewShape(s);
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
						s_new = Shape.makeNewShape(s);
						s_new.appendColumns(-1, j+1);
						if ( !shapeTester.boundPrune(s_new) &&  !shapeTester.dominancePrune(s_new,dpRows[currRow][j+1])) {
							dpRows[currRow][j+1].add(s_new);
							if (lastColWithShape_curr == j)
								lastColWithShape_curr = j+1;

						}
					}  	
					if ( i < M && j < N ) {  
						//*** propagate from D[i][j] to D[i+1][j+1] ***//*
						s_new = Shape.makeNewShape(s);
						s_new.appendColumns(i+1, j+1);
						if ( !shapeTester.boundPrune(s_new) &&  !shapeTester.dominancePrune(s_new,dpRows[nextRow][j+1])) {
							dpRows[nextRow][j+1].add(s_new);
							if (firstColWithShape_next==-1)
								firstColWithShape_next = j+1;
							if (lastColWithShape_next<j+1)
								lastColWithShape_next = j+1;

						}
					}
					s.freeUnusedStructures();
				}
				j++;

			}			
			incrementRow();

		} 
		
		//System.err.println("Finish fill " + config);
			
	}  

	
	protected boolean recover () {
		path = new ArrayList<Direction>(2*Math.max(M,N));

		//find the optimal solution shape in D[m][n]

		Iterator shapeIterator = dpRows[prevRow][N].iterator();
		Shape s = (Shape)(shapeIterator.next());
		Shape bestShape = s;
		while (shapeIterator.hasNext()) {
			s = (Shape)(shapeIterator.next());
			if (s.getCost() < bestShape.getCost() ){
				bestShape = s;
			}
		}

		s = bestShape;
//		s.printInfo();
		while (s.getParent() != null) {
			Direction i = s.getParentDirection();
			if (i==null) return false;
			path.add(i);
			s = s.getParent();
//			s.printInfo();
		}
 
		dpRows = null; // let GC clean up shapes
		Collections.reverse(path);
		
		return true;
	}
	
}

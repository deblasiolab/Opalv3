package opal.align.shapes;

import opal.IO.SequenceConverter;
import opal.align.Aligner;
import opal.align.Alignment;
import opal.IO.Configuration;

public class ShapeLinear extends Shape {

	public ShapeLinear(Configuration c, Alignment A, Alignment B) {
		super(c,A,B);
	}

	public ShapeLinear(Shape s) {
		super(s);
	}
	
	 /* procedure: gapOpenCost
	 *   Input:
	 *     a, b: column indices into A and B respectively
	 *        s: integer array representation of a shape
	 *   Output: gc, a int representing the sum of pairwise
	 *           gap startup penalties between A and B (not within A, nor B)
	 *           resulting from appending column a of A and column b of B 
	 *           to shape s (a < 0 indicates deletion, b < 0 indicates insertion).
	 *           
	 *  This is a O(K + L) version, based on the ideas of O. Gotoh, "Further 
	 *  improvements in methods of group-togroup sequence alignment with 
	 *  generalized profile operations. CABIOS 10(4):379-387, 1994.
	 */
	final protected long gapBoundaryCost (int a, int b, Aligner.Direction dir) {
		int i = 0;
		int[] Ta = new int[maxBlock+1]; //Ta[i] := num gap chars in block i of alignment A
		int[] Tb = new int[maxBlock+1]; //Tb[i] := num gap chars in block i of alignment B
		int[] Ta_left_term = new int[maxBlock+1]; //Ta_term[i] := num gap chars in block i of alignment A that are outside the boundaries of the string for their seqs
		int[] Tb_left_term = new int[maxBlock+1]; //Tb_term[i] := num gap chars in block i of alignment B that are outside the boundaries of the string for their seqs

		int[] Sa = new int[maxBlock+1]; //Sa[i] := (1) num non-gap chars in block i of alignment A
		 								//   then (2) num strings whose final non-gap char occurs at or before block i in A.
		int[] Sb = new int[maxBlock+1];//Sb[i] := same Sa, but for B.

		int[] Ta_right_term = new int[maxBlock+1]; //Ta_term[i] := num gap chars in block i of alignment A that are outside the boundaries of the string for their seqs
		int[] Tb_right_term = new int[maxBlock+1]; //Tb_term[i] := num gap chars in block i of alignment B that are outside the boundaries of the string for their seqs


		
		for (i = 0; i <= maxBlock; i++) {
		    Ta[i] = Tb[i] = Sa[i] = Sb[i] = 0;
		}
		
		if (Aligner.Direction.horiz != dir)  { // VERT or DIAG  ... uses col A
			for (i = 0; i < K; i++) {
				if (a == 0 || SequenceConverter.GAP_VAL == A.seqs[i][a-1]) { 
					Ta[seqBlocks[i]]++;
					if (a < A.firstLetterLoc[i]) 
						Ta_left_term[seqBlocks[i]]++;
					if (a >= A.lastLetterLoc[i]) 
						Ta_right_term[seqBlocks[i]]++;						
				} else { 
					Sa[seqBlocks[i]]++; //(1) 
				}  
		    }
			for (i = 1; i <= maxBlock; i++) { 
				Sa[i] += Sa[i-1]; //(2)
		    }

		}
		if (Aligner.Direction.vert != dir)  { // HORIZ or DIAG  ... uses col B
		    for (i = 0; i < L; i++) {
			      if (b==0 || SequenceConverter.GAP_VAL == B.seqs[i][b-1]) { 
			    	  Tb[seqBlocks[i+K]]++; 
			    	  if (b < B.firstLetterLoc[i]) 
							Tb_left_term[seqBlocks[i+K]]++;	
			    	  if (b >= B.lastLetterLoc[i]) 
							Tb_right_term[seqBlocks[i+K]]++;	
			      } else { 
			    	  Sb[seqBlocks[i+K]]++; 
		    	  }
			}	
		    for (i = 1; i <= maxBlock; i++) { 
			      Sb[i] += Sb[i-1];
			}
		}
		
		if (Aligner.Direction.vert == dir) { 		
			for (i = 0; i < L; i++) { 
				Tb[seqBlocks[i+K]]++;
				if (bPos < B.firstLetterLoc[i]) 
					Tb_left_term[seqBlocks[i+K]]++;
				if (bPos >= B.lastLetterLoc[i]) 
						Tb_right_term[seqBlocks[i+K]]++;
			}
		} else if (Aligner.Direction.horiz == dir) { 
			for (i = 0; i < K; i++) { 
		       Ta[seqBlocks[i]]++;
		       if (aPos < A.firstLetterLoc[i]) 
					Ta_left_term[seqBlocks[i]]++;	
		       if (aPos >= A.lastLetterLoc[i]) 
					Ta_right_term[seqBlocks[i]]++;	
			}
		}
		
		long gammaCount = 0;
		long termGammaCountLeft = 0;
		long termGammaCountRight = 0;
		for (i = 0; i <= maxBlock; i++) {
			gammaCount += ( (Ta[i]-Ta_left_term[i]-Ta_right_term[i]) * Sb[i] + (Tb[i]-Tb_left_term[i]-Tb_right_term[i]) * Sa[i]);
			termGammaCountLeft += ( Ta_left_term[i] * Sb[i] + Tb_left_term[i] * Sa[i]);
			termGammaCountRight += ( Ta_right_term[i] * Sb[i] + Tb_right_term[i] * Sa[i]);
		}
		return gammaCount * config.gamma + termGammaCountLeft * config.leftGammaTerm() + termGammaCountRight * config.rightGammaTerm();
	}


}

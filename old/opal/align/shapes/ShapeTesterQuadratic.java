package opal.align.shapes;

import opal.IO.SequenceConverter;
import opal.IO.StructureFileReader;
import opal.align.Aligner;
import opal.align.Alignment;
import opal.align.StructureAlignment;

public class ShapeTesterQuadratic extends ShapeTester {

	public ShapeTesterQuadratic(long upperBound, long[][] lowerD, long[][] lowerH, long[][] lowerV) {
		super(upperBound, lowerD, lowerH, lowerV);
		// TODO Auto-generated constructor stub
	}

	
	final protected void calcGapBounds(Shape s, long[] H_V_D) {
		/* Hgap, Vgap, and Dgap are upperbounds on the cost (unweighted? aka cost = number) of gaps
		 * that could be open in both the alignment so far and the
		 * computation of the extension (which was done via optimistic
		 * on the reverse of the input )
		 */
	
		Alignment A = Shape.A;
		Alignment B = Shape.B;
		int K = Shape.K;
		int L = Shape.L;
		int M = Shape.M;
		int N = Shape.N;
		int locGamma;
		int Hgap, Vgap, Dgap;
		Hgap = Vgap = Dgap = 0;  
		for (int i = 0; i < K; i++) {
			for (int j = 0; j < L; j++) { // for each pair of seqs
				if (s.seqBlocks[i] < s.seqBlocks[K+j]) {  //underhanging 
					//an extension starting with deletion could count a gap twice 
					locGamma = Aligner.gamma;
					if (s.aPos < A.firstLetterLoc[i] || s.aPos > A.lastLetterLoc[i]) { 
						locGamma = Aligner.gammaTerm;
					} else if (Aligner.useStructure && s.aPos > 0) {
							int pos = ((StructureAlignment)A).origSeqIndices[i][s.aPos-1];
							if (pos>0)
								locGamma += StructureAlignment.gapOpenMods[ StructureFileReader.structureNeighborLevels[A.seqIds[i]][pos] ] ;
					}
					Hgap +=  locGamma; //* weights[i][aa_K+j]
					if (s.aPos != M && A.seqs[i][s.aPos] == SequenceConverter.GAP_VAL) {
					   //cannot determine this case, so add to the upperbound 
					   Vgap += locGamma; // * weights[i][aa_K+j]
					   Dgap += locGamma;
				   }
			   } else if (s.seqBlocks[i] > s.seqBlocks[K+j]) { // overhang
				   //an extension starting with insertion could count a gap twice 
					locGamma = Aligner.gamma;
					if (s.bPos < B.firstLetterLoc[j] || s.bPos > B.lastLetterLoc[j]) {
						locGamma = Aligner.gammaTerm;
					} else if (Aligner.useStructure && s.bPos > 0) {
						int pos = ((StructureAlignment)B).origSeqIndices[j][s.bPos-1];
						if (pos>0)
							locGamma += StructureAlignment.gapOpenMods[ StructureFileReader.structureNeighborLevels[B.seqIds[j]][pos] ] ;
					}
					Vgap += locGamma;
					if (s.bPos != N && B.seqs[j][s.bPos] == SequenceConverter.GAP_VAL) {
					   //cannot determine this case, so add to the upperbound 
					   Hgap += locGamma;
					   Dgap += locGamma;
					}
			   }
		   }
		}
		H_V_D[0] = Hgap;
		H_V_D[1] = Vgap;
		H_V_D[2] = Dgap;
	}
	
	/* procedure: domGapBound
	 *   Input:  shapes s and t
	 *   Output: cost, a int representing an upperbound on the cost of
	 *           gaps that could possibly be started in s but not in t,
	 *           over all possible extensions.
	 *  This is a O(KL) version 
	 */
	final protected long domGapBound (Shape s, Shape t) {
		int K = Shape.K;
		int L = Shape.L;
		Alignment A = Shape.A;
		Alignment B = Shape.B;

		long cost = 0;
		int i,j;
		int locGamma;
		for ( i=0; i<K; i++) {
			for ( j=K; j<K+L; j++) {
				locGamma = Aligner.gamma;
				if 	( t.seqBlocks[i]>t.seqBlocks[j] && s.seqBlocks[i]<=s.seqBlocks[j] ) {
					if (s.bPos < B.firstLetterLoc[j-K] || s.bPos > B.lastLetterLoc[j-K]) {
						locGamma = Aligner.gammaTerm;
					} else if (Aligner.useStructure && s.bPos > 0) {
						int pos = ((StructureAlignment)B).origSeqIndices[j-K][s.bPos-1];
						locGamma += StructureAlignment.gapOpenMods[ StructureFileReader.structureNeighborLevels[B.seqIds[j-K]][pos] ] ;
					}
					cost += locGamma;
				} else if ( t.seqBlocks[i]<t.seqBlocks[j] && s.seqBlocks[i]>=s.seqBlocks[j] ) {
					if (s.aPos < A.firstLetterLoc[i] || s.aPos > A.lastLetterLoc[i]) {
						locGamma = Aligner.gammaTerm;
					} else if (Aligner.useStructure && s.aPos > 0) {
						int pos = ((StructureAlignment)A).origSeqIndices[i][s.aPos-1];
						locGamma += StructureAlignment.gapOpenMods[ StructureFileReader.structureNeighborLevels[A.seqIds[i]][pos] ] ;
					}
					cost += locGamma;		
				}
			}
		}
		return cost;
	}


}

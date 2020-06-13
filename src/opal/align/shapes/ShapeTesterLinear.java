package opal.align.shapes;

import opal.IO.SequenceConverter;
import opal.align.Aligner;
import opal.align.Alignment;
import opal.IO.Configuration;

public class ShapeTesterLinear extends ShapeTester {

	public ShapeTesterLinear(long upperBound, long[][] lowerD, long[][] lowerH, long[][] lowerV, Configuration c) {
		super(upperBound, lowerD, lowerH, lowerV, c);
		// TODO Auto-generated constructor stub
	}

	
	final protected void calcGapBounds(Shape s, long[] H_V_D) {
		/* Hgap, Vgap, and Dgap are upperbounds on the cost (unweighted? aka cost = number) of gaps
		 * that could be open in both the alignment so far and the
		 * computation of the extension (which was done via optimistic
		 * on the reverse of the input )
		 */	
		int K = s.K;
		int L = s.L;

		int i;
		int sBlocknum;  
		
		//precomputation 		
		sBlocknum = 0;  
		for (i = 0; i < K+L; i++) {
			if (s.seqBlocks[i] > sBlocknum) { sBlocknum = s.seqBlocks[i]; }
		}
		sBlocknum++;
	  
		int[] Bas = new int[sBlocknum];  // # seqs in block i from A in shape s
		int[] Bbs = new int[sBlocknum];  // # seqs in block i from B in shape s
		int[] Basc = new int[sBlocknum]; // cumulative: # non-terminal seqs in block <=i from A in shape s
		
		int[] Basc_term_left = new int[sBlocknum]; // cumulative: # terminal seqs in block <=i from A in shape s
		int[] Basc_term_right = new int[sBlocknum]; // cumulative: # terminal seqs in block <=i from A in shape s
		
		int[] Bbsc = new int[sBlocknum]; // cumulative: # non-terminal seqs in block <=i from B in shape s
		
		int[] Bbsc_term_left = new int[sBlocknum]; // cumulative: # terminal seqs in block <=i from B in shape s
		int[] Bbsc_term_right = new int[sBlocknum]; // cumulative: # terminal seqs in block <=i from B in shape s
		
		int[] Basc0 = new int[sBlocknum]; //cumulative: # non-terminal seqs in block <=i from A in shape s, such that the next column of A is a "-"  
		int[] Bbsc0 = new int[sBlocknum];  //cumulative: # non-terminal seqs in block <=i from A in shape s, such that the next column of A is a "-"

		int[] Basc0_term_left = new int[sBlocknum]; //cumulative: # terminal seqs in block <=i from A in shape s, such that the next column of A is a "-"  
		int[] Basc0_term_right = new int[sBlocknum]; //cumulative: # terminal seqs in block <=i from A in shape s, such that the next column of A is a "-"  
		int[] Bbsc0_term_left = new int[sBlocknum];  //cumulative: # terminal seqs in block <=i from A in shape s, such that the next column of A is a "-"
		int[] Bbsc0_term_right = new int[sBlocknum];  //cumulative: # terminal seqs in block <=i from A in shape s, such that the next column of A is a "-"

		buildBlockArrays (s, Bas, Bbs, Basc, Basc_term_left, Basc_term_right, Bbsc, Bbsc_term_left, Bbsc_term_right, 
				Basc0, Basc0_term_left,Basc0_term_right, Bbsc0, Bbsc0_term_left,Bbsc0_term_right);

		
		int Hgap, Vgap, Dgap;
		Hgap = Vgap = Dgap = 0;  
		long val;

		// for each block, check how many seqs in A overhang blocks in B, then vice versa
		for (int block=1; block<sBlocknum; block++){
			//overhangs
			val = Bas[block] * Bbsc[block-1] * config.gamma  
					+ Bas[block] * Bbsc_term_left[block-1] * config.leftGammaTerm()
					+ Bas[block] * Bbsc_term_right[block-1] * config.rightGammaTerm(); 
			Vgap += val;
			val = Bas[block] * Bbsc0[block-1] * config.gamma  
					+ Bas[block] * Bbsc0_term_left[block-1] * config.leftGammaTerm()
					+ Bas[block] * Bbsc0_term_right[block-1] * config.rightGammaTerm();
			Hgap += val;
			Dgap += val;
			
			//underhangs
			val = block==0 ? 0 : 
				Bbs[block] * Basc[block-1] * config.gamma  
				+ Bbs[block] * Basc_term_left[block-1] * config.leftGammaTerm()
				+ Bbs[block] * Basc_term_right[block-1] * config.rightGammaTerm() ; 
			Hgap += val;
			val = block==0 ? 0 : 
				Bbs[block] * Basc0[block-1] * config.gamma  
				+ Bbs[block] * Basc0_term_left[block-1] * config.leftGammaTerm()
				+ Bbs[block] * Basc0_term_right[block-1] * config.rightGammaTerm() ;
			Vgap += val;
			Dgap += val;
		}
		
		H_V_D[0] = Hgap;
		H_V_D[1] = Vgap;
		H_V_D[2] = Dgap;
		
	}
	
	
	
	/* procedure: domGapBoundLinear
	 *   Input:  shapes s and t
	 *   Output: gc, a int representing an upperbound on the cost of
	 *           gaps that could possibly be started in s but not in t,
	 *           over all possible extensions.
	 *  This is a O(K + L) version, based loosely on the ideas in: 
	 *  O. Gotoh, "Further improvements in methods of group-to-group 
	 *  sequence alignment with generalized profile operations. 
	 *  CABIOS 10(4):379-387, 1994."   
	 */
	final protected long domGapBound (Shape s, Shape t) {

		int K = s.K;
		int L = s.L;

		int i;
		int sBlocknum, tBlocknum;  
		
		//precomputation 		
		sBlocknum = tBlocknum = 0;  
		for (i = 0; i < K+L; i++) {
			if (s.seqBlocks[i] > sBlocknum) { sBlocknum = s.seqBlocks[i]; }
			if (t.seqBlocks[i] > tBlocknum) { tBlocknum = t.seqBlocks[i]; }
		}
		sBlocknum++;
		tBlocknum++;  
	  
		int[] Bas = new int[sBlocknum];  // # seqs in block i from A in shape s
		int[] Bbs = new int[sBlocknum];  // # seqs in block i from B in shape s
		int[] Basc = new int[sBlocknum]; // cumulative: # non-terminal seqs in block <=i from A in shape s
		
		int[] Basc_term_left = new int[sBlocknum]; // cumulative: # terminal seqs in block <=i from A in shape s
		int[] Basc_term_right = new int[sBlocknum]; // cumulative: # terminal seqs in block <=i from A in shape s
		
		int[] Bbsc = new int[sBlocknum]; // cumulative: # non-terminal seqs in block <=i from B in shape s
		
		int[] Bbsc_term_left = new int[sBlocknum]; // cumulative: # terminal seqs in block <=i from B in shape s
		int[] Bbsc_term_right = new int[sBlocknum]; // cumulative: # terminal seqs in block <=i from B in shape s
		
		buildBlockArrays (s, Bas, Bbs, Basc, Basc_term_left, Basc_term_right, Bbsc, Bbsc_term_left, Bbsc_term_right);

		
				
		int[] Bat = new int[tBlocknum];  // # seqs in block i from A in shape t
		int[] Bbt = new int[tBlocknum];  // # seqs in block i from B in shape t		
		int[] Batc = new int[tBlocknum]; // cumulative: # non-terminal seqs in block <=i from A in shape t
		
		int[] Batc_term_left = new int[tBlocknum]; // cumulative: # terminal seqs in block <=i from A in shape t
		int[] Batc_term_right = new int[tBlocknum]; // cumulative: # terminal seqs in block <=i from A in shape t
		
		int[] Bbtc = new int[tBlocknum]; // cumulative: # non-terminal seqs in block <=i from B in shape t
		
		int[] Bbtc_term_left = new int[tBlocknum]; // cumulative: # terminal seqs in block <=i from B in shape t
		int[] Bbtc_term_right = new int[tBlocknum]; // cumulative: # terminal seqs in block <=i from B in shape t
		
		buildBlockArrays (t, Bat, Bbt, Batc, Batc_term_left, Batc_term_right, Bbtc, Bbtc_term_left, Bbtc_term_right);
		
		int ot; // number non-terminal of seqs in B that A[i] overhangs (in shape t)
		int ot_term_left; // number terminal of seqs in B that A[i] overhangs (in shape t)
		int ot_term_right; // number terminal of seqs in B that A[i] overhangs (in shape t)			
		int os; // number of non-terminal seqs in B that A[i] overhangs (in shape s)
		int os_term_left; // number of terminal seqs in B that A[i] overhangs (in shape s)
		int os_term_right; // number of terminal seqs in B that A[i] overhangs (in shape s)
		int sBlock, tBlock;
		long cost = 0;
		// For each seq i in A, count number of seqs j in B s.t. i overhangs j in s, but doesn't in t.
		// This gives an upper bound on the number of new gap opens that some extension could
		// cause in s without causing them in t.  
		for (i=0; i<K; i++) {
			tBlock = t.seqBlocks[i];
			if (tBlock == 0 )
				continue;	
			ot = Bbtc[tBlock-1];
			ot_term_left = Bbtc_term_left[tBlock-1];
			ot_term_right = Bbtc_term_right[tBlock-1];

			sBlock = s.seqBlocks[i];
			if (sBlock>0) {
				os = Bbsc[sBlock-1];
				os_term_left = Bbsc_term_left[sBlock-1];
				os_term_right = Bbsc_term_right[sBlock-1];
			} else {
				os = os_term_left = os_term_right = 0; 
			}
					
			int ubNewGaps = ot - os ;// upper bound on non-terminal gap opens hitting s but not t
			if (ubNewGaps<0) ubNewGaps = 0;
			int ubNewGaps_term_left = ot_term_left - os_term_left ;// upper bound on terminal gap opens hitting s but not t
			int ubNewGaps_term_right = ot_term_right - os_term_right ;// upper bound on terminal gap opens hitting s but not t
			if (ubNewGaps_term_left<0) ubNewGaps_term_left = 0;
			if (ubNewGaps_term_right<0) ubNewGaps_term_right = 0;
			
 			cost += (config.gamma * ubNewGaps) + (config.leftGammaTerm() * ubNewGaps_term_left) + (config.rightGammaTerm() * ubNewGaps_term_right); 
		}
		
		// For each seq i in B, count number of seqs j in A s.t. i underhangs j in t, but doesn't in s.
		for (i=0; i<L; i++) {
			tBlock = t.seqBlocks[K+i];
			if (tBlock == 0 )
				continue;	
			ot = Batc[tBlock-1];
			ot_term_left = Batc_term_left[tBlock-1];
			ot_term_right = Batc_term_right[tBlock-1];

			sBlock = s.seqBlocks[K+i];
			if (sBlock>0) {
				os = Basc[sBlock-1];
				os_term_left = Basc_term_left[sBlock-1];
				os_term_right = Basc_term_right[sBlock-1];
			} else {
				os = os_term_left = os_term_right = 0; 
			}		
			
			int ubNewGaps = ot - os ;// upper bound on non-terminal gap opens hitting s but not t
			if (ubNewGaps<0) ubNewGaps = 0;
			int ubNewGaps_term_left = ot_term_left - os_term_left ;// upper bound on terminal gap opens hitting s but not t
			int ubNewGaps_term_right = ot_term_right - os_term_right ;// upper bound on terminal gap opens hitting s but not t
			if (ubNewGaps_term_left<0) ubNewGaps_term_left = 0;
			if (ubNewGaps_term_right<0) ubNewGaps_term_right = 0;
			
 			cost += (config.gamma * ubNewGaps) + (config.leftGammaTerm() * ubNewGaps_term_left) + (config.rightGammaTerm() * ubNewGaps_term_right); 
		}
				
		return  cost;
	}


	/* 
	 * Build block-count arr.ays used by both domGapBound and calcGapBounds
	 */ 
	private void buildBlockArrays (Shape s, int[] Ba, int[] Bb,   
									int[] Bac, int[] Bac_term_left, int[] Bac_term_right, 
									int[] Bbc, int[] Bbc_term_left, int[] Bbc_term_right)  {
		buildBlockArrays (s, Ba, Bb, Bac, Bac_term_left, Bac_term_right, Bbc, Bbc_term_left, Bbc_term_right, null, null, null, null, null, null);
	}

	private void buildBlockArrays (Shape s, int[] Ba, int[] Bb,   
					int[] Bac, int[] Bac_term_left, int[] Bac_term_right, 
					int[] Bbc, int[] Bbc_term_left, int[] Bbc_term_right,
					int[] Bac0, int[] Bac0_term_left, int[] Bac0_term_right, 
					int[] Bbc0, int[] Bbc0_term_left, int[] Bbc0_term_right)  {
		
		int K = s.K;
		int L = s.L;
		int M = s.M;
		int N = s.N;
		Alignment A = s.A;
		Alignment B = s.B;
		
		int i=0;
		int blocknum = Ba.length;  
	
		int[] Ba_term_left = new int[blocknum];  // # terminal seqs in block i from A in shape s  (where "terminal" means the block either preceeds or follows the sequence)
		int[] Ba_term_right = new int[blocknum];  // # terminal seqs in block i from A in shape s  (where "terminal" means the block either preceeds or follows the sequence)
		int[] Bb_term_left = new int[blocknum];  // # terminal seqs in block i from B in shape s
		int[] Bb_term_right = new int[blocknum];  // # terminal seqs in block i from B in shape s
		int[] Ba0_term_left = new int[blocknum];  // # terminal seqs in block i from A in shape s, such that the next column of A is a "-"
		int[] Ba0_term_right = new int[blocknum];  // # terminal seqs in block i from A in shape s, such that the next column of A is a "-"
		int[] Bb0_term_left = new int[blocknum];  // # terminal seqs in block i from B in shape s, such that the next column of B is a "-"
		int[] Bb0_term_right = new int[blocknum];  // # terminal seqs in block i from B in shape s, such that the next column of B is a "-"
		int[] Ba0 = new int[blocknum];  // # non-terminal seqs in block i from A in shape s, such that the next column of A is a "-"
		int[] Bb0 = new int[blocknum];  // # non-terminal seqs in block i from B in shape s, such that the next column of B is a "-"
		
		
		//precomputation 		
		for (i = 0; i < blocknum; i++) {
			Ba[i] = 0;
			Bb[i] = 0;
			Ba_term_left[i] = 0;
			Ba_term_right[i] = 0;
			Bb_term_left[i] = 0;
			Bb_term_right[i] = 0;
			Ba0[i] = 0;
			Bb0[i] = 0;
			Ba0_term_left[i] = 0;
			Ba0_term_right[i] = 0;
			Bb0_term_left[i] = 0;
			Bb0_term_right[i] = 0;
		}
		
		//number of seqs in the A-sequence-set of each block for shape s 
		for (i = 0; i < K; i++) {
			boolean nextIsDash = s.aPos != M && A.seqs[i][s.aPos] == SequenceConverter.GAP_VAL;
			Ba[s.seqBlocks[i]]++;
			if (nextIsDash)	Ba0[s.seqBlocks[i]]++;
			if (s.aPos < A.firstLetterLoc[i]) { 
				Ba_term_left[s.seqBlocks[i]]++;
				if (nextIsDash) Ba0_term_left[s.seqBlocks[i]]++;
			}
			if (s.aPos > A.lastLetterLoc[i]) { 
				Ba_term_right[s.seqBlocks[i]]++;
				if (nextIsDash) Ba0_term_right[s.seqBlocks[i]]++;
			}
		}
		
		//number of seqs in the B-sequence-set of each block for shape s
		for (i = K; i < K + L; i++) {
			boolean nextIsDash = s.bPos != N && B.seqs[i-K][s.bPos] == SequenceConverter.GAP_VAL;
			Bb[s.seqBlocks[i]]++;
			if (nextIsDash)	Bb0[s.seqBlocks[i]]++;
			if (s.bPos < B.firstLetterLoc[i-K]) { 
				Bb_term_left[s.seqBlocks[i]]++;
				if (nextIsDash) Bb0_term_left[s.seqBlocks[i]]++;
			}
			if (s.bPos > B.lastLetterLoc[i-K]) { 
				Bb_term_right[s.seqBlocks[i]]++;
				if (nextIsDash) Bb0_term_right[s.seqBlocks[i]]++;
			}			
		}

		Bac[0] = Ba[0]-Ba_term_left[0]-Ba_term_right[0];
		if (Bac0 != null) Bac0[0] = Ba0[0]-Ba0_term_left[0]-Ba0_term_right[0];

		Bac_term_left[0] = Ba_term_left[0];
		if (Bac0_term_left != null) Bac0_term_left[0] = Ba0_term_left[0];
		Bac_term_right[0] = Ba_term_right[0];
		if (Bac0_term_right != null) Bac0_term_right[0] = Ba0_term_right[0];
		

		Bbc[0] = Bb[0]-Bb_term_left[0]-Bb_term_right[0];
		if (Bbc0 != null) Bbc0[0] = Bb0[0]-Bb0_term_left[0]-Bb0_term_right[0];

		Bbc_term_left[0] = Bb_term_left[0];
		if (Bbc0_term_left != null) Bbc0_term_left[0] = Bb0_term_left[0];
		Bbc_term_right[0] = Bb_term_right[0];
		if (Bbc0_term_right != null) Bbc0_term_right[0] = Bb0_term_right[0];

		for (i = 1; i < blocknum; i++) {
			Bac[i] = Bac[i-1] + (Ba[i]-Ba_term_left[i]-Ba_term_right[i]);
			if (Bac0 != null) Bac0[i] = Bac0[i-1] + (Ba0[i]-Ba0_term_left[i]-Ba0_term_right[i]);
			
			Bac_term_left[i] = Bac_term_left[i-1] + Ba_term_left[i];
			if (Bac0_term_left != null)  Bac0_term_left[i] = Bac0_term_left[i-1] + Ba0_term_left[i];
			Bac_term_right[i] = Bac_term_right[i-1] + Ba_term_right[i];
			if (Bac0_term_right != null)  Bac0_term_right[i] = Bac0_term_right[i-1] + Ba0_term_right[i];

			Bbc[i] = Bbc[i-1] + (Bb[i]-Bb_term_left[i]-Bb_term_right[i]);
			if (Bbc0 != null)  Bbc0[i] = Bbc0[i-1] + (Bb0[i]-Bb0_term_left[i]-Bb0_term_right[i]);

			Bbc_term_left[i] = Bbc_term_left[i-1] + Bb_term_left[i];
			if (Bbc0_term_left != null) Bbc0_term_left[i] = Bbc0_term_left[i-1] + Bb0_term_left[i];
			Bbc_term_right[i] = Bbc_term_right[i-1] + Bb_term_right[i];
			if (Bbc0_term_right != null) Bbc0_term_right[i] = Bbc0_term_right[i-1] + Bb0_term_right[i];
		}
	}
}

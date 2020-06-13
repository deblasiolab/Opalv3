package opal.align;

public class CompressedAlignmentGraph {
	int[] shift;
	int M,N;
	public long[][] H, V, D;
	
	public CompressedAlignmentGraph(Aligner al, int[] nearlyOptimal_horizStart, int[] nearlyOptimal_horizEnd, boolean reverse) {
		this (al, nearlyOptimal_horizStart, nearlyOptimal_horizEnd, reverse, true);
	}

	public CompressedAlignmentGraph(Aligner al, int[] nearlyOptimal_horizStart, int[] nearlyOptimal_horizEnd, boolean reverse, boolean compress) {
		M = al.M;
		N = al.N;
		
		shift = new int[M+1];
		H = new long[M+1][];
		V = new long[M+1][];
		D = new long[M+1][];
		
		for (int i=0; i<=M; i++) {

			int len_i;

			if (!compress) {
				len_i = N+1;
				shift[i]=0;
			} else if (reverse) { 
				
//				int start = nearlyOptimal_horizStart[M-i];
				int start;
				if (i==M)
					start = nearlyOptimal_horizStart[0];
				else
					start = nearlyOptimal_horizStart[M-(i+1)];
				if (start==N) start = N-1;  //do this b/c the range for candidate ks will include this extra position
				if (start>0) start -= 1; // this is to allow edges in the subopt graph to go back one more
				
//				int end = nearlyOptimal_horizEnd[M-i];
				// need to keep more than the value above, because the reverse edges might come up from the row below
				// better tracking of the nearly-optimal paths could reduce this, so I'd only check horiz edges here 
				int end;
				if (i==0)
					end = nearlyOptimal_horizEnd[M];
				else
					end = nearlyOptimal_horizEnd[M-(i-1)]; // end of the next row
				if (end==0) end = 1; //do this b/c the range for candidate ks will include this extra position
				if (end<N) end += 1; // this is to allow edges in the subopt graph to go forward one more

				shift[i] = N - end;
				len_i = end-start + 1;
			} else {

				//int start = nearlyOptimal_horizStart[i];
				int start;
				if (i==0)
					start = nearlyOptimal_horizStart[0];
				else
					start = nearlyOptimal_horizStart[i-1];
				if (start == N) start = N-1; //do this b/c the range for candidate ks will include this extra position
				if (start>0) start -= 1;
				
				//int end = nearlyOptimal_horizEnd[i];
				int end;
				if (i==M) 
					end = nearlyOptimal_horizEnd[M];
				else {
					end = nearlyOptimal_horizEnd[i+1];
				}
				if (end==0) end = 1; //do this b/c the range for candidate ks will include this extra position
				if (end<N) end += 1;
				if (end<N) end += 1; // second time, to account for case in getOpenSubopt2
				
				shift[i] = start;
				len_i = end-start + 1;
			}
			
			H[i] = new long[len_i];
			V[i] = new long[len_i];
			D[i] = new long[len_i];

			for (int j=0; j<len_i; j++) {
				H[i][j] = al.H[i][j+shift[i]];
				V[i][j] = al.V[i][j+shift[i]];
				D[i][j] = al.D[i][j+shift[i]];		
			}
			
		}
	}

}

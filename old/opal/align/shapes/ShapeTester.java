package opal.align.shapes;

import java.util.ArrayList;
import java.util.Iterator;

abstract public class ShapeTester {
/*	ShapeMethods  (accept the upper bound and lower bound table)
	test dominance of one shape over another
	test shape for bound pruning.  requires stitching
	stitching (for bound pruning)
*/	
	long upperBound;
	long[][] lowerD, lowerH, lowerV;

	public ShapeTester (long upperBound, long[][] lowerD, long[][] lowerH, long[][] lowerV) {
		this.upperBound = upperBound;
		this.lowerD = reverseBoundArray(lowerD);
		this.lowerH = reverseBoundArray(lowerH);
		this.lowerV = reverseBoundArray(lowerV);
	}

	final private long[][] reverseBoundArray (long[][] toFlip){
		int X = toFlip.length;
		int Y = toFlip[0].length;
		long[][] flipped = new long[X][Y];
		for (int i=0; i<X; i++){
			for (int j=0; j<Y; j++){
				flipped[i][j] = toFlip[X-i-1][Y-j-1];
			}
		}
		return flipped;
	}
	
	public boolean boundPrune (Shape s) {
		/* see calcGapBounds for def of Hgap, Vgap, and Dgap */
		
		long[] H_V_D = new long[3]; 
		calcGapBounds(s, H_V_D);
		long Hgap = H_V_D[0];
		long Vgap = H_V_D[1];
		long Dgap = H_V_D[2];
	  
		/* if, for any of the three possible ways to begin an extension,
		 * the current cost plus a lowerbound on the extension to a global
		 * alignment (adjusted conservatively for gaps counted in both) is
		 * no worse than the upperbound on an optimal solution, then we
		 * cannot safely prune the shape.
		 */		

		/*
		Logger.stdErrLog(s.aPos + ", " + s.bPos + ": " + s.shapeCost );
		Logger.stdErrLog( " ..  D (" + (0-Dgap) + "+" +  lowerD[s.aPos][s.bPos] + ")" );
		Logger.stdErrLog( " ..  H (" + (0-Hgap) + "+" +  lowerH[s.aPos][s.bPos] + ")" );
		Logger.stdErrLogln( " ..  V (" + (0-Vgap) + "+" +  lowerV[s.aPos][s.bPos] + ")" );
	*/
		
		if ((s.shapeCost - Hgap + lowerH[s.aPos][s.bPos] <= upperBound ) ||
				(s.shapeCost - Vgap + lowerV[s.aPos][s.bPos] <= upperBound ) ||
				(s.shapeCost - Dgap + lowerD[s.aPos][s.bPos] <= upperBound )) {
				// it's possible that some extension of the shape will fall under the upper bound, so can't prune
			return false;
		} else {
			// otherwise the shape cannot be extended to an optimal solution and may be safely pruned.
			return true;
		}
	}
	
	/* Given shape s with associated cost c and a list of shapes, 
	 * dominance returns true if s is dominated by some shape in the 
	 * list and false otherwise.  All shapes in list that are dominated 
	 * by s are removed from list.
	 */
	 
	public boolean dominancePrune(Shape s, ArrayList<Shape> others) {

	  //check to see if s is dominated by a shape already in the list 
	  Iterator<Shape> it = others.iterator();
	  while (it.hasNext()) {
		  Shape t = (Shape)(it.next());
		  //check if s is dominated by current shape (t) from list  
		  if ( t.shapeCost <= s.shapeCost ) {//if that doesn't hold, there's no way t can dominate s 
			  if ( t.shapeCost + domGapBound(t, s) <= s.shapeCost) 
				  return true;
		  }
		
	     //reverse test: check if s dominates current shape (t) from list 
		  if ( s.shapeCost <= t.shapeCost ) {//otherwise, no way s can dominate t 
			  if ( s.shapeCost + domGapBound(s, t) <= t.shapeCost ) { 
				  it.remove();// t is dominated by s, so remove it from list 
				  break; //s dominates another, so nothing dominates it (dom is transitive)
			  }
		  }
	  }
	  
	  while (it.hasNext()) {
		  //If I'm here, it's because s dominated at least one shape t above.
		  //since dom is transitive, no need to check if s is dominated ... 
		  // but still want to see if it dominates anything else
		  Shape t = (Shape)(it.next());
		     //reverse test: check if s dominates current shape (t) from list 
		  if ( s.shapeCost <= t.shapeCost ) {//otherwise, no way s can dominate t 
			  if ( s.shapeCost + domGapBound(s, t) <= t.shapeCost ) {
				  it.remove(); // t is dominated by s, so remove it from list
			  }
		  }
	  }
	  return false;
	}
	
	/* procedure: domGapBound
	 *   Input:  shapes s and t
	 *   Output: cost, a int representing an upperbound on the cost of
	 *           gaps that could possibly be started in s but not in t,
	 *           over all possible extensions.
	 */
	
	abstract protected long domGapBound (Shape s, Shape t) ;
	
	/* Hgap, Vgap, and Dgap are upperbounds on the cost (unweighted? aka cost = number) of gaps
	 * that could be open in both the alignment so far and the
	 * computation of the extension (which was done via optimistic
	 * on the reverse of the input )
	 * This abstract class is implemented in either the Linear or Quadratic ShapeTester
	 * 
	 * The H_V_D array holds the values Hgap, Vgap and Dgap, and allows for easy transport
	 */
	abstract protected void calcGapBounds(Shape s, long[] H_V_D);
		
}


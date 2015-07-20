package opal.tree;


import opal.align.Aligner;
import opal.align.Alignment;

public abstract class Distance {
	//vars in the standard class
	int[][] costs;
	Aligner aligner;
	
	public Distance( Aligner al) {
		this.costs = al.config.cost.costs;
		aligner = al;
	}

	public abstract float calcDistance(Alignment A, Alignment B); 
}

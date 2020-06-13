package opal.align;


public abstract class ExactCountAligner extends Aligner {
		

	public ExactCountAligner() {
		super(true);
	}
	public ExactCountAligner( Alignment A, Alignment B) {
		super(A,B,true);
	}
	
	public ExactCountAligner( boolean pess) {
		super(pess);
	}
	
	public ExactCountAligner( Alignment A, Alignment B, boolean pess) {
		super(A,B,pess);
	}
	
	
}

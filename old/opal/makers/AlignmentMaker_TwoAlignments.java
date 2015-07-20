package opal.makers;

import java.text.NumberFormat;

import opal.IO.AlignmentWriter;
import opal.IO.ClustalWriter;
import opal.IO.FastaWriter;
import opal.IO.OpalLogWriter;

import com.traviswheeler.libs.LogWriter;
import opal.IO.SequenceConverter;
import opal.IO.SequenceFileReader;
import opal.IO.AlignmentWriter.OutputType;
import opal.align.Aligner;
import opal.align.Alignment;
import opal.align.ConsistencyAligner;
import opal.align.ExactCountAligner_Time;
import opal.align.PairAligner;
import opal.align.ProfileAligner;
import opal.align.Aligner.AlignmentType;


public class AlignmentMaker_TwoAlignments extends AlignmentMaker {

	Aligner al;
	Alignment alA;
	Alignment alB; 

	String fileA = "";
	String fileB = "";
	String[] namesA;
	String[] namesB;
	int K;
	int L;

	char [][] charsA;
	char [][] charsB;


	
	/* Align two alignments */
	final public int[][] buildAlignment ( String costName, int verbosity,  boolean toUpper) { 
		
		getAligner();				
		
//		Date date1 = new Date();
		al.align();
//		Date date2 = new Date();
//		long diff1 = date2.getTime() - date1.getTime();

	
		char[][] result = SequenceConverter.convertPathToCharAlignment (al.getPath(), charsA, charsB);

		 
		AlignmentWriter wr;
		if ( AlignmentWriter.outFormat == OutputType.fasta)
			wr = new FastaWriter(namesA, namesB, result, K, L, toUpper);
		else {//clustal
			wr = new ClustalWriter(namesA, namesB, result, K, L, toUpper);
			wr.setPath(al.getPath());
		}
		wr.write(outputWidth);
	
	
		if (verbosity > 0) {
			
			LogWriter.stdErrLogln("================================");
			LogWriter.stdErrLogln("A_file = " + fileA + "(" + K + " sequences)");
			LogWriter.stdErrLogln("B_file = " + fileB + "(" + L + " sequences)");
			if (null != outputName)
				LogWriter.stdErrLogln("output file = " + outputName);
	
			LogWriter.stdErrLogln("Cost matrix is " + costName); 
			
			printParams (alA); 
			
					//LogWriter.stdErrLogln("Estimated cost: " + est);
			if (showCost) {
				long totalCost = al.getTrueCost();
				int[] ids = new int[result.length];
				for (int i=0; i<result.length; i++) ids[i] = i;
				long cost = Aligner.calcCost(result, K, L, ids);
			
				LogWriter.stdErrLogln("A-to-B alignment cost:      " + NumberFormat.getInstance().format( cost ));
				LogWriter.stdErrLogln("All rows cost:  " + NumberFormat.getInstance().format( totalCost ) );
			}
			LogWriter.stdErrLogln("================================");
		}		
		return null;
	}
	
	public void initialize (String AFile, String structAFile, String BFile, String structBFile) { 
		// not implemented, yet
	}
	
	public void initialize (String AFile, String BFile) { 
		fileA = AFile;
		fileB = BFile;
	
		SequenceFileReader seqReaderA = new SequenceFileReader(AFile,false);						
		SequenceFileReader seqReaderB = new SequenceFileReader(BFile,false);
		charsA = seqReaderA.getSeqs();
		charsB = seqReaderB.getSeqs();
				
		initialize( SequenceConverter.convertSeqsToInts(charsA),seqReaderA.getNames(), 
					SequenceConverter.convertSeqsToInts(charsB), seqReaderB.getNames());
	}
	
	public void initialize (int[][] seqsA, String[] namesA, int[][] seqsB, String[] namesB) {

		alA = getAlignment( seqsA );
		alB = getAlignment( seqsB );		

		this.namesA = namesA;
		this.namesB = namesB;
		K = namesA.length;
		L = namesB.length;

	}

	
	protected Alignment getAlignment (int[][] seqs) {
		int[] ids = new int[seqs.length];
		for (int i=0; i<seqs.length; i++) ids[i] = i;
		
		return Alignment.buildNewAlignment(seqs, ids);

	}

	
	protected void getAligner() {		

		if (/*alignMethod == JOpal.ALIGN_EXACT &&*/ alA.K == 1 && alB.K == 1) {
			al = new PairAligner(alA, alB);
		} else if (Aligner.alignmentMethod == Aligner.AlignmentType.profile || 
				(Aligner.alignmentMethod == AlignmentType.mixed && alA.K*alB.K >= Aligner.mixedAlignmentCutoff)) {
			al = new ProfileAligner(alA, alB, true);
		} else if (Aligner.alignmentMethod == Aligner.AlignmentType.exact || 
				(Aligner.alignmentMethod == AlignmentType.mixed && alA.K*alB.K < Aligner.mixedAlignmentCutoff)) {
			al = new ExactCountAligner_Time(alA, alB);		
		} 
		
	}
	
	protected void printParams (Alignment example /* used in subclass*/) {
		LogWriter.stdErrLogln("gamma is " + Aligner.gamma + " and lambda is " + Aligner.lambda);
		if (Aligner.gammaTerm != Aligner.gamma  ||  Aligner.lambdaTerm != Aligner.lambda)
			LogWriter.stdErrLogln("gamma_term is " + Aligner.gammaTerm + " and lambda_term is " + Aligner.lambdaTerm);
		LogWriter.stdErrLogln("Solution alignment length is " + al.getPath().size());

		LogWriter.stdErrLogln("Alignment method is : ");
			if (Aligner.alignmentMethod == Aligner.AlignmentType.profile) {
				LogWriter.stdErrLogln("Profile alignment (pessimistic heuristic)");
			} else if (Aligner.alignmentMethod == Aligner.AlignmentType.exact) {
				LogWriter.stdErrLogln("Exact alignment");
			} else if (Aligner.alignmentMethod == Aligner.AlignmentType.mixed) {
				LogWriter.stdErrLogln("Exact for small alignments (" 
						+ Aligner.mixedAlignmentCutoff 	+ " pairs), profile for large");
			}
	}
	
}

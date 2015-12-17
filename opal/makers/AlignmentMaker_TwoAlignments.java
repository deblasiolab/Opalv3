package opal.makers;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.text.NumberFormat;

import opal.IO.AlignmentWriter;
import opal.IO.ClustalWriter;
import opal.IO.FastaWriter;
import com.traviswheeler.libs.LogWriter;

import opal.IO.SequenceFileReader;
import opal.IO.AlignmentWriter.OutputType;
import opal.align.Aligner;
import opal.align.Alignment;
import opal.align.ExactCountAligner_Time;
import opal.align.PairAligner;
import opal.align.ProfileAligner;
import opal.align.Aligner.AlignmentType;
import opal.exceptions.GenericOpalException;
import opal.IO.Configuration;
import opal.IO.Inputs;


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
	Configuration conf;
	Inputs in;

	public char[][] result;
	
	/* Align two alignments */
	final public int[][] buildAlignment () { 
		
		getAligner();				
		
//		Date date1 = new Date();
		al.align();
//		Date date2 = new Date();
//		long diff1 = date2.getTime() - date1.getTime();

	
		result = conf.sc.convertPathToCharAlignment (al.getPath(), charsA, charsB);
	
		//al = null;
		return conf.sc.convertSeqsToInts(result);
		//return null;
	}
	
	public boolean printOutput(int[][] result_i, String fname){
		PrintStream stdout = System.out;
		if(fname!=null){
			try{
				PrintStream out = new PrintStream(new FileOutputStream(fname));
				System.setOut(out);
			}catch(FileNotFoundException e){
				throw new GenericOpalException(e.toString());
			}
		}
		
		AlignmentWriter wr;
		if ( AlignmentWriter.outFormat == OutputType.fasta)
			wr = new FastaWriter(namesA, namesB, result, K, L, in.toUpper);
		else {//clustal
			wr = new ClustalWriter(namesA, namesB, result, K, L, in.toUpper);
			wr.setPath(al.getPath());
		}
		wr.write(outputWidth);
	
	
		if (in.verbosity > 0) {
			
			LogWriter.stdErrLogln("================================");
			LogWriter.stdErrLogln("A_file = " + fileA + "(" + K + " sequences)");
			LogWriter.stdErrLogln("B_file = " + fileB + "(" + L + " sequences)");
			if (null != fname)
				LogWriter.stdErrLogln("output file = " + fname);
	
			LogWriter.stdErrLogln("Cost matrix is " + conf.cost.costName); 
			
			printParams (alA); 
			
					//LogWriter.stdErrLogln("Estimated cost: " + est);
			if (showCost) {
				long totalCost = al.getTrueCost();
				int[] ids = new int[result.length];
				for (int i=0; i<result.length; i++) ids[i] = i;
				long cost = Aligner.calcCost(result, K, L, ids, conf);
			
				LogWriter.stdErrLogln("A-to-B alignment cost:      " + NumberFormat.getInstance().format( cost ));
				LogWriter.stdErrLogln("All rows cost:  " + NumberFormat.getInstance().format( totalCost ) );
			}
			LogWriter.stdErrLogln("================================");
		}
		
		if(fname!=null){
			System.setOut(stdout);
		}
		
		return true;
	}
	
	public void initialize (String AFile, String structAFile, String BFile, String structBFile) { 
		// not implemented, yet
	}
	
	public void initialize (Configuration config, Inputs inp) { 
		conf = config;
		in = inp;
		fileA = in.fileA;
		fileB = in.fileB;
	
		SequenceFileReader seqReaderA = new SequenceFileReader(in.fileA,false);						
		SequenceFileReader seqReaderB = new SequenceFileReader(in.fileB,false);
		charsA = seqReaderA.getSeqs();
		charsB = seqReaderB.getSeqs();
				
		initialize( conf.sc.convertSeqsToInts(charsA),seqReaderA.getNames(), 
					conf.sc.convertSeqsToInts(charsB), seqReaderB.getNames(),
					conf, in);
	}
	
	public void initialize (int[][] seqsA, String[] namesA, int[][] seqsB, String[] namesB, Configuration c, Inputs i) {

		this.conf = c;
		this.in = i;
		
		alA = getAlignment( seqsA );
		alB = getAlignment( seqsB, seqsA.length );		

		this.namesA = namesA;
		this.namesB = namesB;
		K = namesA.length;
		L = namesB.length;

	}

	protected Alignment getAlignment (int[][] seqs) {
		return getAlignment(seqs, 0);
	}
	
	protected Alignment getAlignment (int[][] seqs, int startIndex) {
		int[] ids = new int[seqs.length];
		for (int i=0; i<seqs.length; i++) ids[i] = i;
		
		return Alignment.buildNewAlignment(seqs, ids, conf, startIndex);

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
		
		al.config = conf;
		
	}
	
	protected void printParams (Alignment example /* used in subclass*/) {
		LogWriter.stdErrLogln("gamma is " + conf.gamma + " and lambda is " + conf.lambda);
		if (conf.useLeftTerminal && conf.useRightTerminal && (conf.leftGammaTerm() != conf.gamma  ||  conf.leftLambdaTerm() != conf.lambda))
			LogWriter.stdErrLogln("gamma_term is " + conf.leftGammaTerm() + " and lambda_term is " + conf.leftLambdaTerm());
		else if (conf.useLeftTerminal && (conf.leftGammaTerm() != conf.gamma  ||  conf.leftLambdaTerm() != conf.lambda))
			LogWriter.stdErrLogln("left gamma_term is " + conf.leftGammaTerm() + " and left lambda_term is " + conf.leftLambdaTerm());
		else if (conf.useRightTerminal && (conf.rightGammaTerm() != conf.gamma  ||  conf.rightLambdaTerm() != conf.lambda))
			LogWriter.stdErrLogln("right gamma_term is " + conf.rightGammaTerm() + " and right lambda_term is " + conf.rightLambdaTerm());
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

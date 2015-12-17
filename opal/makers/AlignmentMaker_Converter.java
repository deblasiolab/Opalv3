package opal.makers;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;

import opal.IO.AlignmentWriter;
import opal.IO.ClustalWriter;
import opal.IO.FastaWriter;
import opal.IO.OpalLogWriter;

import com.traviswheeler.libs.LogWriter;

import opal.IO.SequenceConverter;
import opal.IO.SequenceFileReader;
import opal.IO.StructureFileReader;
import opal.IO.AlignmentWriter.OutputType;
import opal.align.Aligner;
import opal.align.ProfileAligner;
import opal.exceptions.GenericOpalException;
import opal.IO.Configuration;
import opal.IO.Inputs;


public class AlignmentMaker_Converter extends AlignmentMaker {

	String[] names;
	char[][] seqs;
	int K;
	String file;
	Configuration conf;
	Inputs in;
			
	public void initialize (String fileA, String structFileA, String fileB, String structFileB) { 
		/* nothing */
	}

	public void initialize (Configuration c, Inputs inp) { 
		in = inp;
		conf = c;
		SequenceFileReader seqReader = new SequenceFileReader(in.fileA, false);		
		names = seqReader.getNames();
		seqs = seqReader.getSeqs();
		K = names.length;
		this.file = in.fileA;
	}

	
	final public int[][] buildAlignment () {
		return conf.sc.convertSeqsToInts(seqs);
		//return null;
	}

	public boolean printOutput(int[][] num, String fname){
		PrintStream stdout = System.out;
		if(fname!=null){
			try{
				PrintStream out = new PrintStream(new FileOutputStream(fname));
				System.setOut(out);
			}catch(FileNotFoundException e){
				throw new GenericOpalException(e.toString());
			}
		}
		
		if (in.verbosity > 1) {
			AlignmentWriter wr;
			if ( AlignmentWriter.outFormat == OutputType.fasta)
				wr = new FastaWriter(names, seqs, K, in.toUpper);
			else {//CLUSTAL 
				wr = new ClustalWriter(names, seqs, K, in.toUpper);
			}			
			wr.write(outputWidth);			
		}		

		int[] ids = new int[names.length];
		for (int i=0; i<names.length; i++) ids[i] = i;

		long cost = Aligner.calcCost(seqs, ids, conf); 

		
		LogWriter.stdErrLogln("================================");
		LogWriter.stdErrLogln("input file = " + file + "(" + K + " sequences)");
		LogWriter.stdErrLogln("Cost matrix is " + conf.cost.costName); 
		LogWriter.stdErrLogln("gamma is " + conf.gamma + " and lambda is " + conf.lambda);
		if (conf.gammaTerm != conf.gamma  ||  conf.lambdaTerm != conf.lambda)
			LogWriter.stdErrLogln("gamma_term is " + conf.gammaTerm + " and lambda_term is " + conf.lambdaTerm);
		LogWriter.stdErrLogln("Alignment length is " + seqs[0].length);
		LogWriter.stdErrLogln("Alignment cost:      " + cost);
		LogWriter.stdErrLogln("================================");
		
		if(fname!=null){
			System.setOut(stdout);
		}
		return true;
	}
	
	/* Align two alignments */
	final public int[][] buildAlignment ( String fileA, String structFileA, 
			String fileB, String structFileB, String costName, boolean quiet, boolean toUpper) { 
		//bad
		return null;
	}
	
	protected Aligner getAligner() {
		return new ProfileAligner();
	}


}

package opal.makers;

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


public class AlignmentMaker_Converter extends AlignmentMaker {

	String[] names;
	char[][] seqs;
	int K;
	String file;
			
	public void initialize (String fileA, String structFileA, String fileB, String structFileB) { 
		/* nothing */
	}

	public void initialize (String file, String structFile) { 
		SequenceFileReader seqReader = new SequenceFileReader(file, false);		
		names = seqReader.getNames();
		seqs = seqReader.getSeqs();
		K = names.length;
		this.file = file;
	}

	
	final public int[][] buildAlignment ( String costName, int verbosity, boolean toUpper) {
		if (verbosity > 0) {
			AlignmentWriter wr;
			if ( AlignmentWriter.outFormat == OutputType.fasta)
				wr = new FastaWriter(names, seqs, K, toUpper);
			else {//CLUSTAL 
				wr = new ClustalWriter(names, seqs, K, toUpper);
			}			
			wr.write(outputWidth);			
		}		

		int[] ids = new int[names.length];
		for (int i=0; i<names.length; i++) ids[i] = i;

		long cost = Aligner.calcCost(seqs, ids); 

		
		LogWriter.stdErrLogln("================================");
		LogWriter.stdErrLogln("input file = " + file + "(" + K + " sequences)");
		LogWriter.stdErrLogln("Cost matrix is " + costName); 
		LogWriter.stdErrLogln("gamma is " + Aligner.gamma + " and lambda is " + Aligner.lambda);
		if (Aligner.gammaTerm != Aligner.gamma  ||  Aligner.lambdaTerm != Aligner.lambda)
			LogWriter.stdErrLogln("gamma_term is " + Aligner.gammaTerm + " and lambda_term is " + Aligner.lambdaTerm);
		LogWriter.stdErrLogln("Alignment length is " + seqs[0].length);
		LogWriter.stdErrLogln("Alignment cost:      " + cost);
		LogWriter.stdErrLogln("================================");

		return null;
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

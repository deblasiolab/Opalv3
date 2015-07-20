package opal.makers;

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
import opal.align.PairAligner_SplitGamma;
import opal.align.PairSuboptimalityMatrices;

import opal.IO.Configuration;
import opal.IO.Inputs;


public class AlignmentMaker_SuboptimalityTester extends AlignmentMaker {
		
	String[] names;
	char [][] chars;
	int [][] seqs;
	int K;
	String file;
	Configuration conf;
	Inputs in;

	public AlignmentMaker_SuboptimalityTester() {
		super();
	}

	public void initialize (String fileA, String structFileA, String fileB, String structFileB) { 
		/* nothing */
	}

	public void initialize (Configuration c, Inputs inp) { 
		conf = c;
		in = inp;
		SequenceFileReader seqReader = new SequenceFileReader(in.fileA, false );		
		names = seqReader.getNames();
		chars = seqReader.getSeqs();
		seqs = conf.sc.convertSeqsToInts( chars );				
		K = names.length;
		this.file = in.fileA;
	}

	
	final public int[][] buildAlignment () {


		//Aligner al = getAligner();
		Alignment alA, alB;
		char[][] charsA = new char[1][];
		char[][] charsB = new char[1][];
		String[] namesA = new String[1];
		String[] namesB = new String[1];
		
		PairSuboptimalityMatrices cont;
		
		for (int i=0; i<K-1; i++ ) {
			alA = Alignment.buildNewAlignment(seqs[i], i, conf);
			charsA[0] = chars[i];
			namesA[0] = names[i];
			for (int j=i+1; j<K; j++){
				alB = Alignment.buildNewAlignment(seqs[j], j, conf);
				charsB[0] = chars[j];
				namesB[0] = names[j];

				LogWriter.stdOutLogln(names[i] + ", " + names[j] + " :");
				cont = new PairSuboptimalityMatrices (alA, alB, conf);
				
				if (in.verbosity>1) {
					char[][] result = conf.sc.convertPathToCharAlignment (cont.getAlignment().getPath(), charsA, charsB);
					AlignmentWriter wr;
					if ( AlignmentWriter.outFormat == OutputType.fasta)
						wr = new FastaWriter(namesA, namesB, result, 1, 1, in.toUpper);
					else {//CLUSTAL 
						wr = new ClustalWriter(namesA, namesB, result, 1, 1, in.toUpper);
						wr.setPath(cont.getAlignment().getPath());
					}
					wr.write(outputWidth);
				}		

				cont = null;

				LogWriter.stdOutLogln("\n------------------------\n");
			}
		}

		return null;
	}

	public boolean printOutput(int[][] input, String fname){
		return true;
	}
	
	/* Align two alignments */
	final public int[][] buildAlignment ( String fileA, String structFileA, 
			String fileB, String structFileB, String costName, boolean quiet, boolean toUpper) { 
		//bad
		return null;
	}
	
	protected Aligner getAligner() {
		return new PairAligner_SplitGamma();
	}

}

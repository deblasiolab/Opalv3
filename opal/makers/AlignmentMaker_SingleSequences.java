package opal.makers;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.io.FileOutputStream;

import opal.polish.*;
import opal.polish.Polisher.PolishType;
import opal.tree.*;
import opal.tree.Tree.DistanceType;
import opal.IO.*;
import opal.IO.AlignmentWriter.OutputType;
import opal.align.*;
import opal.align.Aligner.AlignmentType;
import opal.exceptions.GenericOpalException;

import com.traviswheeler.libs.LogWriter;

import opal.IO.SequenceConverter;

public class AlignmentMaker_SingleSequences extends AlignmentMaker {

	static String treeOutFile;
	public static boolean consistency = false;
	public static int consistency_subs_increase = 40;
	public static int kmerTreeCutoff = 20; 
	public static boolean writeResult = true;
	
	int[][] seqs;
	char [][] chars;

	Aligner al = null;
	
	int verbosity;
	String costName;
	public String file;
	boolean toUpper;
	TreeNode latestNode = null; 
	String[] names = null;
	public Polisher polisher = null;
	public Tree tree = null;
	private Distance dist = null;
	
	private int currIteration;
	
	public boolean alignmentComplete = false;
	
	private Configuration conf;
	private Inputs in;
	
	public AlignmentMaker_SingleSequences( ) {
		super();
		currIteration = 1;
	}
	
	public AlignmentMaker_SingleSequences(int iter) {
		super();
		currIteration = iter;
	}	
	
	public static void setTreeOutFile (String s) {
		treeOutFile = s;
	}

	final public void buildTree ( Configuration c, Inputs inp) { 
		conf = c;
		in = inp;
		Alignment[] alignments = getAlignments(seqs);				

		getAligner();

		if (Tree.distanceType == DistanceType.kmer ||
				(currIteration == Tree.iterations && Tree.distanceType == DistanceType.kmer_normcost)) {
			for (int i=0; i<alignments.length; i++)
				alignments[i].countKmers();
		}
		getDistance(currIteration, conf);
		Tree tree = getTree(alignments, currIteration);
		
		if (Tree.justPWDists) 
			return;
		
		if (treeOutFile != null ){
			StringBuffer sb = new StringBuffer();			
			tree.buildTreeString(sb, names);
						
			PrintStream out = null;
	    	try {
				out = new PrintStream(new FileOutputStream(treeOutFile));
	    	} catch (Exception e) {
	    		LogWriter.stdErrLogln("Error creating file for tree output: " +treeOutFile);
	    		throw new GenericOpalException(e.getMessage());		            		
	    	}
			
			out.print(sb);
			out.close();

			if (verbosity>1) LogWriter.stdErrLogln("output tree file = " + treeOutFile);
		}

	}	
	
	/* Align a set of sequences*/
	final public int[][] buildAlignment () { 

		this.verbosity = in.verbosity;
		this.costName = conf.cost.costName;
		this.toUpper = in.toUpper;
		
		Alignment[] alignments;

//		if (Tree.distanceType == DistanceType.kmer_normcost && currIteration<Tree.iterations) {
		if ((Tree.distanceType == DistanceType.kmer_normcost || Tree.distanceType == DistanceType.normcost) && currIteration<Tree.iterations) {
			
			AlignmentMaker_SingleSequences maker = new AlignmentMaker_SingleSequences(currIteration+1);
			maker.initialize(seqs, names, conf, in);
			int[][] res = maker.buildAlignment(); 
			alignments = new Alignment[1];
			
			int[] ids = new int[res.length];
			for (int i=0; i<res.length; i++) ids[i] = i;
			alignments[0] = Alignment.buildNewAlignment(res, ids, conf, in);

			if (verbosity > 1) {
				if (initAlignmentProvided) {
					LogWriter.stdErrLog("\nRead in input alignment. Treated as draft #" + (Tree.iterations-currIteration) + " of " + Tree.iterations + "  " );					
				} else {
					LogWriter.stdErrLog("\nFinished draft alignment #" + (Tree.iterations-currIteration) + " of " + Tree.iterations + "  " );
				}
				if (showCost) {
					long cost = Aligner.calcCost(res, ids, conf, in); 
					LogWriter.stdErrLogln("Alignment cost:      " + NumberFormat.getInstance().format( cost ));
				}
				LogWriter.stdErrLogln("Entering next phase\n======================\n");
			}
		} else {
			alignments = getAlignments(seqs);
		}
		
		if (currIteration==Tree.iterations && initAlignmentProvided && !Polisher.justPolish) {
			return alignments[0].seqs;  
		}
		
		getAligner();
		if (Polisher.justPolish) {
			tree = getTree(alignments, currIteration);
			latestNode = tree.getRoots().get(0);
			latestNode.setAlignment(alignments[0]);
		} else {
			
			if (Tree.distanceType == DistanceType.kmer ||
					(currIteration == Tree.iterations && Tree.distanceType == DistanceType.kmer_normcost)) {
				for (int i=0; i<alignments.length; i++)
					alignments[i].countKmers();
			}
			getDistance(currIteration, conf);
				
			int numSeqsForLog =  currIteration<Tree.iterations ? alignments[0].K : alignments.length;
			
			if (verbosity > 1) {
				LogWriter.stdErrLogln("Calculating pairwise distances for " + numSeqsForLog + " sequences");
			}
			tree = getTree(alignments, currIteration);
	
			if (verbosity > 1)  {
				LogWriter.stdErrLogln("\nPairwise distances calculated.\n");
				
				LogWriter.stdErrLogln("Opal will now merge sequences into increasingly large alignments.");
				LogWriter.stdErrLogln("A total of " + (numSeqsForLog-1) + " merges will be performed for the current iteration.");
				if (verbosity == 2) 
					LogWriter.stdErrLogln("After each merge, a dot will be printed:");
				
			}
	
	
			tree.prepareToMerge();
			
			if (treeOutFile != null ){
				StringBuffer sb = new StringBuffer();
				tree.buildTreeString(sb, names);
				
				PrintStream out = null;
	        	try {
	    			out = new PrintStream(new FileOutputStream(treeOutFile));
	        	} catch (Exception e) {
	        		LogWriter.stdErrLogln("Error creating file for tree output: " +treeOutFile);
	        		throw new GenericOpalException(e.getMessage());            		
	        	}
				
				out.print(sb);
				out.close();
			}
	
	
			int cnt = 1;
			while (tree.mergesRemaining() > 0) {
				latestNode = tree.mergeNext();
				
				if (verbosity == 2) {	
					LogWriter.stdErrLog(".");
					if ( cnt%20 == 0) LogWriter.stdErrLogln("  " + cnt + " complete");
					cnt++;
				}
	
				//possibly polish the just-formed group, using tree.getLastAlignment();  
	
	/*			LogWriter.stdOutLogln("*********\nmultiple alignment: \n**************\n");
				ArrayList inputOrder = latestNode.leafOrderFromInput;
				int K = inputOrder.size();
				int[][] seqs = latestNode.alignment.seqs;
				int[][] reorderedSeqs = new int[K][];
				for (int i=0; i<K; i++) {
					reorderedSeqs[((Integer)inputOrder.get(i)).intValue()] = seqs[i];
				}		
				char[][] result = Aligner.seqConv.convertIntArrayToCharAlignment(reorderedSeqs,seqReader.getSeqs());
				FastaWriter wr = new FastaWriter(names, result, K, toUpper);			
				wr.write(outputWidth);
				LogWriter.stdOutLogln("\n**************\n");
	*/		
			}
			
			if (verbosity > 1)  
				LogWriter.stdErrLogln("\n");

			alignmentComplete = true;

		}
		
		
		//move the costs back to originals, for polishing and output
		if (AlignmentMaker_SingleSequences.consistency) {
			conf.cost.multiplyCosts((float)0.5);
			conf.cost.increaseCosts(-consistency_subs_increase); // arbitrary number ... just trying to make w-w identities non-zero cost, and high enough to avoid rounding effects
			conf.multiplyGapCosts((float)0.5);
			conf.increaseGapCosts(-consistency_subs_increase/2);
		}

		
		//then polish latestNode? 
		if (Polisher.polishMethod != null && 
				(Polisher.polishIterations >0 || Polisher.polishMethod == Polisher.PolishType.exhaust_twocut)) {
			
			int polishIters = Polisher.polishIterations;
			if (currIteration > 1) {
				if (Polisher.polishMethod == PolishType.exhaust_twocut) 
					polishIters = (int)Math.ceil(polishIters/Math.pow(2, currIteration-1));
				else
					polishIters = (int)Math.ceil(polishIters/Math.pow(5, currIteration-1));
			}
			
			
			getPolishAligner();
	
			if (Polisher.polishMethod == Polisher.PolishType.exhaust_threecut) {
				if (verbosity>1) LogWriter.stdErrLogln("First: exhaustive polishing phase");
				polisher = new ExhaustiveTwoCutPolisher(latestNode, al, verbosity);
				//polisher.polish(Polisher.polishIterations);
				polisher.polish(polishIters);
				if (verbosity>1) LogWriter.stdErrLogln("Next: three-cut polishing phase");
				polisher = new RandomThreeCutPolisher(latestNode, al, verbosity);
			} else if (Polisher.polishMethod == Polisher.PolishType.exhaust_twocut) 
				polisher = new ExhaustiveTwoCutPolisher(latestNode, al, verbosity);
			else if (Polisher.polishMethod == Polisher.PolishType.random_tree_twocut) 
				polisher = new RandomTreeTwoCutPolisher(latestNode, al, verbosity);
			else if (Polisher.polishMethod == Polisher.PolishType.random_twocut) 
				polisher = new RandomTreeTwoCutPolisher(latestNode, al, verbosity);
			else // if (Polisher.polishMethod == Polisher.PolishType.random_threecut) 
				polisher = new RandomThreeCutPolisher(latestNode, al, verbosity);
		
//			polisher.polish(Polisher.polishIterations);
			polisher.polish(polishIters);
		}		
		
		
		
		int[][] reordered = reorderSeqs();
		//if (currIteration == 1 &&  writeResult) printOutput(reordered);

		alignments = null;
		return reordered;
	}

	public int[][] reorderSeqs () {
		if (null == seqs || null == latestNode )
			return null;
		
		// because of merge tree, sequences will be scrambled up
		int[][] reorderedSeqs;
		int K = seqs.length;
		if (AlignmentMaker.outputOrderMethod == OutputOrderType.input) {
			ArrayList<Integer> inputOrder = latestNode.leafOrderFromInput;
			seqs = latestNode.alignment.seqs;
			reorderedSeqs = new int[K][];
			for (int i=0; i<K; i++) {
				reorderedSeqs[inputOrder.get(i).intValue()] = seqs[i];
			}
		} else {
			reorderedSeqs = seqs;
		}
		
		return reorderedSeqs;
	}
	
	public boolean printOutput (int[][] reorderedSeqs, String fname, boolean printRealignmentLines) {

		PrintStream stdout = System.out;
		PrintStream out = null;
		if(fname!=null){
			try{
				java.io.File file = new java.io.File(fname);
				if(file.getParentFile()!=null)file.getParentFile().mkdirs();
				out = new PrintStream(new FileOutputStream(file));
				System.setOut(out);
			}catch(FileNotFoundException e){
				throw new GenericOpalException(e.toString());
			}
		}
		
		if (seqs == null || (Polisher.polishMethod != null && polisher == null && seqs.length > 2))
			return false;
		
		int K = seqs.length;
		char[][] result = conf.sc.convertIntArrayToCharAlignment(reorderedSeqs,chars);
		 
		AlignmentWriter wr;
		if ( AlignmentWriter.outFormat == OutputType.fasta)
			wr = new FastaWriter(names, result, K, toUpper);
		else {//CLUSTAL 
			wr = new ClustalWriter(names, result, K, toUpper);
			//wr.setPath(al.getPath());
		}			
		wr.write(outputWidth);
			
		if (verbosity>0) {
		
			LogWriter.stdErrLogln("\n================================");
			LogWriter.stdErrLogln("input file = " + file + " (" + K + " sequences)");
		}
		if (verbosity>-1) {
			if (null != fname)
				LogWriter.stdErrLogln("output file = " + fname);
		}
		
		if (verbosity>0) {
			if (null != EnteredTree.treeFile  )
				LogWriter.stdErrLogln("input tree file = " + EnteredTree.treeFile);
	
			if (null != treeOutFile  )
				LogWriter.stdErrLogln("output tree file = " + treeOutFile);
			
			LogWriter.stdErrLogln("Input is treated as a " + (conf.cost.isDNA ? "DNA" : "protein") + " sequence");
				
			LogWriter.stdErrLogln("Cost matrix is " + costName); 
			printParams(result[0].length, latestNode.alignment);
						
			int[] ids = new int[result.length];
			for (int i=0; i<result.length; i++) ids[i] = i;
			if (showCost) {
				long cost = Aligner.calcCost(result, ids, conf, in); 
				LogWriter.stdErrLogln("Alignment cost:      " + NumberFormat.getInstance().format( cost ));
			}
			if(printRealignmentLines){
				LogWriter.stdErrLog(conf.realignmentLog);
			}
			LogWriter.stdErrLogln("================================");
		}
		if(fname!=null){
			out.close();
			System.setOut(stdout);
		}
		return true;
	}
	
	protected void printParams (int length, Alignment example /* used in subclass*/) {
		if (verbosity>0) {
			LogWriter.stdErrLogln("gamma is " + conf.gamma + " and lambda is " + conf.lambda);
			if (conf.useLeftTerminal && conf.useRightTerminal)
				LogWriter.stdErrLogln("gamma_term is " + conf.leftGammaTerm() + " and lambda_term is " + conf.leftLambdaTerm());
			else if (conf.useLeftTerminal)
				LogWriter.stdErrLogln("left gamma_term is " + conf.leftGammaTerm() + " and left lambda_term is " + conf.leftLambdaTerm());
			else if (conf.useRightTerminal)
				LogWriter.stdErrLogln("right gamma_term is " + conf.rightGammaTerm() + " and right lambda_term is " + conf.rightLambdaTerm());
				
			LogWriter.stdErrLogln("Solution alignment length is " + length);
	
			if (!Polisher.justPolish) {
				LogWriter.stdErrLog("Alignment method is : ");
				if (Aligner.alignmentMethod == Aligner.AlignmentType.profile) {
					LogWriter.stdErrLog("Profile alignment (pessimistic heuristic)");
				} else if (Aligner.alignmentMethod == Aligner.AlignmentType.exact) {
					LogWriter.stdErrLog("Exact alignment");
				} else if (Aligner.alignmentMethod == Aligner.AlignmentType.mixed) {
					LogWriter.stdErrLog("Exact for small alignments (K*L<" 
							+ Aligner.mixedAlignmentCutoff 	+ "), profile for large");
				}
				if (consistency) {
					LogWriter.stdErrLog(" (with consistency modified costs)");
				}
				LogWriter.stdErrLog("\n");
		
				LogWriter.stdErrLog("Distance type: ");
				if (Tree.distanceType == DistanceType.normcost) {
					LogWriter.stdErrLogln(" pairwise normalized alignment cost");
				} else if (Tree.distanceType == DistanceType.pctid) {
					LogWriter.stdErrLogln(" pairwise alignment percent identity ");
				} else if (Tree.distanceType == DistanceType.kmer_normcost) {
					LogWriter.stdErrLog(" kmer-tree first pass, normalized alignment cost for " + 
							(Tree.iterations - 1 ) + " extra ");
					if (Tree.iterations==2)
						LogWriter.stdErrLogln("pass");
					else 
						LogWriter.stdErrLogln("passes");
	
				}
			}
			
			if (Polisher.polishMethod != null ) {
				if (Polisher.polishMethod == Polisher.PolishType.exhaust_twocut) { 
					LogWriter.stdErrLogln("Exhaustive two-cut polishing was performed, making " + 
							((ExhaustiveTwoCutPolisher)polisher).passCnt + " pass" + 
							(((ExhaustiveTwoCutPolisher) polisher).passCnt>1?"es":"") +" over the tree");
				} else if (Polisher.polishMethod == Polisher.PolishType.random_tree_twocut) {
					LogWriter.stdErrLogln("Random tree-based two-cut polishing was iterated " + 
							Polisher.polishIterations + " times");
				} else if (Polisher.polishMethod == Polisher.PolishType.random_twocut) {
					LogWriter.stdErrLogln("Random two-cut polishing was iterated " + 
							Polisher.polishIterations + " times");
				} else {// if (Polisher.polishMethod == Polisher.PolishType.random_threecut) 
					LogWriter.stdErrLogln("Random three-cut polishing was iterated " + 
													Polisher.polishIterations + " times");
				}
				
				LogWriter.stdErrLog("Polish method is : ");
				if (Polisher.polishAligmentMethod == Aligner.AlignmentType.profile) {
					LogWriter.stdErrLog("Profile alignment (pessimistic heuristic)");
				} else if (Polisher.polishAligmentMethod == Aligner.AlignmentType.exact) {
					if (Polisher.polishMethod == PolishType.exhaust_twocut) {
						LogWriter.stdErrLog("Exact alignment on final pass, after intial passes that are \n"+
								"exact for small alignments (K*L<" 	+ Aligner.mixedAlignmentCutoff 	+ 
								") and profile for large." );
					} else {	
						if (Polisher.polishIterations_exact > 0 )
							LogWriter.stdErrLog("Exact alignment (" + Polisher.polishIterations_exact + " reps) after " +
								(Polisher.polishIterations - Polisher.polishIterations_exact ) + " reps of \n" +
								"exact for small alignments (K*L<" 	+ Aligner.mixedAlignmentCutoff 	+ 
								") and profile for large.");
					}
				} else if (Polisher.polishAligmentMethod == Aligner.AlignmentType.mixed) {
					if (Polisher.polishMethod == PolishType.exhaust_twocut) {
						LogWriter.stdErrLog("Exact for small alignments (K*L<" 
								+ Aligner.mixedAlignmentCutoff 	+ "), profile for large, for all but the final pass" );
					} else {
						LogWriter.stdErrLog("Exact for small alignments (K*L<" 
								+ Aligner.mixedAlignmentCutoff 	+ "), profile for large." );
/*						
						if (Polisher.polishIterations_exact > 0 )
							LogWriter.stdErrLog("Exact for small alignments (K*L<" 
									+ Aligner.mixedAlignmentCutoff 	+ "), profile for large (" 
									+ Polisher.polishIterations_exact + " reps) after " );
						LogWriter.stdErrLog( "profile/exact (" + (Polisher.polishIterations - Polisher.polishIterations_exact ) + " reps)");
						*/
					}
				}
				LogWriter.stdErrLog("\n");				
			}
			if (conf.useStructure) {
				LogWriter.stdErrLogln("Using structure model: " + conf.modelType);
			}

		}

	}

	public void initialize (String fileA, String structFileA, String fileB, String structFileB) { 
			/* nothing */
	}
	
	public void initialize (int[][] seqs, String[] names, Configuration c, Inputs inp) {
		this.in = inp;
		this.conf = c;
		this.seqs = seqs;
		this.names = names;
		
		if (Polisher.polishMethod != null &&  Polisher.polishAligmentMethod == null) 
			Polisher.polishAligmentMethod = Aligner.alignmentMethod;

		int K = seqs.length;		
		if ( K==2 ) { 
			Polisher.polishIterations = 0;
		} else if (Polisher.polishIterations == -1) {

			//default number of polish iterations ... not too many
			int mult = 5;
			if (Polisher.polishAligmentMethod == AlignmentType.profile ) {
				mult = 10;
			} else if (Polisher.polishAligmentMethod == AlignmentType.mixed ) {
				mult = 8;
			}
			
			if (Polisher.polishMethod == PolishType.random_tree_twocut ||
					Polisher.polishMethod == PolishType.random_twocut) {
				mult *= 2;				
			}

			if (Polisher.polishMethod != PolishType.exhaust_twocut) 
				Polisher.polishIterations = mult*(K-3) + 10;
			
			if ( Polisher.polishIterations > 20*mult)  // EXACT or EXACT_GREEDY 
				Polisher.polishIterations = 20*mult;			
		}
		if (Polisher.polishIterations_exact == -1 && Polisher.polishMethod != PolishType.exhaust_twocut){
			if (Polisher.polishAligmentMethod == AlignmentType.profile ) {
				Polisher.polishIterations_exact = 0;
			} else {
				Polisher.polishIterations_exact = Math.round(Polisher.polishIterations/5);
			}
		}	

		if (consistency && PairwiseAlignmentsContainer.neighborCount > K-2) {
			PairwiseAlignmentsContainer.neighborCount = K-2;
		}
		
		int maxLen = 0;
		for (int i=0; i<seqs.length; i++) {
			if (seqs[i].length > maxLen) maxLen = seqs[i].length; 
		}
		
		
	}
	
	public void initialize (Configuration c, Inputs i) { 
		in = i;
		conf = c;

		SequenceFileReader seqReader = new SequenceFileReader(in.fileA, !initAlignmentProvided);		
		chars = seqReader.getSeqs();
		
		this.file = in.fileA;

		if (Tree.distanceType == null ){ 
			if (chars.length < kmerTreeCutoff ) {
				Tree.distanceType = DistanceType.normcost;
				Tree.iterations = 1;
			} else {
				Tree.distanceType = DistanceType.kmer_normcost;
				if (Tree.iterations < 0) {
					Tree.iterations = 2;
				}
			}
		}

		
		initialize (conf.sc.convertSeqsToInts(chars), seqReader.getNames(), conf, in);
	}
	
	protected Alignment[] getAlignments (int[][] seqs) {
		Alignment[] alignments;
		
		if (initAlignmentProvided) {
			alignments = new Alignment[1];
			int[] ids = new int[seqs.length];
			for (int i=0; i<seqs.length; i++) ids[i] = i;
			alignments[0] = Alignment.buildNewAlignment(seqs, ids, conf, in);
		} else {
			alignments = new Alignment[seqs.length];
			for (int i=0; i<seqs.length; i++) {
				alignments[i] = Alignment.buildNewAlignment(seqs[i], i, conf, in );
			}
		}
		return alignments;
	}
	
	protected void getAligner() {		

		if (Aligner.alignmentMethod == Aligner.AlignmentType.profile) {
			al = new ProfileAligner();
		} else if (Aligner.alignmentMethod == AlignmentType.exact || 
				Aligner.alignmentMethod == AlignmentType.mixed ) {
//			if (Aligner.useStructure)
//				al = new StructureExactCountAligner_Time();		
//			else
				al = new ExactCountAligner_Time();		
		}
		if (consistency) {
			al = new ConsistencyAligner(al);
		}
		al.config = conf;
	}

	protected void getPolishAligner() {		
	//	if (Aligner.alignmentMethod == AlignmentType.profile || 
	//				Polisher.polishAligmentMethod == AlignmentType.profile) {
		al = new ProfileAligner();
	//	} else {
		//} else if (Aligner.alignmentMethod == Aligner.AlignmentType.exact ||
		//			Aligner.alignmentMethod == Aligner.AlignmentType.mixed	) {
	//		al = new ExactCountAligner_Time();		
	//	}
		al.config = conf;
	}
		
	protected Tree getTree (Alignment[] als, int currIteration) {
		if (Tree.treeType == Tree.TreeType.mst)
			return new MST_Tree(als, al, dist, currIteration, verbosity, conf);
		else if (Tree.treeType == Tree.TreeType.entered) {
			EnteredTree.names = names;
			return new EnteredTree(als, al, dist, currIteration, verbosity, conf);
		} else {
			LogWriter.stdErrLogln("unrecognized tree type");
			throw new GenericOpalException("unrecognized tree type");
		}
	}
	
	protected void getDistance(int currIteration, Configuration config) {
		Aligner pairAligner;
		pairAligner = new PairAligner(al);	
		pairAligner.config = conf;

		
//		if (Tree.distanceType == DistanceType.normcost || Tree.distanceType == DistanceType.kmer_normcost) 
//		if (Tree.distanceType == DistanceType.normcost) 
		if (Tree.distanceType == DistanceType.normcost || 
				(currIteration<Tree.iterations && Tree.distanceType == DistanceType.kmer_normcost)) 		
			dist = new NormCostDistance(pairAligner);
		else if (Tree.distanceType == DistanceType.kmer || 
				(currIteration==Tree.iterations && Tree.distanceType == DistanceType.kmer_normcost)) 
			dist = new KmerDistance(pairAligner);
		else if (Tree.distanceType == DistanceType.pctid)
			dist = new PctIDDistance(pairAligner);
		else {
			LogWriter.stdErrLogln("unrecognized choice of distance!");
			throw new GenericOpalException("unrecognized choice of distance!");
		}
	}
	
}

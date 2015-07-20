package opal;

import java.util.Date;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import opal.IO.CostMatrix;
import opal.IO.OpalLogWriter;

import com.traviswheeler.libs.*;

import opal.makers.*;
import opal.polish.Polisher;
import opal.tree.Tree;
import opal.align.*;
import opal.exceptions.GenericOpalException;
import opal.IO.*;

public class Opal {


	public static void main(String[] argv) {
		OpalLogWriter.setLogger(new DefaultLogger());
		
		ArgumentHandler argHandler = new ArgumentHandler(argv);
				
		if (argHandler.getVerbosity() > 0) {
			OpalLogWriter.printVersion();
			OpalLogWriter.stdErrLogln("");
		}
		
		String costName = argHandler.getCostName(); 
		String fileA = argHandler.getFileA();
		String fileB = argHandler.getFileB();
		String structFileA = argHandler.getStructFileA();
		String structFileB = argHandler.getStructFileB();
		int gamma = argHandler.getGamma();
		int lambda = argHandler.getLambda();
		int gammaTerm = argHandler.getGammaTerm();
		int lambdaTerm = argHandler.getLambdaTerm();
			
		int verbosity = argHandler.getVerbosity();
		boolean toUpper = argHandler.isToUpper(); 
		boolean justDoConvert= argHandler.isJustDoConvert();
		boolean justDoSubOpt= argHandler.isJustDoSubOpt();
		boolean justTree= argHandler.isJustTree();
		
		
				
		Date start = new Date();
 
		

		try {		
			if (Aligner.useStructure) {
				CostMatrix.isDNA = false;
				costName = "BLOSUM62"; 
				CostMatrix.initialize(costName);
			} else {
				if (costName.equals("")) { 
					//need to figure out if it's DNA or protein
					costName = getCostName(fileA);
				}
				
				if (costName.equals("DNA")) {
					CostMatrix.isDNA = true;
					if (gamma == -1 ) {
						gamma = CostMatrix.dnaDefaultGamma;
					}
					if (lambda == -1 ) {
						lambda = CostMatrix.dnaDefaultLambda;   
					}
				} else { //blosum62
					CostMatrix.isDNA = false;
					if (gamma == -1 ) {
						gamma = CostMatrix.protDefaultGamma;
						gammaTerm = CostMatrix.protDefaultGammaTerm;
					}
					if (lambda == -1 ) {
						lambda = CostMatrix.protDefaultLambda;
						lambdaTerm = CostMatrix.protDefaultLambdaTerm;   
					}
				}		
	
				if (gammaTerm == -1) gammaTerm = gamma;
				if (lambdaTerm == -1) lambdaTerm = lambda;
				
				CostMatrix.initialize(costName, gamma, gammaTerm, lambda, lambdaTerm);
			}
			
			
			
			Alignment.setAlphabetLength(CostMatrix.getChars().length);
			SequenceConverter sc = new SequenceConverter(CostMatrix.getChars());
			Aligner.setSequenceConverter(sc);		
	
			Aligner.setParams();
	
			if (AlignmentMaker_SingleSequences.consistency) {
				int increase = AlignmentMaker_SingleSequences.consistency_subs_increase;
				CostMatrix.increaseCosts(increase); // arbitrary number ... just trying to make w-w identities non-zero cost, and high enough to avoid rounding effects
				CostMatrix.multiplyCosts(2);
				
				Aligner.increaseGapCosts(increase/2);
				Aligner.multiplyGapCosts(2);
			}
			
			AlignmentMaker am;
			
			if (Tree.treeType == Tree.TreeType.entered) {
				Tree.iterations = 1;
			}

		
			if (justDoSubOpt) {
				am = new AlignmentMaker_SuboptimalityTester();
				am.initialize(fileA, structFileA);
				am.buildAlignment( costName, verbosity, toUpper);
			} else if (justDoConvert) {
				am = new AlignmentMaker_Converter();
				am.initialize(fileA, structFileA);
				am.buildAlignment(costName, verbosity, toUpper);
			} else if (fileB != null) { // alignalign call
				am = new AlignmentMaker_TwoAlignments();
				am.initialize(fileA, fileB);
				am.buildAlignment (costName, verbosity, toUpper);
			} else if (justTree || Tree.justPWDists) { // 
				am = new AlignmentMaker_SingleSequences();
				am.initialize(fileA, structFileA);
				((AlignmentMaker_SingleSequences)am).buildTree(costName, verbosity); 
			} else { // multiple alignment call
				am = new AlignmentMaker_SingleSequences();
				am.initialize(fileA, structFileA);
				//final Timer queueRunner = new Timer(); 
				//ShutdownHandler sh = new ShutdownHandler(queueRunner, (AlignmentMaker_SingleSequences)am);
				//Runtime.getRuntime().addShutdownHook(sh);
				//Thread.setDefaultUncaughtExceptionHandler(sh);  // this should be available for all calls ... but this is the one I use a lot, so it's all I've implemented it for
	
				am.buildAlignment( costName, verbosity, toUpper); 
			}
		} catch (GenericOpalException e ) {
			System.exit(1);
		}
		
		if (verbosity>0) {
			Date now = new Date();
			long diff = now.getTime() - start.getTime();
			System.err.printf("Total time for job: %.1f seconds\n\n", ((float)diff/1000 + .05));
		}
		
		System.exit(0);
	}

	private static String getCostName (String file) {
		SequenceFileReader seqReader = new SequenceFileReader(file, true);
		char[][] seqs = seqReader.getSeqs();

		//this is a list of characters that are amino acids, but not in the DNA
		// ambiguity code list
		Pattern forceAA_pattern = Pattern.compile("[QEILFPZ]",Pattern.CASE_INSENSITIVE);	
        Matcher matcher;       
		for (int i=0; i<seqs.length; i++) {
			String str = String.valueOf(seqs[i]);
			matcher = forceAA_pattern.matcher(str);
			if (matcher.find() ) {
				return "BLOSUM62"; 			   
			}
			//if I wanted to be really anal, I'd account for the exceedingly 
			//rare case where the sequences are proteins, but contain none of those
			//protein-only characters, and maybe count the number of ACGTs, but
			//that could still be wrong.
		}
		return "DNA";
	}

		
}
	

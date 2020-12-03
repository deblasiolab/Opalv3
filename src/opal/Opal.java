package opal;

import java.util.Date;

import opal.IO.Configuration;
import opal.IO.OpalLogWriter;

import com.traviswheeler.libs.*;

import opal.makers.*;
import opal.tree.Tree;
import opal.exceptions.GenericOpalException;
import opal.IO.*;
import facet.FacetAlignment;
import facet.Facet;
import opal.realignment.realignmentDriver;

import java.util.concurrent.Executors;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
class runAlignment extends Thread{
	Configuration conf;
	Inputs in;
	AlignmentMaker am;
	int[][] alignmentInstance;
	int[][] preRealignmentAlignmentInstance;
	double facetScore =-1;
	double preRealignmentFacetScore =-1;
	Configuration[] realignmentConfigList = null;
	
	
	runAlignment(Configuration c, Inputs i){
		conf = c;
		in = new Inputs(i);
	}
	runAlignment(Configuration c, Inputs i, Configuration[] rcList){
		this(c,i);
		realignmentConfigList = rcList;
	}
	
	public void run(){
		if(in.verbosity>1) System.err.println("Config : " + conf);
		
		try {		
			
			
			if (AlignmentMaker_SingleSequences.consistency) {
				int increase = AlignmentMaker_SingleSequences.consistency_subs_increase;
				conf.cost.increaseCosts(increase); // arbitrary number ... just trying to make w-w identities non-zero cost, and high enough to avoid rounding effects
				conf.cost.multiplyCosts(2);
				
				conf.increaseGapCosts(increase/2);
				conf.multiplyGapCosts(2);
			}
			
			
			
			if (Tree.treeType == Tree.TreeType.entered) {
				Tree.iterations = 1;
			}

			FacetAlignment fa = null;
			if (in.justDoSubOpt) {
				am = new AlignmentMaker_SuboptimalityTester();
				//am.initialize(in.fileA, in.structFileA);
				am.initialize(conf, in);
				alignmentInstance = am.buildAlignment();
				if(in.structFileA != null){
					fa = new FacetAlignment(conf.sc.convertIntsToSeqs(alignmentInstance),in.structure.structure);
				}
			} else if (in.justDoConvert) {
				am = new AlignmentMaker_Converter();
				//am.initialize(in.fileA, in.structFileA);
				am.initialize(conf, in);
				alignmentInstance = am.buildAlignment();
				if(in.structFileA != null){
					fa = new FacetAlignment(conf.sc.convertIntsToSeqs(alignmentInstance),in.structure.structure);
				}
			} else if (in.fileB != null) { // alignalign call
				am = new AlignmentMaker_TwoAlignments();
				//am.initialize(in.fileA, in.fileB);
				//Aligner.useStructure = false;
				am.initialize(conf, in);
				alignmentInstance = am.buildAlignment ();
				if(in.structFileA != null){
					fa = new FacetAlignment(((AlignmentMaker_TwoAlignments)am).result,in.structure.structure);
				}
				
			} else if (in.justTree || Tree.justPWDists) { // 
				am = new AlignmentMaker_SingleSequences();
				//am.initialize(in.fileA, in.structFileA);
				am.initialize(conf, in);
				((AlignmentMaker_SingleSequences)am).buildTree(conf, in);
			} else { // multiple alignment call
				am = new AlignmentMaker_SingleSequences();
				//am.initialize(in.fileA, in.structFileA);
				am.initialize(conf, in);
				//final Timer queueRunner = new Timer(); 
				//ShutdownHandler sh = new ShutdownHandler(queueRunner, (AlignmentMaker_SingleSequences)am);
				//Runtime.getRuntime().addShutdownHook(sh);
				//Thread.setDefaultUncaughtExceptionHandler(sh);  // this should be available for all calls ... but this is the one I use a lot, so it's all I've implemented it for

				alignmentInstance = am.buildAlignment();
				
				
				if(in.structFileA != null){
					
					if(realignmentConfigList != null){
						Configuration newRealignmentConfigList[] = new Configuration[realignmentConfigList.length];
						for(int rNum = 0; rNum < realignmentConfigList.length; rNum++){
							newRealignmentConfigList[rNum] = new Configuration(realignmentConfigList[rNum]);
						}
						preRealignmentAlignmentInstance = alignmentInstance.clone();
						preRealignmentFacetScore = Facet.defaultValue(
								new FacetAlignment(conf.sc.convertIntsToSeqs(preRealignmentAlignmentInstance),in.structure.structure),
								conf.useLegacyFacetFunction
						);
						if(conf.useTCSforAdvising)
							preRealignmentFacetScore = TCS.TCSValue(preRealignmentAlignmentInstance, am, conf);
						realignmentDriver realigner = new realignmentDriver(conf.sc.convertIntsToSeqs(preRealignmentAlignmentInstance),in.structure.structure, newRealignmentConfigList, conf, (float) preRealignmentFacetScore, am);
						if(conf.realignment_window_type == Configuration.WINDOW_SIZE.VALUE) realigner.simpleRealignment((int)conf.realignment_window_value);
						else if(conf.realignment_window_type == Configuration.WINDOW_SIZE.PERCENTAGE){
							int window_size = (int)(conf.realignment_window_value * (float)preRealignmentAlignmentInstance[0].length);
							if(window_size>conf.maximum_realignment_window_value) window_size = (int)conf.maximum_realignment_window_value;
							if(window_size<conf.minimum_realignment_window_value) window_size = (int)conf.minimum_realignment_window_value;
							realigner.simpleRealignment(window_size);
						}
						alignmentInstance = realigner.newAlignment();
						newRealignmentConfigList = null;
						
					}
					fa = new FacetAlignment(conf.sc.convertIntsToSeqs(alignmentInstance),in.structure.structure);
					
				}
			}
			
			if(fa != null){
				if(in.featureOutputFile!=null){
					String fname = in.featureOutputFile.replace("__CONFIG__", conf.toString());
      fname = fname.replace("__MATRIX__",conf.cost.costName);
      fname = fname.replace("__LAMBDA__",String.valueOf(conf.getLambda()));
      fname = fname.replace("__LAMBDATERM__",String.valueOf(conf.getLambdaTerm()));
      fname = fname.replace("__GAMMA__",String.valueOf(conf.getGamma()));
      fname = fname.replace("__GAMMATERM__",String.valueOf(conf.getGammaTerm()));
					if(conf.repetition>=0) fname = fname.replace("__ITTERATION__", Integer.toString(conf.repetition));
					Facet.outputDefaultFeatures(fname, fa);
				}
				facetScore = Facet.defaultValue(fa,conf.useLegacyFacetFunction);
				if(conf.useTCSforAdvising)
					facetScore = TCS.TCSValue(alignmentInstance, am, conf);
			}
			//in.structFileA = "";
			//System.err.println("File: " + in.structFileA);
			
			
			//System.err.printf("facet value --  %.6f\n\n",facetScore);
		} catch (GenericOpalException e ) {
			System.exit(1);
		} /*catch (FileNotFoundException e){
			System.err.println("File Error:" + e.toString());
			System.exit(15);
		}*/
    
		/* Moved here to allow for print as you go on long runs */
		if(in.configOutputFile != null) print();
    if(in.configOutputFile != null && realignmentConfigList != null) printPreRealignment();
	}// end run
	
	public boolean print(){
		//if(facetScore==-1){
		//	am = null;
		//}
		if (in.verbosity>0 && facetScore>-1) {
			System.err.printf("\nfacet score: %.6f",facetScore);
		}
		
		if(in.configOutputFile != null){
			String fname = in.configOutputFile.replace("__CONFIG__", conf.toString());
      fname = fname.replace("__MATRIX__",conf.cost.costName);
      fname = fname.replace("__LAMBDA__",String.valueOf(conf.getLambda()));
      fname = fname.replace("__LAMBDATERM__",String.valueOf(conf.getLambdaTerm()));
      fname = fname.replace("__GAMMA__",String.valueOf(conf.getGamma()));
      fname = fname.replace("__GAMMATERM__",String.valueOf(conf.getGammaTerm()));
			if(conf.repetition>=0) fname = fname.replace("__ITTERATION__", Integer.toString(conf.repetition));
			if(facetScore>=0) fname = fname.replace("__FACETSCORE__", "facetScore" + Double.toString(facetScore));
			return am.printOutput(alignmentInstance, fname);
			
		}
		else return am.printOutput(alignmentInstance, null);
		
	}
	
	public boolean printBest(){
		if (facetScore != -1 && in.verbosity>-1) {
			System.err.printf("\nbest facet value --  %.6f (%s)",facetScore, conf );
		}

		String fname = null; 
    if(in.bestOutputFile != null){
      fname = in.bestOutputFile.replace("__CONFIG__", conf.toString());
      fname = fname.replace("__MATRIX__",conf.cost.costName);
      fname = fname.replace("__LAMBDA__",String.valueOf(conf.getLambda()));
      fname = fname.replace("__LAMBDATERM__",String.valueOf(conf.getLambdaTerm()));
      fname = fname.replace("__GAMMA__",String.valueOf(conf.getGamma()));
      fname = fname.replace("__GAMMATERM__",String.valueOf(conf.getGammaTerm()));
		  if(conf.repetition>=0) fname = fname.replace("__ITTERATION__", Integer.toString(conf.repetition));
		  if(facetScore>=0) fname = fname.replace("__FACETSCORE__", "facetScore" + Double.toString(facetScore));
    }

		return am.printOutput(alignmentInstance, fname);
	}
	
	public boolean printBestIncludePreRealignment(){
		if(facetScore>preRealignmentFacetScore){
			if (in.verbosity>-1) {
				System.err.printf("\nbest facet value --  %.6f (%s)",facetScore, conf );
			}
			String fname = null;
      if(in.bestOutputFileIncludePreRealignment != null){
        fname = in.bestOutputFileIncludePreRealignment.replace("__CONFIG__", conf.toString());
      fname = fname.replace("__MATRIX__",conf.cost.costName);
      fname = fname.replace("__LAMBDA__",String.valueOf(conf.getLambda()));
      fname = fname.replace("__LAMBDATERM__",String.valueOf(conf.getLambdaTerm()));
      fname = fname.replace("__GAMMA__",String.valueOf(conf.getGamma()));
      fname = fname.replace("__GAMMATERM__",String.valueOf(conf.getGammaTerm()));
			  if(conf.repetition>=0) fname = fname.replace("__ITTERATION__", Integer.toString(conf.repetition));
			  if(facetScore>=0) fname = fname.replace("__FACETSCORE__", "facetScore" + Double.toString(facetScore));
      }
      return am.printOutput(alignmentInstance, fname);
		}else{
			if (in.verbosity>0 && facetScore>-1) {
				System.err.printf("\npre-realignment facet score: %.6f",preRealignmentFacetScore);
			}
			String fname = null;
      if(in.bestOutputFileIncludePreRealignment != null){
        fname = in.bestOutputFileIncludePreRealignment.replace("__CONFIG__", conf.toString());
      fname = fname.replace("__MATRIX__",conf.cost.costName);
      fname = fname.replace("__LAMBDA__",String.valueOf(conf.getLambda()));
      fname = fname.replace("__LAMBDATERM__",String.valueOf(conf.getLambdaTerm()));
      fname = fname.replace("__GAMMA__",String.valueOf(conf.getGamma()));
      fname = fname.replace("__GAMMATERM__",String.valueOf(conf.getGammaTerm()));
			  if(conf.repetition>=0) fname = fname.replace("__ITTERATION__", Integer.toString(conf.repetition));
			  if(facetScore>=0) fname = fname.replace("__FACETSCORE__", "facetScore" + Double.toString(facetScore));
      }
			return am.printOutput(preRealignmentAlignmentInstance, fname, false);
			
		}
	}
	
	public boolean printPreRealignment(){
		if (in.verbosity>0 && facetScore>-1) {
			System.err.printf("\npre-realignment facet score: %.6f",preRealignmentFacetScore);
		}
		if(in.preRealignmentOutputFile != null){
			String fname = null; 
      if(in.preRealignmentOutputFile != null){
        in.preRealignmentOutputFile.replace("__CONFIG__", conf.toString());
      fname = fname.replace("__MATRIX__",conf.cost.costName);
      fname = fname.replace("__LAMBDA__",String.valueOf(conf.getLambda()));
      fname = fname.replace("__LAMBDATERM__",String.valueOf(conf.getLambdaTerm()));
      fname = fname.replace("__GAMMA__",String.valueOf(conf.getGamma()));
      fname = fname.replace("__GAMMATERM__",String.valueOf(conf.getGammaTerm()));
			  if(conf.repetition>=0) fname = fname.replace("__ITTERATION__", Integer.toString(conf.repetition));
			  if(facetScore>=0) fname = fname.replace("__FACETSCORE__", "facetScore" + Double.toString(preRealignmentFacetScore));
      }
			return am.printOutput(preRealignmentAlignmentInstance, fname, false);
		}else{
			return am.printOutput(preRealignmentAlignmentInstance, null, false);
		}
		
	}
	
	public boolean printBestPreRealignment(){
		if (in.verbosity>-1) {
			System.err.printf("\nbest pre-realignment facet value --  %.6f (%s)",preRealignmentFacetScore, conf );
		}
		String fname = null;
    if(in.bestPreRealignmentOutputFile != null){
      in.bestPreRealignmentOutputFile.replace("__CONFIG__", conf.toString());
      fname = fname.replace("__MATRIX__",conf.cost.costName);
      fname = fname.replace("__LAMBDA__",String.valueOf(conf.getLambda()));
      fname = fname.replace("__LAMBDATERM__",String.valueOf(conf.getLambdaTerm()));
      fname = fname.replace("__GAMMA__",String.valueOf(conf.getGamma()));
      fname = fname.replace("__GAMMATERM__",String.valueOf(conf.getGammaTerm()));
		  if(conf.repetition>=0) fname = fname.replace("__ITTERATION__", Integer.toString(conf.repetition));
		  if(facetScore>=0) fname = fname.replace("__FACETSCORE__", "facetScore" + Double.toString(preRealignmentFacetScore));
    }
    return am.printOutput(preRealignmentAlignmentInstance, fname, false);
	}	
	
	public boolean printBestPreRealignmentsRealignment(){
		if (in.verbosity>-1) {
			System.err.printf("\nbest pre-realignments realignment facet value --  %.6f (%s)",facetScore, conf );
		}
		String fname = null;
    if(in.bestPreRealignmentsRealignmentOutputFile != null){
      in.bestPreRealignmentsRealignmentOutputFile.replace("__CONFIG__", conf.toString());
      fname = fname.replace("__MATRIX__",conf.cost.costName);
      fname = fname.replace("__LAMBDA__",String.valueOf(conf.getLambda()));
      fname = fname.replace("__LAMBDATERM__",String.valueOf(conf.getLambdaTerm()));
      fname = fname.replace("__GAMMA__",String.valueOf(conf.getGamma()));
      fname = fname.replace("__GAMMATERM__",String.valueOf(conf.getGammaTerm()));
		  if(conf.repetition>=0) fname = fname.replace("__ITTERATION__", Integer.toString(conf.repetition));
		  if(facetScore>=0) fname = fname.replace("__FACETSCORE__", "facetScore" + Double.toString(preRealignmentFacetScore));	
    }
    return am.printOutput(alignmentInstance, fname);
	}	
	
	public boolean printBestPreRealignmentsRealignmentIncludePreRealignment(){
		if(facetScore>preRealignmentFacetScore){
			if (in.verbosity>-1) {
				System.err.printf("\nbest pre-realignments realignment facet value --  %.6f (%s)",facetScore, conf );
			}
			String fname = null;
      if(in.bestPreRealignmentsRealignmentOutputFileIncludePreRealignment != null){
        fname = in.bestPreRealignmentsRealignmentOutputFileIncludePreRealignment.replace("__CONFIG__", conf.toString());
      fname = fname.replace("__MATRIX__",conf.cost.costName);
      fname = fname.replace("__LAMBDA__",String.valueOf(conf.getLambda()));
      fname = fname.replace("__LAMBDATERM__",String.valueOf(conf.getLambdaTerm()));
      fname = fname.replace("__GAMMA__",String.valueOf(conf.getGamma()));
      fname = fname.replace("__GAMMATERM__",String.valueOf(conf.getGammaTerm()));
			  if(conf.repetition>=0) fname = fname.replace("__ITTERATION__", Integer.toString(conf.repetition));
			  if(facetScore>=0) fname = fname.replace("__FACETSCORE__", "facetScore" + Double.toString(preRealignmentFacetScore));	
      }
			return am.printOutput(alignmentInstance, fname);
		}else{
			if (in.verbosity>-1) {
				System.err.printf("\nbest pre-realignment facet value --  %.6f (%s)",preRealignmentFacetScore, conf );
			}
			String fname = null;
      if(in.bestPreRealignmentsRealignmentOutputFileIncludePreRealignment != null){
        fname = in.bestPreRealignmentsRealignmentOutputFileIncludePreRealignment.replace("__CONFIG__", conf.toString());
      	fname = fname.replace("__MATRIX__",conf.cost.costName);
     	 	fname = fname.replace("__LAMBDA__",String.valueOf(conf.getLambda()));
      	fname = fname.replace("__LAMBDATERM__",String.valueOf(conf.getLambdaTerm()));
      	fname = fname.replace("__GAMMA__",String.valueOf(conf.getGamma()));
      	fname = fname.replace("__GAMMATERM__",String.valueOf(conf.getGammaTerm()));
			  if(conf.repetition>=0) fname = fname.replace("__ITTERATION__", Integer.toString(conf.repetition));
			  if(facetScore>=0) fname = fname.replace("__FACETSCORE__", "facetScore" + Double.toString(preRealignmentFacetScore));		
      }
			return am.printOutput(preRealignmentAlignmentInstance, fname, false);
		}
	}
}

public class Opal {


	public static void main(String[] argv) {

		OpalLogWriter.setLogger(new DefaultLogger());

		ArgumentHandler argHandler = new ArgumentHandler(argv);
				
		if (argHandler.getVerbosity() > 0) {
			OpalLogWriter.printVersion();
			OpalLogWriter.stdErrLogln("");
		}
		
		/*String fileA = argHandler.getFileA();
		String fileB = argHandler.getFileB();
		String structFileA = argHandler.getStructFileA();
		String structFileB = argHandler.getStructFileB();
		
		
		
		String costName = argHandler.getCostName(); 
		int gamma = argHandler.getGamma();
		int lambda = argHandler.getLambda();
		int gammaTerm = argHandler.getGammaTerm();
		int lambdaTerm = argHandler.getLambdaTerm();
			
		int verbosity = argHandler.getVerbosity();
		boolean toUpper = argHandler.isToUpper(); 
		boolean justDoConvert= argHandler.isJustDoConvert();
		boolean justDoSubOpt= argHandler.isJustDoSubOpt();
		boolean justTree= argHandler.isJustTree();*/
		
		
				
		Date start = new Date();
 
		
		Configuration[] advising_config = argHandler.getAdvisingConfigs();
		Configuration[] realignment_config = argHandler.getRealignmentConfigs();
		Inputs input = argHandler.getInputs();
		runAlignment[] thread = new runAlignment[advising_config.length];
		
		int last_joined = -1;
		int max_threads = Runtime.getRuntime().availableProcessors();
		if(argHandler.getMaxThreads() > -1) max_threads = argHandler.getMaxThreads();
		if(input.verbosity>0){
			OpalLogWriter.stdErrLogln("The number of available threads is " + Runtime.getRuntime().availableProcessors() + ", Opal will use " + max_threads);
		}
    ThreadPoolExecutor executor = (ThreadPoolExecutor) Executors.newFixedThreadPool(max_threads);
		//ThreadPoolExecutor executor = (ThreadPoolExecutor) Executors.newWorkStealingPool(max_threads);
		int maxIndex = 0;
    int maxPreRealignmentIndex = 0;
		for(int i=0;i<advising_config.length;i++){
			thread[i] = new runAlignment(advising_config[i],input, realignment_config);
			//thread[i] = new printLine(config[i],i);
			//thread[i].start();
      executor.execute(thread[i]);
    }
    executor.shutdown();
    try{
      executor.awaitTermination(10,TimeUnit.HOURS);
    }catch(InterruptedException e){
      System.out.println("Command interupted: " + e);
    }
		/*  	if(i-last_joined>=max_threads){
				last_joined++;
				try{
					thread[last_joined].join();
					if(input.configOutputFile != null) thread[last_joined].print();
					if(input.configOutputFile != null && realignment_config != null) thread[last_joined].printPreRealignment(); 
					if(thread[last_joined].facetScore > thread[maxIndex].facetScore){
						if(maxPreRealignmentIndex != maxIndex) thread[maxIndex] = null;
						maxIndex = last_joined;
					} 
					if(thread[last_joined].preRealignmentFacetScore > thread[maxPreRealignmentIndex].preRealignmentFacetScore){
						if(maxPreRealignmentIndex != maxIndex) thread[maxPreRealignmentIndex] = null;
						maxPreRealignmentIndex = last_joined;
					}
					if(last_joined != maxIndex && last_joined != maxPreRealignmentIndex) thread[last_joined] = null;
				}catch(InterruptedException e){
					OpalLogWriter.stdErrLogln("InterruptedException "+e.toString());
					throw new GenericOpalException("InterruptedException "+e.toString());
				}
			//}
		}*/
		
		for(last_joined++;last_joined<advising_config.length;last_joined++){
			try{
				thread[last_joined].join();
				if(thread[last_joined].facetScore > thread[maxIndex].facetScore){
					if(maxPreRealignmentIndex != maxIndex) thread[maxIndex] = null;
					maxIndex = last_joined;
				} 
				if(thread[last_joined].preRealignmentFacetScore > thread[maxPreRealignmentIndex].preRealignmentFacetScore){
					if(maxPreRealignmentIndex != maxIndex) thread[maxPreRealignmentIndex] = null;
					maxPreRealignmentIndex = last_joined;
				}
				if(last_joined != maxIndex && last_joined != maxPreRealignmentIndex) thread[last_joined] = null;
			}catch(InterruptedException e){
				OpalLogWriter.stdErrLogln("InterruptedException "+e.toString());
				throw new GenericOpalException("InterruptedException "+e.toString());
			}
		}
		
		
		
		//if(advising_config.length>1 && thread[maxIndex]!=null && thread[maxIndex].facetScore>=0)
		if(!thread[maxIndex].printBest()) 
			System.err.println("Print returned false");
		
		if(advising_config.length>1 && realignment_config != null  && thread[maxPreRealignmentIndex]!=null && thread[maxPreRealignmentIndex].facetScore>=0){
				if(input.bestPreRealignmentOutputFile != null) 
					if(!thread[maxPreRealignmentIndex].printBestPreRealignment()) 
						System.err.println("Print returned false");
				if(input.bestPreRealignmentsRealignmentOutputFile != null)
					if(!thread[maxPreRealignmentIndex].printBestPreRealignmentsRealignment()) 
						System.err.println("Print returned false");
				if(input.bestPreRealignmentsRealignmentOutputFileIncludePreRealignment != null)
					if(!thread[maxPreRealignmentIndex].printBestPreRealignmentsRealignmentIncludePreRealignment()) 
						System.err.println("Print returned false");
		}
		
		if(advising_config.length>1 && realignment_config != null  && thread[maxIndex]!=null && thread[maxIndex].facetScore>=0 
				&& thread[maxPreRealignmentIndex]!=null && thread[maxPreRealignmentIndex].facetScore>=0 
				&& thread[maxPreRealignmentIndex].preRealignmentFacetScore < thread[maxIndex].facetScore){
			if(input.bestOutputFileIncludePreRealignment != null)
				if(!thread[maxIndex].printBestIncludePreRealignment()) 
					System.err.println("Print returned false");
		}else if(advising_config.length>1 && thread[maxPreRealignmentIndex]!=null && thread[maxPreRealignmentIndex].facetScore>=0){
			if(input.bestOutputFileIncludePreRealignment != null)
				if(!thread[maxPreRealignmentIndex].printBestIncludePreRealignment()) 
					System.err.println("Print returned false");
		}
		
		if (input.verbosity>0) {
			Date now = new Date();
			long diff = now.getTime() - start.getTime();
			System.err.printf("Total time for job: %.1f seconds\n\n", ((float)diff/1000 + .05));
		}
		
		System.exit(0);
	}



		
}
	

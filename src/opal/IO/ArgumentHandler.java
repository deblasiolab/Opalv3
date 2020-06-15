package opal.IO;

import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.FileNotFoundException;

import opal.IO.AlignmentWriter.OutputType;
import opal.IO.Configuration;
import opal.align.Aligner;
import opal.align.ConsistencyAligner;
import opal.align.PairSuboptimalityMatrices;
import opal.align.PairwiseAlignmentsContainer;
import opal.align.ProfileAligner;
import opal.align.StructureAlignment;
import opal.align.Aligner.AlignmentType;
import opal.exceptions.GenericOpalException;
import opal.makers.AlignmentMaker;
import opal.makers.AlignmentMaker_SingleSequences;
import opal.makers.AlignmentMaker.OutputOrderType;
import opal.polish.Polisher;
import opal.tree.EnteredTree;
import opal.tree.Tree;
import gnu.getopt.Getopt;
import gnu.getopt.LongOpt;

import com.traviswheeler.libs.LogWriter;

public class ArgumentHandler {

	String costName = ""; 
	String fileA = null;
	String fileB = null; // will remain null for a multiple alignment request.  If a value is assigned, it's a request to align two alignment files
	String structFileA = null;
	String structFileB = null;
	int gamma = -1;
	int lambda = -1;
	int gammaTerm = -1;
	int lambdaTerm = -1;
	StructureAlignment.ParamModel structmodel = null;
	Configuration[] advising_configs;
	Configuration[] realignment_configs;
	int repeat_config = 1;
	int max_threads = -1;
	boolean useLegacyFacetFunction = false;
	boolean doReverse = false;
	
	public String configOutputFile = null;
	public String bestOutputFile = null;
	public String bestOutputFileIncludePreRealignment = null;
	public String featureOutputFile = null;
	public String preRealignmentOutputFile = null;
	public String bestPreRealignmentOutputFile = null;
	public String bestPreRealignmentsRealignmentOutputFile = null;
	public String bestPreRealignmentsRealignmentOutputFileIncludePreRealignment = null;
		
	int verbosity = 1;
	boolean toUpper = false; 
	int c;
	String arg;
	boolean justDoConvert= false;
	boolean justDoSubOpt= false;
	boolean justTree= false;

	Configuration.THRESHOLD_TYPE temp_threshold_type = null;
	Configuration.WINDOW_SIZE_MINIMUM temp_window_minimum_type = null;
	Configuration.WINDOW_SIZE temp_window_type = null;
	Configuration.REALIGNMENT_TERMINALS temp_realignment_terminals = null;
	float temp_realign_threshold = -1;
	float temp_realign_threshold_lower = (float)-10000.0;
	float temp_realign_window_size = -1;
	float temp_maximum_realign_window_size = -1;
	float temp_minimum_realign_window_size = -1;
	float temp_realign_minimum_value = -1;
	float temp_realignmentWindowWeightDecay = -1;
	int temp_realignment_ittertations = -1;
	boolean temp_realignment_save_threshold = false;
	
	public String temp_temporaryFileDirectory = null;
	public String temp_tcoffeeDirectory = null;
	public boolean temp_useTCSforRealignment = false;
	public boolean temp_useTCSforAdvising = false;

	public boolean temp_useLeftTerminal = true;
	public boolean temp_useRightTerminal = true;
	
	
	public ArgumentHandler (String argString) {
		this(argString.split("\\s"));
	}
	
	public ArgumentHandler (String[] argv) {
		/*
		LogWriter.stdErrLogln("Got these arguments: ");
		for (String s : argv) {
			LogWriter.stdErrLogln(s);
		}
		*/
		
		LongOpt[] longopts = new LongOpt[93];
		int longopts_index=0;
		longopts[longopts_index++] = new LongOpt("help", LongOpt.NO_ARGUMENT, null, 'h');
		longopts[longopts_index++] = new LongOpt("usage", LongOpt.NO_ARGUMENT, null, 'u');
		
		longopts[longopts_index++] = new LongOpt("cost", LongOpt.REQUIRED_ARGUMENT, null, 'c');
		longopts[longopts_index++] = new LongOpt("in", LongOpt.REQUIRED_ARGUMENT, null, 'a');
		longopts[longopts_index++] = new LongOpt("in2", LongOpt.REQUIRED_ARGUMENT, null, 'b'); // if this is included, it's a request to align two alignments
		longopts[longopts_index++] = new LongOpt("facet_structure", LongOpt.REQUIRED_ARGUMENT, null, 's');
		longopts[longopts_index++] = new LongOpt("init_alignment", LongOpt.REQUIRED_ARGUMENT, null, 'a');
		longopts[longopts_index++] = new LongOpt("distance_type", LongOpt.REQUIRED_ARGUMENT, null, 'd'); // normcost, pctid
		
		longopts[longopts_index++] = new LongOpt("gamma", LongOpt.REQUIRED_ARGUMENT, null, 'g');
		longopts[longopts_index++] = new LongOpt("lambda", LongOpt.REQUIRED_ARGUMENT, null, 'l');
		longopts[longopts_index++] = new LongOpt("gamma_term", LongOpt.REQUIRED_ARGUMENT, null, 'e');
		longopts[longopts_index++] = new LongOpt("lambda_term", LongOpt.REQUIRED_ARGUMENT, null, 'f');
		longopts[longopts_index++] = new LongOpt("configuration_file", LongOpt.REQUIRED_ARGUMENT, null, 'j');

		longopts[longopts_index++] = new LongOpt("realignment_configuration_file", LongOpt.REQUIRED_ARGUMENT, null, 'j');
		longopts[longopts_index++] = new LongOpt("advising_configuration_file", LongOpt.REQUIRED_ARGUMENT, null, 'j');

		longopts[longopts_index++] = new LongOpt("repeat_configurations", LongOpt.REQUIRED_ARGUMENT, null, 'j');
		longopts[longopts_index++] = new LongOpt("max_threads", LongOpt.REQUIRED_ARGUMENT, null, 'j');
		longopts[longopts_index++] = new LongOpt("use_legacy_facet", LongOpt.NO_ARGUMENT, null, 'j');
		
		longopts[longopts_index++] = new LongOpt("quiet", LongOpt.NO_ARGUMENT, null, 'q');
		longopts[longopts_index++] = new LongOpt("silent", LongOpt.NO_ARGUMENT, null, 'q');
		longopts[longopts_index++] = new LongOpt("nomrmal", LongOpt.NO_ARGUMENT, null, 'q');
		longopts[longopts_index++] = new LongOpt("verbose", LongOpt.NO_ARGUMENT, null, 'q');
		longopts[longopts_index++] = new LongOpt("noisy", LongOpt.NO_ARGUMENT, null, 'q');
		
		longopts[longopts_index++] = new LongOpt("align_method", LongOpt.REQUIRED_ARGUMENT, null, 'm'); // exact, profile
		longopts[longopts_index++] = new LongOpt("out", LongOpt.REQUIRED_ARGUMENT, null, 'o'); // fasta, clustalw
		longopts[longopts_index++] = new LongOpt("out_format", LongOpt.REQUIRED_ARGUMENT, null, 'o'); // fasta, clustalw 
		longopts[longopts_index++] = new LongOpt("output_width", LongOpt.REQUIRED_ARGUMENT, null, 'w'); // default depends on output type: 55 for clustal, 80 for fasta
		longopts[longopts_index++] = new LongOpt("upper_case", LongOpt.NO_ARGUMENT, null, 'u'); // default retain input case
		longopts[longopts_index++] = new LongOpt("tree_order", LongOpt.NO_ARGUMENT, null, 'o'); //<<<<<<<<<
		longopts[longopts_index++] = new LongOpt("input_order", LongOpt.NO_ARGUMENT, null, 'o'); //<<<<<<<<<
		longopts[longopts_index++] = new LongOpt("show_cost", LongOpt.NO_ARGUMENT, null, 's'); //<<<<<<<<<
		longopts[longopts_index++] = new LongOpt("out_best", LongOpt.REQUIRED_ARGUMENT, null, 'o');
		longopts[longopts_index++] = new LongOpt("out_best_include_pre", LongOpt.REQUIRED_ARGUMENT, null, 'o'); 
		longopts[longopts_index++] = new LongOpt("out_config", LongOpt.REQUIRED_ARGUMENT, null, 'o'); 
		longopts[longopts_index++] = new LongOpt("out_feature", LongOpt.REQUIRED_ARGUMENT, null, 'o'); 

		longopts[longopts_index++] = new LongOpt("out_prerealignment", LongOpt.REQUIRED_ARGUMENT, null, 'o'); 
		longopts[longopts_index++] = new LongOpt("out_prerealignment_best", LongOpt.REQUIRED_ARGUMENT, null, 'o'); 
		longopts[longopts_index++] = new LongOpt("out_prerealignment_best_realignment", LongOpt.REQUIRED_ARGUMENT, null, 'o'); 
		longopts[longopts_index++] = new LongOpt("out_prerealignment_best_realignment_include_pre", LongOpt.REQUIRED_ARGUMENT, null, 'o');
		
		longopts[longopts_index++] = new LongOpt("linear_cutoff", LongOpt.REQUIRED_ARGUMENT, null, 'z'); // default = 20
		longopts[longopts_index++] = new LongOpt("mixed_alignment_cutoff", LongOpt.REQUIRED_ARGUMENT, null, 'm');
		
		longopts[longopts_index++] = new LongOpt("just_polish", LongOpt.NO_ARGUMENT, null, 'p'); //<<<<<<<<<
		longopts[longopts_index++] = new LongOpt("seed", LongOpt.REQUIRED_ARGUMENT, null, 's'); // random seed for polishing  
		longopts[longopts_index++] = new LongOpt("polish_reps", LongOpt.REQUIRED_ARGUMENT, null, 'r'); // # repetitions of polishing
		longopts[longopts_index++] = new LongOpt("polish", LongOpt.REQUIRED_ARGUMENT, null, 'p'); //<<<<<<<<<
		longopts[longopts_index++] = new LongOpt("convert", LongOpt.NO_ARGUMENT, null, 1);
		longopts[longopts_index++] = new LongOpt("polish_align_method", LongOpt.REQUIRED_ARGUMENT, null, 'p'); //<<<<<<<<<
		longopts[longopts_index++] = new LongOpt("polish_reps_exhaustive", LongOpt.REQUIRED_ARGUMENT, null, 'r'); // # repetitions of exhaustive polishing (only worth doing if first polish is with heuristic) 
		
		longopts[longopts_index++] = new LongOpt("subopt", LongOpt.REQUIRED_ARGUMENT, null, 's');
		longopts[longopts_index++] = new LongOpt("just_subopt", LongOpt.REQUIRED_ARGUMENT, null, 's');
		
		longopts[longopts_index++] = new LongOpt("treein", LongOpt.REQUIRED_ARGUMENT, null, 't');
		longopts[longopts_index++] = new LongOpt("treeout", LongOpt.REQUIRED_ARGUMENT, null, 't');
		longopts[longopts_index++] = new LongOpt("just_tree", LongOpt.REQUIRED_ARGUMENT, null, 't');
		
		longopts[longopts_index++] = new LongOpt("consistency_badscore_mult", LongOpt.REQUIRED_ARGUMENT, null, 'y');
		longopts[longopts_index++] = new LongOpt("consistency_use_avg", LongOpt.NO_ARGUMENT, null, 'y'); //<<<<<<<<<
		longopts[longopts_index++] = new LongOpt("consistency_use_neighbor_weights", LongOpt.NO_ARGUMENT, null, 'y'); //<<<<<<<<<
		longopts[longopts_index++] = new LongOpt("consistency_neighbors", LongOpt.REQUIRED_ARGUMENT, null, 'y');
		longopts[longopts_index++] = new LongOpt("consistency_maxsubtree", LongOpt.REQUIRED_ARGUMENT, null, 'y');
		longopts[longopts_index++] = new LongOpt("consistency_flatten_abc_subopt", LongOpt.REQUIRED_ARGUMENT, null, 'y');
		longopts[longopts_index++] = new LongOpt("consistency_weight", LongOpt.REQUIRED_ARGUMENT, null, 'y');
		longopts[longopts_index++] = new LongOpt("consistency_other_seqs_weight", LongOpt.REQUIRED_ARGUMENT, null, 'y');
		longopts[longopts_index++] = new LongOpt("consistency_blend_type", LongOpt.REQUIRED_ARGUMENT, null, 'y');
		longopts[longopts_index++] = new LongOpt("use_consistency", LongOpt.NO_ARGUMENT, null, 'y'); //<<<<<<<<<
		longopts[longopts_index++] = new LongOpt("consistency_neighbor_dist_thresh", LongOpt.REQUIRED_ARGUMENT, null, 'y');
		longopts[longopts_index++] = new LongOpt("consistency_align_method", LongOpt.REQUIRED_ARGUMENT, null, 'y');
		
		longopts[longopts_index++] = new LongOpt("just_pair_dists", LongOpt.NO_ARGUMENT, null, 'p'); 
		longopts[longopts_index++] = new LongOpt("pess_do_reverse", LongOpt.NO_ARGUMENT, null, 'r'); //<<<<<<<<<
		longopts[longopts_index++] = new LongOpt("tree_iterations", LongOpt.REQUIRED_ARGUMENT, null, 'i'); 
		longopts[longopts_index++] = new LongOpt("structure_file", LongOpt.REQUIRED_ARGUMENT, null, 's');
		longopts[longopts_index++] = new LongOpt("structure_file2", LongOpt.REQUIRED_ARGUMENT, null, 's');
		longopts[longopts_index++] = new LongOpt("structure_model", LongOpt.REQUIRED_ARGUMENT, null, 's');
		
		longopts[longopts_index++] = new LongOpt("dna", LongOpt.NO_ARGUMENT, null, 'd'); //<<<<<<<<<
		longopts[longopts_index++] = new LongOpt("protein", LongOpt.NO_ARGUMENT, null, 'p');
		longopts[longopts_index++] = new LongOpt("dnaAG", LongOpt.REQUIRED_ARGUMENT, null, 'd');
		longopts[longopts_index++] = new LongOpt("dnaCT", LongOpt.REQUIRED_ARGUMENT, null, 'd');
		longopts[longopts_index++] = new LongOpt("dnaCU", LongOpt.REQUIRED_ARGUMENT, null, 'd');

		longopts[longopts_index++] = new LongOpt("realignmentMinimumSizeType", LongOpt.REQUIRED_ARGUMENT, null, 'k');
		longopts[longopts_index++] = new LongOpt("realignmentMinimumSizeValue", LongOpt.REQUIRED_ARGUMENT, null, 'k');
		longopts[longopts_index++] = new LongOpt("realignmentWindowSizeType", LongOpt.REQUIRED_ARGUMENT, null, 'k');
		longopts[longopts_index++] = new LongOpt("realignmentWindowSizeValue", LongOpt.REQUIRED_ARGUMENT, null, 'k');
		longopts[longopts_index++] = new LongOpt("realignmentMinimumWindowSizeValue", LongOpt.REQUIRED_ARGUMENT, null, 'k');
		longopts[longopts_index++] = new LongOpt("realignmentMaximumWindowSizeValue", LongOpt.REQUIRED_ARGUMENT, null, 'k');
		longopts[longopts_index++] = new LongOpt("realignmentThresholdType", LongOpt.REQUIRED_ARGUMENT, null, 'k');
		longopts[longopts_index++] = new LongOpt("realignmentWindowWeightDecay", LongOpt.REQUIRED_ARGUMENT, null, 'k');
		longopts[longopts_index++] = new LongOpt("realignmentThresholdValue", LongOpt.REQUIRED_ARGUMENT, null, 'k');
		longopts[longopts_index++] = new LongOpt("realignmentThresholdLowerValue", LongOpt.REQUIRED_ARGUMENT, null, 'k');
		longopts[longopts_index++] = new LongOpt("realignmentThresholdItterationMethod", LongOpt.REQUIRED_ARGUMENT, null, 'k');
		//realignmentAlwaysTerminals
		longopts[longopts_index++] = new LongOpt("realignmentTerminals", LongOpt.REQUIRED_ARGUMENT, null, 'k');
		longopts[longopts_index++] = new LongOpt("terminals", LongOpt.REQUIRED_ARGUMENT, null, 'k');
		longopts[longopts_index++] = new LongOpt("realignmentItterations", LongOpt.REQUIRED_ARGUMENT, null, 'k');


		longopts[longopts_index++] = new LongOpt("useTCS", LongOpt.REQUIRED_ARGUMENT, null, 'v');
		longopts[longopts_index++] = new LongOpt("TcoffeePath", LongOpt.REQUIRED_ARGUMENT, null, 'v');
		longopts[longopts_index++] = new LongOpt("TCSTempFilePath", LongOpt.REQUIRED_ARGUMENT, null, 'v');
		
		Getopt g = new Getopt("opal", argv, "k:a:b:c:d:e:f:g:hi:j:l:m:n:o:pqr:s:t:uv:w:y:z:1", longopts);		
		 
		String optName;
		while ((c = g.getopt()) != -1) {
			//System.out.println((char)c + " " + longopts_index);
            arg = g.getOptarg();
            
			switch (c)  {
	        	case 1:
            		if (g.getLongind() == -1)  optName = "";	
            		else  optName = longopts[g.getLongind()].getName(); 
            		
		            if (optName.equals("convert")) {
		            	justDoConvert = true;
		            }
		            break;
			
				case 'a':
	        		fileA = arg.toString();

	        		if (g.getLongind() == -1)  optName = "";	
            		else  optName = longopts[g.getLongind()].getName(); 

		            if (optName.equals("init_alignment")) {
		            	AlignmentMaker.initAlignmentProvided = true;
		            }
	        		break;	
				case 'b':
            		fileB = arg.toString();
            		break;	
            	case 'c' :
            		if(advising_configs != null){
            			LogWriter.stdErrLogln("If the configuration list file is specified, then you cannot also specify a single parameter.");
        				throw new GenericOpalException("If the configuration list file is specified, then you cannot also specify a single parameter.");
            		}
            		costName = arg.toString();
            		break;
            	case 'd':
            		if (g.getLongind() == -1)  optName = "";	
            		else  optName = longopts[g.getLongind()].getName(); 
            		
		            if (optName.equals("dna")) {
            			costName = "DNA"; 
		            } else if (optName.equals("dnaAG") ) {
//		            	costName = "DNA";
		            	CostMatrix.dnaSubAG = Integer.parseInt(arg.toString());
		            } else if (optName.equals("dnaCT") || optName.equals("dnaCU")) {
//		            	costName = "DNA";
		            	CostMatrix.dnaSubCT = Integer.parseInt(arg.toString());
		            } else { //distance_type
	            		if ( arg.toString().toLowerCase().equals("normcost"))  {
	    	            	Tree.distanceType = Tree.DistanceType.normcost;
	    	            	if (Tree.iterations==-1) Tree.iterations = 1;
	            		} else if ( arg.toString().toLowerCase().equals("pctid")) {
	    	            	Tree.distanceType = Tree.DistanceType.pctid;
	    	            	if (Tree.iterations==-1) Tree.iterations = 1;
	            		} else if ( arg.toString().toLowerCase().equals("kmer")) {
	    	            	Tree.distanceType = Tree.DistanceType.kmer;
	    	            	if (Tree.iterations==-1) Tree.iterations = 1;
	            		} else if ( arg.toString().toLowerCase().equals("kmer_normcost")) {
	    	            	Tree.distanceType = Tree.DistanceType.kmer_normcost;
	            			if (Tree.iterations==-1) Tree.iterations = 2;
	    	            } else {
	    	            	LogWriter.stdErrLogln("unrecognized distance type. Try 'normcost' or 'pctid'");
	    	            	System.exit(1);
	    	            }	
	            		if ( arg.toString().toLowerCase().equals("kmer") ||  
	            				arg.toString().toLowerCase().equals("kmer_normcost"))  ; 
	            			//conf.sc.fillCompressedAlph();
            		}
		            break;
            	case 'e' :
            		if(advising_configs != null){
            			LogWriter.stdErrLogln("If the configuration list file is specified, then you cannot also specify a single parameter.");
        				throw new GenericOpalException("If the configuration list file is specified, then you cannot also specify a single parameter.");
            		}
            		gammaTerm = Integer.parseInt(arg.toString());
            		break;            
            	case 'f' :
            		if(advising_configs != null){
            			LogWriter.stdErrLogln("If the configuration list file is specified, then you cannot also specify a single parameter.");
        				throw new GenericOpalException("If the configuration list file is specified, then you cannot also specify a single parameter.");
            		}
            		lambdaTerm = Integer.parseInt(arg.toString());
            		break;            
            	case 'h':
//		            optName = longopts[g.getLongind()].getName();
            		OpalLogWriter.printUsage();
            		System.exit(1);
  		          	break;

            	case 'g' :
            		if(advising_configs != null){
            			LogWriter.stdErrLogln("If the configuration list file is specified, then you cannot also specify a single parameter.");
        				throw new GenericOpalException("If the configuration list file is specified, then you cannot also specify a single parameter.");
            		}
            		gamma = Integer.parseInt(arg.toString());
            		break;            
            	case 'l' :
            		if(advising_configs != null){
            			LogWriter.stdErrLogln("If the configuration list file is specified, then you cannot also specify a single parameter.");
        				throw new GenericOpalException("If the configuration list file is specified, then you cannot also specify a single parameter.");
            		}
            		lambda = Integer.parseInt(arg.toString());
            		break;            
            	case 'i' :
            		//tree_iterations
            		Tree.iterations = Integer.parseInt(arg.toString());
            		break;         
            	case 'm':
            		if (g.getLongind() == -1)  optName = "";	
            		else  optName = longopts[g.getLongind()].getName(); 
            		
            		if (optName.equals("mixed_alignment_cutoff")) {
		            	Aligner.mixedAlignmentCutoff = Integer.parseInt(arg.toString());
		            } else { //		if (optName.equals("align_method")) {
	            		if ( arg.toString().toLowerCase().startsWith("exact"))
	            			Aligner.alignmentMethod = AlignmentType.exact;
	    	            else if ( arg.toString().toLowerCase().startsWith("profile"))
	    	            	Aligner.alignmentMethod = AlignmentType.profile;
	    	            else if ( arg.toString().toLowerCase().startsWith("mixed"))
	    	            	Aligner.alignmentMethod = AlignmentType.mixed;
	    	            else {
	    	            	LogWriter.stdErrLogln("unrecognized alignment method. Try 'exact', 'profile', or 'mixed'");
	    	            	System.exit(1);
	    	            }	
		            } 
    	            break;
            	case 'o':  
            		if (g.getLongind() == -1)  optName = "";	
            		else  optName = longopts[g.getLongind()].getName(); 
            		
            		if (optName.equals("out_best")) {
            			bestOutputFile = arg.toString();
            		} else if (optName.equals("out_best_include_pre")) {
            			bestOutputFileIncludePreRealignment = arg.toString();
            		} else if (optName.equals("out_config")) {
            			configOutputFile = arg.toString();
            		} else if (optName.equals("out_feature")) {
            			featureOutputFile = arg.toString();
            		} else if (optName.equals("out_prerealignment")) {
            			preRealignmentOutputFile = arg.toString();
            		} else if (optName.equals("out_prerealignment_best")) {
            			bestPreRealignmentOutputFile = arg.toString();
            		} else if (optName.equals("out_prerealignment_best_realignment")) {
            			bestPreRealignmentsRealignmentOutputFile = arg.toString();
            		}  else if (optName.equals("out_prerealignment_best_realignment_include_pre")) {
            			bestPreRealignmentsRealignmentOutputFileIncludePreRealignment = arg.toString();
        			} else if (optName.equals("out")) {
            			//if(configOutputFile == null) configOutputFile = arg.toString();
            			if(bestOutputFile == null) bestOutputFile = arg.toString();
		            	/*try {
		            		AlignmentMaker.setOutputName(arg.toString());
		            	} catch (Exception E) {
		            		LogWriter.stdErrLogln("Error creating file for output: " + arg.toString());
	                		System.exit(1);            		
		            	}*/
            		} else if (optName.equals("tree_order")) {
    		            AlignmentMaker.outputOrderMethod = OutputOrderType.tree;
		            } else if (optName.equals("input_order")) {
		            	AlignmentMaker.outputOrderMethod = OutputOrderType.input;
		            } else { //  -o, or --out_format
	            		if ( arg.toString().toLowerCase().startsWith("clustal"))
	            			AlignmentWriter.outFormat = OutputType.clustal;
	    	            else if ( arg.toString().toLowerCase().startsWith("fasta"))
	    	            	AlignmentWriter.outFormat = OutputType.fasta;
	    	            else {
	    	            	LogWriter.stdErrLogln("unrecognized output format. Try 'clustalw' or 'fasta'");
	    	            	System.exit(1);
	    	            }	
		            }
    	            break;
            	case 'p':
            		if (g.getLongind() == -1)  optName = "";	
            		else  optName = longopts[g.getLongind()].getName(); 
            		
		            if (optName.equals("protein")) {
            			costName = "VTML200"; 
		            } else if (optName.equals("just_pair_dists")) {
		            	Tree.justPWDists = true;
		            	if (Tree.iterations==-1) Tree.iterations = 1;
		            } else if (optName.equals("polish_align_method")) {
		            	if (arg.toString().equals("profile")) {
		            		Polisher.polishAligmentMethod = AlignmentType.profile;
		            	} else if (arg.toString().equals("exact")) {
		            		Polisher.polishAligmentMethod = AlignmentType.exact;
		            	} else if (arg.toString().equals("mixed")) {
		            		Polisher.polishAligmentMethod = AlignmentType.mixed;
		            	} else {
		            		LogWriter.stdErrLogln("unknown polish aligner");
		            		System.exit(1);		            		
		            	}
		            } else if (optName.equals("just_polish")) {
		            	Polisher.justPolish = true;
		            } else { //(optName.equals("polish"))  or just "p"
		            	if (arg.toString().equals("exhaust_twocut")) {
		            		Polisher.polishMethod = Polisher.PolishType.exhaust_twocut;
		            	} else if (arg.toString().equals("exhaust_threecut")) {
		            		Polisher.polishMethod = Polisher.PolishType.exhaust_threecut;
		            	} else if (arg.toString().equals("random_tree_twocut")) {
			            	Polisher.polishMethod = Polisher.PolishType.random_tree_twocut;
		            	} else if (arg.toString().equals("random_twocut")) {
		            		Polisher.polishMethod = Polisher.PolishType.random_twocut;
		            	} else if (arg.toString().equals("random_threecut")) { 
		            		Polisher.polishMethod = Polisher.PolishType.random_threecut;
		            	} else if (arg.toString().equals("none")) { 
		            		Polisher.polishMethod = null;
		            	} else {
		            		LogWriter.stdErrLogln("Unknown polish method. Try exhaust_twocut, random_tree_twocut, random_twocut, random_threecut");
		            		System.exit(1);
		            	}
		            	
		            }
    	            break;
            	case 'q':
            		if (g.getLongind() == -1)  optName = "";	
            		else  optName = longopts[g.getLongind()].getName(); 
            		
            		if (optName.equals("silent")) {
		            	verbosity = -1;
		            }else if (optName.equals("quiet")) {
		            	verbosity = 0;
		            } else if (optName.equals("normal")) {
		            	verbosity = 1;
		            } else if (optName.equals("verbose")) {
		            	verbosity = 2;
		            } else if (optName.equals("noisy")) {
		            	verbosity = 3;
		            }
		            break;
            	case 'r' :
            		if (g.getLongind() == -1)  optName = "";	
            		else  optName = longopts[g.getLongind()].getName(); 
            		
		            if (optName.equals("polish_reps_exhaustive")) {
	            		Polisher.polishIterations_exact = Integer.parseInt(arg.toString());
		            } else if (optName.equals("pess_do_reverse")) {
	            		doReverse = true;
		            } else {// reps
		            	Polisher.polishIterations = Integer.parseInt(arg.toString());
		            }
		            if ( ! optName.equals("pess_do_reverse")) {
		            	if (Polisher.polishMethod==null) Polisher.polishMethod = Polisher.PolishType.random_tree_twocut;
		            }
            		break;            
            	case 's' :
            		if (g.getLongind() == -1)  optName = "";	
            		else  optName = longopts[g.getLongind()].getName(); 
            		
		            if (optName.equals("seed")) {
	            		Polisher.setRandomSeed( Long.parseLong(arg.toString()));
		            } else if (optName.equals("show_cost")) {
	            		AlignmentMaker.showCost = true;
		            } else if (optName.endsWith("subopt")) {
		            	PairSuboptimalityMatrices.setDelta( (Integer.valueOf(arg.toString())).intValue());
		            	if (optName.equals("just_subopt")) {
		            		justDoSubOpt = true;
		            	}
		            } else if (optName.startsWith("structure_file")) {
		            	
		            	/*if(advising_configs == null){
		            		advising_configs = new Configuration[1];
		            		advising_configs[0] = new Configuration();
		            	}
		            	if (!advising_configs[0].useStructure) // if model hasn't already been set
		            		StructureAlignment.setParams(StructureAlignment.ParamModel.G8,advising_configs[0]);
		            		
                  advising_configs[0].useStructure = true;
                  */
		            	
		              if (structmodel == null) // if model hasn't already been set
		            		structmodel = StructureAlignment.ParamModel.G8;
		            	

	            		Aligner.linearCutoff = Integer.MAX_VALUE;
		            	if (optName.equals("structure_file")) {
		            		if(structFileA!=null){
			            		LogWriter.stdErrLogln("Only 'facet_structure' or 'structure_file' can be specified. Alignment structure is used for Facet.");
		        				throw new GenericOpalException("Only 'facet_structure' or 'structure_file' can be specified. Alignment structure is used for Facet");
			            	}
		            		
		            		structFileA = arg.toString();
		            	} else if (optName.equals("structure_file2")) {
			            	structFileB = arg.toString();
		            	}
		            } else if (optName.equals("structure_model")) {
		            	/*if(advising_configs == null){
		            		advising_configs = new Configuration[1];
		            		advising_configs[0] = new Configuration();
		            	}
		            	
		            	advising_configs[0].useStructure = true;*/
		            	String s = arg.toString();
		            	boolean isSet = false;
		            	for (StructureAlignment.ParamModel type : StructureAlignment.ParamModel.values()) {
		            		if (s.equalsIgnoreCase(type.toString())) {
		            			structmodel = type;
		            			//StructureAlignment.setParams(type, advising_configs[0]);
		            			isSet = true;
		            			break;
		            		}
		            	}
		            	if (!isSet) {
		            		LogWriter.stdErrLogln("Unknown structure model");
		            		System.exit(1);
		            	}
		            } else if (optName.equals("facet_structure")){
		            	if(structFileA!=null){
		            		LogWriter.stdErrLogln("Only 'facet_structure' or 'structure_file' can be specified. Alignment structure is used for Facet.");
	        				throw new GenericOpalException("Only 'facet_structure' or 'structure_file' can be specified. Alignment structure is used for Facet");
		            	}
		            	structFileA = arg.toString();
		            	//configs[0].useStructure = false;
		            }
            		break;            
            	case 't' :
            		if (g.getLongind() == -1)  optName = "";	
            		else  optName = longopts[g.getLongind()].getName(); 
            		
            		if ( optName.equals("treein")) { 
            			EnteredTree.setTreeFromFile(arg.toString());
            			Tree.treeType = Tree.TreeType.entered; 
            		} else if ( optName.equals("treeout")) {
            			AlignmentMaker_SingleSequences.setTreeOutFile(arg.toString());
            		} else if (optName.equals("just_tree")) {
		            	justTree = true;
		            	Tree.iterations = 1;
		            	AlignmentMaker_SingleSequences.setTreeOutFile(arg.toString());
		            }
            		break;            
            	case 'u':
            		if (g.getLongind() == -1)  optName = "";	
            		else  optName = longopts[g.getLongind()].getName(); 
            		
            		if ( optName.equals("usage")) { 
            			OpalLogWriter.printUsage();
            			System.exit(1);
            		} else {
            			toUpper = true;
            		}
            		break;
            	case 'y':
            		AlignmentMaker_SingleSequences.consistency = true;

            		if (g.getLongind() == -1)  optName = "";	
            		else  optName = longopts[g.getLongind()].getName(); 
            		
		            if (optName.equals("consistency_weight")) {
		            	PairwiseAlignmentsContainer.consistency_weight = Float.parseFloat(arg.toString());
		            } else if (optName.equals("consistency_other_seqs_weight")) {
			            PairwiseAlignmentsContainer.consistency_other_seqs_weight = Float.parseFloat(arg.toString());
		            } else if (optName.equals("consistency_neighbors")) {
		            	PairwiseAlignmentsContainer.neighborCount = Integer.parseInt(arg.toString());
		            } else if (optName.equals("consistency_maxsubtree")) {
		            	PairwiseAlignmentsContainer.maxSubtreeSize = Integer.parseInt(arg.toString());
		            } else if (optName.equals("consistency_badscore_mult")) {
		            	PairwiseAlignmentsContainer.badScoreMult = Double.parseDouble(arg.toString());
		            } else if (optName.equals("consistency_use_neighbor_weights")) {
		            	PairwiseAlignmentsContainer.useWeights = true;
		            } else if (optName.equals("consistency_use_avg")) {
		            	PairwiseAlignmentsContainer.useMax = false;
		            } else if (optName.equals("consistency_flatten_abc_subopt")) {
		            	PairwiseAlignmentsContainer.flattenABC = Float.parseFloat(arg.toString());
		            } else if (optName.equals("consistency_neighbor_dist_thresh")) {
		            	PairwiseAlignmentsContainer.neighborDistThreshold = Float.parseFloat(arg.toString());
		            } else if (optName.equals("consistency_blend_type")) {
		            	if (arg.toString().equals("simple")) {
		            		PairwiseAlignmentsContainer.blendMethod = PairwiseAlignmentsContainer.BlendType.simple; 
		            	} else if (arg.toString().equals("asym_oneparam")) {
		            		PairwiseAlignmentsContainer.blendMethod = PairwiseAlignmentsContainer.BlendType.asym_oneparam;
		            	} else if (arg.toString().equals("asym_twoparam")) {
		            		PairwiseAlignmentsContainer.blendMethod = PairwiseAlignmentsContainer.BlendType.asym_twoparam;
		            	} else if (arg.toString().equals("symmetric")) {
		            		PairwiseAlignmentsContainer.blendMethod = PairwiseAlignmentsContainer.BlendType.symmetric;
		            	}
		            } else if (optName.equals("consistency_align_method")) {
		            	if ( arg.toString().toLowerCase().startsWith("exact"))
	            			ConsistencyAligner.alignmentMethod = AlignmentType.exact;
	    	            else if ( arg.toString().toLowerCase().startsWith("profile"))
	    	            	ConsistencyAligner.alignmentMethod = AlignmentType.profile;
	    	            else {
	    	            	LogWriter.stdErrLogln("unknown consistnecy aligner method");
	    	            	System.exit(1);
	    	            }
		            }	            
		            
		            break;                        	
            	case 'w' :
            		if (g.getLongind() == -1)  optName = "";	
            		else  optName = longopts[g.getLongind()].getName(); 
            		
		            if (optName.equals("output_width")) {
		            	AlignmentMaker.outputWidth = Integer.parseInt(arg.toString());
		            }
		            break;            
            	case 'z' :
            		Aligner.linearCutoff = Integer.parseInt(arg.toString());
            		break; 
            		
            	case 'k':	
            		if (g.getLongind() == -1)  optName = "";	
            		else  optName = longopts[g.getLongind()].getName(); 
            		// --realignmentMinimumSizeTypeValue --realignmentMinimumSizeValue 10 --realignmentWindowSizeValue 3 --realignmentThresholdTypeAverage --realignmentThresholdValue 0.5
            		if (optName.equals("realignmentMinimumSizeType")){
            			if (arg.toString().equals("multiplier")){ temp_window_minimum_type = Configuration.WINDOW_SIZE_MINIMUM.WINDOW_MULTIPLIER; }
            			if (arg.toString().equals("none")){ temp_window_minimum_type = Configuration.WINDOW_SIZE_MINIMUM.NONE; }
						if (arg.toString().equals("value")){ temp_window_minimum_type = Configuration.WINDOW_SIZE_MINIMUM.VALUE; }
            		}
					if (optName.equals("realignmentMinimumSizeValue")){ temp_realign_minimum_value = Float.parseFloat(arg.toString()); }
					
					
					if (optName.equals("realignmentWindowSizeType")){ 
						if (arg.toString().equals("value")){ temp_window_type = Configuration.WINDOW_SIZE.VALUE; }
						if (arg.toString().equals("percentage")){ temp_window_type = Configuration.WINDOW_SIZE.PERCENTAGE; }
					}
					if (optName.equals("realignmentWindowSizeValue")){ temp_realign_window_size = Float.parseFloat(arg.toString()); }
					if (optName.equals("realignmentMaximumWindowSizeValue")){ temp_maximum_realign_window_size = Float.parseFloat(arg.toString()); }
					if (optName.equals("realignmentMinimumWindowSizeValue")){ temp_minimum_realign_window_size = Float.parseFloat(arg.toString()); }
					
					if (optName.equals("realignmentThresholdType")){
						if (arg.toString().equals("value")){ temp_threshold_type = Configuration.THRESHOLD_TYPE.VALUE; }
						if (arg.toString().equals("average")){ temp_threshold_type = Configuration.THRESHOLD_TYPE.AVERAGE_WINDOW; }
            			if (arg.toString().equals("whole")){temp_threshold_type = Configuration.THRESHOLD_TYPE.WHOLE_ALIGNMENT;}
            			if (arg.toString().equals("two_value")){ temp_threshold_type = Configuration.THRESHOLD_TYPE.TWO_VALUE; }
            			if (arg.toString().equals("two_average")){ temp_threshold_type = Configuration.THRESHOLD_TYPE.TWO_AVERAGE; }
            			if (arg.toString().equals("two_sd")){ temp_threshold_type = Configuration.THRESHOLD_TYPE.TWO_SD; }
            			if (arg.toString().equals("two_whole")){temp_threshold_type = Configuration.THRESHOLD_TYPE.TWO_WHOLE;}
            			if (arg.toString().equals("two_percentage")){temp_threshold_type = Configuration.THRESHOLD_TYPE.TWO_PERCENTAGE;}
					}
            		if (optName.equals("realignmentThresholdValue")){ temp_realign_threshold = Float.parseFloat(arg.toString());}
            		if (optName.equals("realignmentThresholdLowerValue")){ temp_realign_threshold_lower = Float.parseFloat(arg.toString());}

            		if (optName.equals("realignmentWindowWeightDecay")){ temp_realignmentWindowWeightDecay = Float.parseFloat(arg.toString());}
            		
            		if (optName.equals("realignmentTerminals")){
            			if (arg.toString().equals("always")){ temp_realignment_terminals = Configuration.REALIGNMENT_TERMINALS.ALWAYS; /*LogWriter.stdErrLogln("always and positional terminals are dissabled due to an error"); System.exit(1);*/ }
            			if (arg.toString().equals("never")){ temp_realignment_terminals = Configuration.REALIGNMENT_TERMINALS.NEVER; }
            			if (arg.toString().equals("positional")){ temp_realignment_terminals = Configuration.REALIGNMENT_TERMINALS.POSITIONAL; /*LogWriter.stdErrLogln("always and positional terminals are dissabled due to an error"); System.exit(1);*/ }
            		}
            		if (optName.equals("realignmentItterations")){
            			temp_realignment_ittertations = Integer.parseInt(arg.toString());
            		}
            		//realignmentThresholdItterationMethod keepFirstItteration
            		if (optName.equals("realignmentThresholdItterationMethod")){
						if (arg.toString().equals("keepFirstItteration")){ temp_realignment_save_threshold = true; }
						if (arg.toString().equals("newThresholdEachItteration")){ temp_realignment_save_threshold = false; }
					}
            		if (optName.equals("terminals")){
						if (arg.toString().equals("all")){ temp_useLeftTerminal = true; temp_useRightTerminal = true; }
						if (arg.toString().equals("noLeft") || arg.toString().equals("Right")){ temp_useLeftTerminal = false; temp_useRightTerminal = true; }
						if (arg.toString().equals("noRight") || arg.toString().equals("Left")){ temp_useLeftTerminal = true; temp_useRightTerminal = false; }
						if (arg.toString().equals("none") || arg.toString().equals("noTerm")){ temp_useLeftTerminal = false; temp_useRightTerminal = false; }
					}
            		break; 
            	case 'v':	

            		if (g.getLongind() == -1)  optName = "";	
            		else  optName = longopts[g.getLongind()].getName(); 
            		
            		//longopts[longopts_index++] = new LongOpt("useTCS", LongOpt.REQUIRED_ARGUMENT, null, 'v');
            		//longopts[longopts_index++] = new LongOpt("TcoffeePath", LongOpt.REQUIRED_ARGUMENT, null, 'v');
            		//longopts[longopts_index++] = new LongOpt("TCSTempFilePath", LongOpt.REQUIRED_ARGUMENT, null, 'v');
            		if (optName.equals("useTCS")) {
		            	if ( arg.toString().toLowerCase().startsWith("both")){
		            		temp_useTCSforRealignment = true;
		            		temp_useTCSforAdvising = true;
		            	}
		            	else if ( arg.toString().toLowerCase().startsWith("realignment"))
		            		temp_useTCSforRealignment = true;
		            	else if ( arg.toString().toLowerCase().startsWith("advisisng"))
		            		temp_useTCSforAdvising = true;
		            }else if (optName.equals("TcoffeePath")) {
		            	temp_tcoffeeDirectory = arg.toString();
		            }else if (optName.equals("TCSTempFilePath")) {
		            	temp_temporaryFileDirectory = arg.toString();
		            }
            		
            		break;
            	case 'j':
            		if (g.getLongind() == -1)  optName = "";	
            		else  optName = longopts[g.getLongind()].getName(); 
            		
        			if(optName.equals("configuration_file") || optName.equals("advising_configuration_file") || optName.equals("realignment_configuration_file")){
        				if(optName.equals("configuration_file") || optName.equals("advising_configuration_file")){
		            		if(gamma != -1 || gammaTerm != -1 || lambda != -1 || lambdaTerm != -1 || !costName.equals("")){
		            			LogWriter.stdErrLogln("If the advisisng configuration list file is specified, then you cannot also specify a single parameter.");
		        				throw new GenericOpalException("If the advisisng configuration list file is specified, then you cannot also specify a single parameter.");
		            		}
		            		
		            		InputStream is = null;
		        			
		        			try {
		        				is = new FileInputStream(arg.toString());
		        			} catch (FileNotFoundException e) {
		        				LogWriter.stdErrLogln("The file '" + arg.toString() + "' cannot be found.  Qutting");
		        				throw new GenericOpalException(e.getMessage());
		        			}
		        			
		        			byte b[] = null;
		        			try {
		        				int x = is.available();
		        				b = new byte[x];
		        				is.read(b);
		        				is.close();
		        			} catch ( IOException e) {
		        				LogWriter.stdErrLogln("Error reading file '" + arg.toString() + "'");
		        				throw new GenericOpalException(e.getMessage());
		        				
		        			}
		        			String content = new String(b);
		        			//System.err.println("File:\n" + content);
		        			
		        			String[] cstrings = content.split("\n");
		        		
		        			advising_configs = new Configuration[cstrings.length];
		        			for(int i=0;i<cstrings.length;i++){
		        				//System.err.println("Configuration: " + cstrings[i]);
		        				advising_configs[i] = new Configuration(cstrings[i]);
		        			}
		        			
        				}
        				if(optName.equals("configuration_file") || optName.equals("realignment_configuration_file")){
		            		
		            		InputStream is = null;
		        			
		        			try {
		        				is = new FileInputStream(arg.toString());
		        			} catch (FileNotFoundException e) {
		        				LogWriter.stdErrLogln("The file '" + arg.toString() + "' cannot be found.  Qutting");
		        				throw new GenericOpalException(e.getMessage());
		        			}
		        			
		        			byte b[] = null;
		        			try {
		        				int x = is.available();
		        				b = new byte[x];
		        				is.read(b);
		        				is.close();
		        			} catch ( IOException e) {
		        				LogWriter.stdErrLogln("Error reading file '" + arg.toString() + "'");
		        				throw new GenericOpalException(e.getMessage());
		        			}
		        			String content = new String(b);
		        			//System.err.println("File:\n" + content);
		        			
		        			String[] cstrings = content.split("\n");
		        		
		        			realignment_configs = new Configuration[cstrings.length];
		        			for(int i=0;i<cstrings.length;i++){
		        				//System.err.println("Configuration: " + cstrings[i]);
		        				realignment_configs[i] = new Configuration(cstrings[i]);
		        			}
		        		
        				}
        			}else if(optName.equals("repeat_configurations")){
        				repeat_config = Integer.parseInt(arg.toString());
        			}
        			else if(optName.equals("max_threads")){
        				max_threads = Integer.parseInt(arg.toString());
        			}
        			else if(optName.equals("use_legacy_facet")){
        				useLegacyFacetFunction = true;
        			}
        			else if(optName.equals("use_updated_facet")){
        				useLegacyFacetFunction = false;
        			}
        			
        			
            		break;
            	default:
            		LogWriter.stdErrLogln("unrecognized option in command line");
            		System.exit(1);
            		break;
		          
			}

		}
		if (null == fileA) {
			if (argv.length>0  && g.getOptind()<argv.length && argv[g.getOptind()] != null) {
				fileA = argv[g.getOptind()].toString();
			}
		}


	}



	public String getFileA() {
		return fileA;
	}

	public String getFileB() {
		return fileB;
	}

	public String getStructFileA() {
		return structFileA;
	}

	public String getStructFileB() {
		return structFileB;
	}

	/*public String getCostName() {
		return costName;
	}
	
	public int getGamma() {
		return gamma;
	}

	public int getLambda() {
		return lambda;
	}

	public int getGammaTerm() {
		return gammaTerm;
	}

	public int getLambdaTerm() {
		return lambdaTerm;
	}*/
	
	/*Configuration.THRESHOLD_TYPE temp_threshold_type = null;
	Configuration.WINDOW_SIZE_MINIMUM temp_window_minimum_type = null;
	Configuration.WINDOW_SIZE temp_window_type = null;
	float temp_realign_threshold = -1;
	float temp_realign_window_size = -1;
	float temp_realign_minimum_value = -1;*/
	
	public Configuration[] getAdvisingConfigs(){
		if(advising_configs != null){
			for(int i=0;i<advising_configs.length;i++){
				if(temp_threshold_type != null) advising_configs[i].realignment_threshold_type = temp_threshold_type;
				if(temp_realign_threshold != -1) advising_configs[i].realignment_threshold_value = temp_realign_threshold;
				if(temp_realign_threshold_lower != (float)-10000.0) advising_configs[i].realignment_threshold_value_lower = temp_realign_threshold_lower;
				if(temp_window_minimum_type != null) advising_configs[i].realignment_minimum_type = temp_window_minimum_type;
				if(temp_realign_minimum_value != -1) advising_configs[i].realignment_minimum_window_value = temp_realign_minimum_value;
				if(temp_window_type != null) advising_configs[i].realignment_window_type = temp_window_type;
				if(temp_realign_window_size != -1) advising_configs[i].realignment_window_value = temp_realign_window_size;
				if(temp_realignment_terminals != null) advising_configs[i].realignment_use_terminals = temp_realignment_terminals;
				if(temp_realignmentWindowWeightDecay != -1) advising_configs[i].realignmentWindowWeightDecay = temp_realignmentWindowWeightDecay;
				advising_configs[i].useLegacyFacetFunction = useLegacyFacetFunction;
				advising_configs[i].doReverse = doReverse;
				advising_configs[i].useTCSforAdvising = temp_useTCSforAdvising;
				advising_configs[i].useTCSforRealignment = temp_useTCSforRealignment;
				if(temp_realignment_ittertations != -1) advising_configs[i].realignment_itterations = temp_realignment_ittertations;
				if(temp_tcoffeeDirectory != null) advising_configs[i].tcoffeeDirectory = temp_tcoffeeDirectory;
				if(temp_temporaryFileDirectory != null) advising_configs[i].temporaryFileDirectory = temp_temporaryFileDirectory;
				advising_configs[i].realignment_save_threshold = temp_realignment_save_threshold;
				if(temp_maximum_realign_window_size != -1) advising_configs[i].maximum_realignment_window_value = temp_maximum_realign_window_size;
				if(temp_minimum_realign_window_size != -1) advising_configs[i].minimum_realignment_window_value = temp_minimum_realign_window_size;
				advising_configs[i].useLeftTerminal = temp_useLeftTerminal;
				advising_configs[i].useRightTerminal = temp_useRightTerminal;
			}
			if(repeat_config==1){
				return advising_configs;
			}else{
				Configuration temp[] = new Configuration[advising_configs.length * repeat_config];
				for(int i=0;i<advising_configs.length;i++){
					for(int j=0;j<repeat_config;j++){
						temp[repeat_config*i + j] = new Configuration(advising_configs[i]);
						temp[repeat_config*i + j].repetition = j+1;
					}
				}
				return temp;
			}
		}
		
		// advising configs is null, so use the input information
		advising_configs = new Configuration[repeat_config];
		for(int j=0;j<repeat_config;j++){
		  
		  if(structmodel != null) {
		    //System.out.println("Structure Model Not Null");
        if(gamma == -1) gamma = StructureAlignment.getDefaultGamma(structmodel);
        if(gammaTerm == -1) gammaTerm = StructureAlignment.getDefaultGammaTerm(structmodel);
        if(lambda == -1) lambda = StructureAlignment.getDefaultLambda(structmodel);
        if(lambdaTerm == -1) lambdaTerm = StructureAlignment.getDefaultLambdaTerm(structmodel);
		  }else {

        //System.out.println("Structure Model IS Null");
		  }
			advising_configs[j] = new Configuration(costName, gamma, gammaTerm, lambda, lambdaTerm, fileA);
			if(structmodel != null) {
			  advising_configs[j].useStructure = true;
			  StructureAlignment.setParams(structmodel, advising_configs[j]);
			}
			
			if(temp_threshold_type != null) advising_configs[j].realignment_threshold_type = temp_threshold_type;
			if(temp_realign_threshold != -1) advising_configs[j].realignment_threshold_value = temp_realign_threshold;
			if(temp_realign_threshold_lower != (float)-10000.0) advising_configs[j].realignment_threshold_value_lower = temp_realign_threshold_lower;
			if(temp_window_minimum_type != null) advising_configs[j].realignment_minimum_type = temp_window_minimum_type;
			if(temp_realign_minimum_value != -1) advising_configs[j].realignment_minimum_window_value = temp_realign_minimum_value;
			if(temp_window_type != null) advising_configs[j].realignment_window_type = temp_window_type;
			if(temp_realign_window_size != -1) advising_configs[j].realignment_window_value = temp_realign_window_size;
			if(temp_realignment_terminals != null) advising_configs[j].realignment_use_terminals = temp_realignment_terminals;
			if(temp_realignmentWindowWeightDecay != -1) advising_configs[j].realignmentWindowWeightDecay = temp_realignmentWindowWeightDecay;
			advising_configs[j].useLegacyFacetFunction = useLegacyFacetFunction;
			advising_configs[j].doReverse = doReverse;
			advising_configs[j].useTCSforAdvising = temp_useTCSforAdvising;
			advising_configs[j].useTCSforRealignment = temp_useTCSforRealignment;
			if(temp_realignment_ittertations != -1) advising_configs[j].realignment_itterations = temp_realignment_ittertations;
			if(temp_tcoffeeDirectory != null) advising_configs[j].tcoffeeDirectory = temp_tcoffeeDirectory;
			if(temp_temporaryFileDirectory != null) advising_configs[j].temporaryFileDirectory = temp_temporaryFileDirectory;
			advising_configs[j].realignment_save_threshold = temp_realignment_save_threshold;
			advising_configs[j].useLeftTerminal = temp_useLeftTerminal;
			advising_configs[j].useRightTerminal = temp_useRightTerminal;

			if(temp_maximum_realign_window_size != -1) advising_configs[j].maximum_realignment_window_value = temp_maximum_realign_window_size;
			if(temp_minimum_realign_window_size != -1) advising_configs[j].minimum_realignment_window_value = temp_minimum_realign_window_size;
			
			if(repeat_config>1){
				advising_configs[j].repetition = j+1;
			}
		}

		return advising_configs;
	}
	
	public Configuration[] getRealignmentConfigs(){
		if(realignment_configs != null){
			for(int i=0;i<realignment_configs.length;i++){
				realignment_configs[i].useLegacyFacetFunction = useLegacyFacetFunction;
			}
			
			if(repeat_config==1){
				return realignment_configs;
			}else{
				Configuration temp[] = new Configuration[realignment_configs.length * repeat_config];
				for(int i=0;i<realignment_configs.length;i++){
					for(int j=0;j<repeat_config;j++){
						temp[repeat_config*i + j] = new Configuration(realignment_configs[i]);
						temp[repeat_config*i + j].repetition = j+1;
					}
				}
				return temp;
			}
		}
		return realignment_configs;
	}
	

	public Inputs getInputs(){
		Inputs in = new Inputs();
		in.verbosity = getVerbosity();
		in.toUpper = isToUpper(); 
		in.justDoConvert = isJustDoConvert();
		in.justDoSubOpt= isJustDoSubOpt();
		in.justTree= isJustTree();
		in.fileA = getFileA();
		in.fileB = getFileB();
		in.structFileA = getStructFileA();
		in.structFileB = getStructFileB();
		in.configOutputFile = configOutputFile;
		in.bestOutputFile = bestOutputFile;
		in.bestOutputFileIncludePreRealignment = bestOutputFileIncludePreRealignment;
		in.featureOutputFile = featureOutputFile;
		in.preRealignmentOutputFile = preRealignmentOutputFile;
		in.bestPreRealignmentOutputFile = bestPreRealignmentOutputFile;
		in.bestPreRealignmentsRealignmentOutputFile = bestPreRealignmentsRealignmentOutputFile;
		in.bestPreRealignmentsRealignmentOutputFileIncludePreRealignment = bestPreRealignmentsRealignmentOutputFileIncludePreRealignment;
        if (in.structFileA != null)
			in.structure = new StructureFileReader(in.structFileA, in.structFileB);
		return in;
	}
	
	public int getMaxThreads(){
		return max_threads;
	}
	
	public int getVerbosity() {
		return verbosity;
	}

	public boolean isToUpper() {
		return toUpper;
	}

	public boolean isJustDoConvert() {
		return justDoConvert;
	}

	public boolean isJustDoSubOpt() {
		return justDoSubOpt;
	}

	public boolean isJustTree() {
		return justTree;
	}

	
	
}

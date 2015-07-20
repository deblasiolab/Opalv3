package opal.IO;

import java.io.FileOutputStream;
import java.io.PrintStream;

import opal.IO.AlignmentWriter.OutputType;
import opal.align.Aligner;
import opal.align.ConsistencyAligner;
import opal.align.PairSuboptimalityMatrices;
import opal.align.PairwiseAlignmentsContainer;
import opal.align.ProfileAligner;
import opal.align.StructureAlignment;
import opal.align.Aligner.AlignmentType;
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
		
	int verbosity = 1;
	boolean toUpper = false; 
	int c;
	String arg;
	boolean justDoConvert= false;
	boolean justDoSubOpt= false;
	boolean justTree= false;

	
	
	public ArgumentHandler (String argString) {
		this(argString.split("\\s"));
	}
	
	public ArgumentHandler (String[] argv) {
		
		LogWriter.stdErrLogln("Got these arguments: ");
		for (String s : argv) {
			LogWriter.stdErrLogln(s);
		}
		
		
		
		LongOpt[] longopts = new LongOpt[58];
		longopts[0] = new LongOpt("help", LongOpt.NO_ARGUMENT, null, 'h');
		longopts[1] = new LongOpt("cost", LongOpt.REQUIRED_ARGUMENT, null, 'c');
		longopts[2] = new LongOpt("in", LongOpt.REQUIRED_ARGUMENT, null, 'a');
		longopts[3] = new LongOpt("in2", LongOpt.REQUIRED_ARGUMENT, null, 'b'); // if this is included, it's a request to align two alignments
		longopts[4] = new LongOpt("gamma", LongOpt.REQUIRED_ARGUMENT, null, 'g');
		longopts[5] = new LongOpt("lambda", LongOpt.REQUIRED_ARGUMENT, null, 'l');
		longopts[6] = new LongOpt("gamma_term", LongOpt.REQUIRED_ARGUMENT, null, 'e');
		longopts[7] = new LongOpt("lambda_term", LongOpt.REQUIRED_ARGUMENT, null, 'f');
		longopts[8] = new LongOpt("quiet", LongOpt.NO_ARGUMENT, null, 'q');
		longopts[9] = new LongOpt("align_method", LongOpt.REQUIRED_ARGUMENT, null, 'm'); // exact, profile
		longopts[10] = new LongOpt("out", LongOpt.REQUIRED_ARGUMENT, null, 'o'); // fasta, clustalw
		longopts[11] = new LongOpt("out_format", LongOpt.REQUIRED_ARGUMENT, null, 'o'); // fasta, clustalw 
		longopts[12] = new LongOpt("linear_cutoff", LongOpt.REQUIRED_ARGUMENT, null, 'z'); // default = 20
		longopts[13] = new LongOpt("output_width", LongOpt.REQUIRED_ARGUMENT, null, 'w'); // default depends on output type: 55 for clustal, 80 for fasta
		longopts[14] = new LongOpt("upper_case", LongOpt.NO_ARGUMENT, null, 'u'); // default retain input case 
		longopts[15] = new LongOpt("seed", LongOpt.REQUIRED_ARGUMENT, null, 's'); // random seed for polishing  
		longopts[16] = new LongOpt("polish_reps", LongOpt.REQUIRED_ARGUMENT, null, 'r'); // # repetitions of polishing
		longopts[17] = new LongOpt("polish", LongOpt.REQUIRED_ARGUMENT, null, 'p');
		longopts[18] = new LongOpt("convert", LongOpt.NO_ARGUMENT, null, 1);
		longopts[19] = new LongOpt("subopt", LongOpt.REQUIRED_ARGUMENT, null, 's');
		longopts[20] = new LongOpt("just_subopt", LongOpt.REQUIRED_ARGUMENT, null, 's');
		longopts[21] = new LongOpt("treein", LongOpt.REQUIRED_ARGUMENT, null, 't');
		longopts[22] = new LongOpt("treeout", LongOpt.REQUIRED_ARGUMENT, null, 't');
		longopts[23] = new LongOpt("just_tree", LongOpt.REQUIRED_ARGUMENT, null, 't');
		longopts[24] = new LongOpt("distance_type", LongOpt.REQUIRED_ARGUMENT, null, 'd'); // normcost, pctid
		longopts[25] = new LongOpt("usage", LongOpt.NO_ARGUMENT, null, 'u');
		longopts[26] = new LongOpt("consistency_badscore_mult", LongOpt.REQUIRED_ARGUMENT, null, 'y');
		longopts[27] = new LongOpt("consistency_use_avg", LongOpt.NO_ARGUMENT, null, 'y');
		longopts[28] = new LongOpt("consistency_use_neighbor_weights", LongOpt.NO_ARGUMENT, null, 'y');
		longopts[29] = new LongOpt("consistency_neighbors", LongOpt.REQUIRED_ARGUMENT, null, 'y');
		longopts[30] = new LongOpt("consistency_maxsubtree", LongOpt.REQUIRED_ARGUMENT, null, 'y');
		longopts[31] = new LongOpt("consistency_flatten_abc_subopt", LongOpt.REQUIRED_ARGUMENT, null, 'y');
		longopts[32] = new LongOpt("consistency_weight", LongOpt.REQUIRED_ARGUMENT, null, 'y');
		longopts[33] = new LongOpt("consistency_other_seqs_weight", LongOpt.REQUIRED_ARGUMENT, null, 'y');
		longopts[34] = new LongOpt("consistency_blend_type", LongOpt.REQUIRED_ARGUMENT, null, 'y');
		longopts[35] = new LongOpt("use_consistency", LongOpt.NO_ARGUMENT, null, 'y');
		longopts[36] = new LongOpt("consistency_neighbor_dist_thresh", LongOpt.REQUIRED_ARGUMENT, null, 'y');
		longopts[37] = new LongOpt("just_pair_dists", LongOpt.NO_ARGUMENT, null, 'p');
		longopts[38] = new LongOpt("tree_order", LongOpt.NO_ARGUMENT, null, 'o');
		longopts[39] = new LongOpt("input_order", LongOpt.NO_ARGUMENT, null, 'o');
		longopts[40] = new LongOpt("polish_align_method", LongOpt.REQUIRED_ARGUMENT, null, 'p');
		longopts[41] = new LongOpt("polish_reps_exhaustive", LongOpt.REQUIRED_ARGUMENT, null, 'r'); // # repetitions of exhaustive polishing (only worth doing if first polish is with heuristic) 
		longopts[42] = new LongOpt("consistency_align_method", LongOpt.REQUIRED_ARGUMENT, null, 'y');
		longopts[43] = new LongOpt("mixed_alignment_cutoff", LongOpt.REQUIRED_ARGUMENT, null, 'm');
		longopts[44] = new LongOpt("pess_do_reverse", LongOpt.NO_ARGUMENT, null, 'r');
		longopts[45] = new LongOpt("show_cost", LongOpt.NO_ARGUMENT, null, 's');
		longopts[46] = new LongOpt("tree_iterations", LongOpt.REQUIRED_ARGUMENT, null, 'i');
		longopts[47] = new LongOpt("structure_file", LongOpt.REQUIRED_ARGUMENT, null, 's');
		longopts[48] = new LongOpt("structure_model", LongOpt.REQUIRED_ARGUMENT, null, 's');
		longopts[49] = new LongOpt("dna", LongOpt.NO_ARGUMENT, null, 'd');
		longopts[50] = new LongOpt("protein", LongOpt.NO_ARGUMENT, null, 'p');
		longopts[51] = new LongOpt("dnaAG", LongOpt.REQUIRED_ARGUMENT, null, 'd');
		longopts[52] = new LongOpt("dnaCT", LongOpt.REQUIRED_ARGUMENT, null, 'd');
		longopts[53] = new LongOpt("dnaCU", LongOpt.REQUIRED_ARGUMENT, null, 'd');
		longopts[54] = new LongOpt("verbose", LongOpt.NO_ARGUMENT, null, 'q');
		longopts[55] = new LongOpt("noisy", LongOpt.NO_ARGUMENT, null, 'q');
		longopts[56] = new LongOpt("init_alignment", LongOpt.REQUIRED_ARGUMENT, null, 'a');
		longopts[57] = new LongOpt("just_polish", LongOpt.NO_ARGUMENT, null, 'p');
		
		Getopt g = new Getopt("opal", argv, "a:b:c:d:e:f:g:hi:jl:m:o:p:r:s:t:w:y:z:huq", longopts);		
		 
		String optName;
		while ((c = g.getopt()) != -1) {
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
	            				arg.toString().toLowerCase().equals("kmer_normcost")) 
	            			SequenceConverter.fillCompressedAlph();
            		}
		            break;
            	case 'e' :
            		gammaTerm = Integer.parseInt(arg.toString());
            		break;            
            	case 'f' :
            		lambdaTerm = Integer.parseInt(arg.toString());
            		break;            
            	case 'h':
//		            optName = longopts[g.getLongind()].getName();
            		OpalLogWriter.printUsage();
            		System.exit(1);
  		          	break;

            	case 'g' :
            		gamma = Integer.parseInt(arg.toString());
            		break;            
            	case 'l' :
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
            		
            		if (optName.equals("out")) {
		            	try {
		            		AlignmentMaker.setOutputName(arg.toString());
		            		PrintStream out = new PrintStream(new FileOutputStream(arg.toString()));
			                System.setOut(out);
		            	} catch (Exception E) {
		            		LogWriter.stdErrLogln("Error creating file for output: " + arg.toString());
	                		System.exit(1);            		
		            	}
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
            			costName = "BLOSUM62"; 
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
            		
            		if (optName.equals("quiet")) {
		            	verbosity = 0;
		            } else if (optName.equals("verbose")) {
		            	verbosity = 1;
		            } else if (optName.equals("noisy")) {
		            	verbosity = 2;
		            }
		            break;
            	case 'r' :
            		if (g.getLongind() == -1)  optName = "";	
            		else  optName = longopts[g.getLongind()].getName(); 
            		
		            if (optName.equals("polish_reps_exhaustive")) {
	            		Polisher.polishIterations_exact = Integer.parseInt(arg.toString());
		            } else if (optName.equals("pess_do_reverse")) {
	            		ProfileAligner.doReverse = true;
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
	            		Polisher.setRandomSeed( (new Long(arg.toString())).longValue());
		            } else if (optName.equals("show_cost")) {
	            		AlignmentMaker.showCost = true;
		            } else if (optName.endsWith("subopt")) {
		            	PairSuboptimalityMatrices.setDelta( (new Integer(arg.toString())).intValue());
		            	if (optName.equals("just_subopt")) {
		            		justDoSubOpt = true;
		            	}
		            } else if (optName.startsWith("structure_file")) {
		            	if (!Aligner.useStructure) // if model hasn't already been set
		            		StructureAlignment.setParams(StructureAlignment.ParamModel.UN4);

	            		Aligner.linearCutoff = Integer.MAX_VALUE;
		            	
		            	Aligner.useStructure = true;
		            	if (optName.equals("structure_file")) {
		            		structFileA = arg.toString();
		            	} else if (optName.equals("structure_file2")) {
			            	structFileB = arg.toString();
		            	}
		            } else if (optName.equals("structure_model")) {
		            	Aligner.useStructure = true;
		            	String s = arg.toString();
		            	boolean isSet = false;
		            	for (StructureAlignment.ParamModel type : StructureAlignment.ParamModel.values()) {
		            		if (s.equalsIgnoreCase(type.toString())) {
		            			StructureAlignment.setParams(type);
		            			isSet = true;
		            			break;
		            		}
		            	}
		            	if (!isSet) {
		            		LogWriter.stdErrLogln("Unknown structure model");
		            		System.exit(1);
		            	}
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
		            	PairwiseAlignmentsContainer.consistency_weight = (new Float(arg.toString())).floatValue();
		            } else if (optName.equals("consistency_other_seqs_weight")) {
			            PairwiseAlignmentsContainer.consistency_other_seqs_weight = (new Float(arg.toString())).floatValue();
		            } else if (optName.equals("consistency_neighbors")) {
		            	PairwiseAlignmentsContainer.neighborCount = (new Integer(arg.toString())).intValue();
		            } else if (optName.equals("consistency_maxsubtree")) {
		            	PairwiseAlignmentsContainer.maxSubtreeSize = (new Integer(arg.toString())).intValue();
		            } else if (optName.equals("consistency_badscore_mult")) {
		            	PairwiseAlignmentsContainer.badScoreMult = (new Double(arg.toString())).doubleValue();
		            } else if (optName.equals("consistency_use_neighbor_weights")) {
		            	PairwiseAlignmentsContainer.useWeights = true;
		            } else if (optName.equals("consistency_use_avg")) {
		            	PairwiseAlignmentsContainer.useMax = false;
		            } else if (optName.equals("consistency_flatten_abc_subopt")) {
		            	PairwiseAlignmentsContainer.flattenABC = (new Float(arg.toString())).floatValue();
		            } else if (optName.equals("consistency_neighbor_dist_thresh")) {
		            	PairwiseAlignmentsContainer.neighborDistThreshold = (new Float(arg.toString())).floatValue();            
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
		            	AlignmentMaker.outputWidth = (new Integer(arg.toString())).intValue();
		            }
		            break;            
            	case 'z' :
            		Aligner.linearCutoff = (new Integer(arg.toString())).intValue();
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

	public String getCostName() {
		return costName;
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

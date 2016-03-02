package opal.realignment;

import java.util.Arrays;

import facet.Facet;
import facet.FacetAlignment;
import opal.IO.Configuration;
import opal.IO.TCS;
import opal.makers.AlignmentMaker;
import opal.makers.AlignmentMaker_SingleSequences;
import opal.IO.Inputs;

public class realignmentDriver {
	char[][] sequence;
	float[][][] structure_prob;
	Configuration[] configList;
	Configuration globalConfiguration;
	float wholeAlignmentScore;
	AlignmentMaker globalAligmentMaker;
	float[][][] originalStrucutre;
	
	float fixed_good_threshold = -1;
	float fixed_bad_threshold = -1;
	float fixed_threshold = -1;
	
	int[][] alignmentInstance;
	boolean changed_one_regions;
	
	public static float mean(float[] nums) {
		float sum = 0;
		 
		for (int i = 0; i < nums.length; i++) {
			sum += nums[i];
		}
	
		return sum / nums.length;
	} 
	private float standardDeviation(float[] numbers){
		float mean = mean(numbers);
		float squareSum = 0;

		for (int i = 0; i < numbers.length; i++) {
			squareSum += Math.pow(numbers[i] - mean, 2);
		}

		return (float)Math.sqrt((squareSum) / (numbers.length - 1));
	}
	
	private float sum(float[] nums){
		float sum = 0;
		 
		for (int i = 0; i < nums.length; i++) {
			sum += nums[i];
		}
	
		return sum;
	}
	
	private void derive_structure_array_from_sequence(){
		structure_prob = new float[sequence.length][sequence[0].length][3];
		for(int i=0;i<sequence.length;i++){
            int k=0;
            for(int j=0;j<sequence[0].length;j++){
                    if(sequence[i][j] == '-'){
                            structure_prob[i][j][0]=-1;
                            structure_prob[i][j][1]=-1;
                            structure_prob[i][j][2]=-1;
                    }
                    else{
                            structure_prob[i][j][0] = originalStrucutre[i][k][0];
                            structure_prob[i][j][1] = originalStrucutre[i][k][1];
                            structure_prob[i][j][2] = originalStrucutre[i][k][2];
                            k++;
                    }
            }
		}
	}
	
	public realignmentDriver(char[][] se, float[][][] sp, Configuration[] cList, Configuration gConfig, float score, AlignmentMaker gAM) {
		sequence = se.clone();
		wholeAlignmentScore = score;
		globalConfiguration = gConfig;
		globalAligmentMaker = gAM;
		originalStrucutre = sp;
		
		derive_structure_array_from_sequence();
		
		configList = cList.clone();
		for(int i=0;i<configList.length;i++){
			configList[i].useLeftTerminal = false;
			configList[i].useRightTerminal = false;
		}
	}
	
	private float evaluateWindow(int startIndex, int endIndex){
		char[][] windowSequences = new char[sequence.length][endIndex-startIndex+1];
		float[][][] windowStructure = new float[sequence.length][][];
		for(int i=0; i<sequence.length; i++){
			int numberOfNonGaps = 0;
			for(int j=startIndex; j<=endIndex; j++){
				if(sequence[i][j] != '-') numberOfNonGaps++;
				windowSequences[i][j-startIndex] = sequence[i][j];
			}
			windowStructure[i] = new float[numberOfNonGaps][3];
			int k=0;
			for(int j=startIndex; j<=endIndex; j++){
				if(sequence[i][j] != '-'){
					windowStructure[i][k][0] = structure_prob[i][j][0];
					windowStructure[i][k][1] = structure_prob[i][j][1];
					windowStructure[i][k][2] = structure_prob[i][j][2];
					k++;
				}
			}
		}
		float value = Facet.defaultValue(new FacetAlignment(windowSequences, windowStructure), globalConfiguration.useLegacyFacetFunction);
		if(globalConfiguration.useTCSforRealignment)
			value = (float)TCS.TCSValue(globalConfiguration.sc.convertSeqsToInts(windowSequences), globalAligmentMaker, globalConfiguration);
		return (Float.isNaN(value))?0:value;
	}

	private char[][] realignWindow(int startIndex, int endIndex){
		
		// because the region may be larger than the windows we tested earlier
		float originalScore = evaluateWindow(startIndex, endIndex);
		
		
		char[][] originalWindowSequences = new char[sequence.length][endIndex-startIndex+1];		
		int[] numberOfNonGaps = new int[sequence.length];
		int numberOfNonBlankSequences = 0;
		String realignmentIndexes = "";
		
		String[] names = new String[sequence.length];
		for(int i=0; i<sequence.length; i++){
			names[i] = Integer.toString(i);
			int stringStartChar = 0;
			int stringEndChar = 0;
			for(int j=0; j<startIndex; j++){
				if(sequence[i][j] != '-'){
					stringStartChar++;
					stringEndChar++;
				}
			}
			for(int j=startIndex; j<=endIndex; j++){
				if(sequence[i][j] != '-'){
					numberOfNonGaps[i]++;
					stringEndChar++;
				}
				originalWindowSequences[i][j-startIndex] = Character.toUpperCase(sequence[i][j]);;
			}
			if(numberOfNonGaps[i] > 0) numberOfNonBlankSequences++;
			realignmentIndexes += "("+stringStartChar+","+stringEndChar+");";
		}
		
		char[][] windowSequences = new char[numberOfNonBlankSequences][];
		float[][][] windowStructure = new float[numberOfNonBlankSequences][][];
		
		
		int putInIndex = 0;
		for(int i=0; i<sequence.length; i++){
			int location = 0;
			if(numberOfNonGaps[i] > 0){
				windowSequences[putInIndex] = new char[numberOfNonGaps[i]];
				for(int j=startIndex; j<=endIndex; j++){
					if(sequence[i][j] != '-') windowSequences[putInIndex][location++] = Character.toUpperCase(sequence[i][j]);
				}
				windowStructure[putInIndex] = new float[numberOfNonGaps[i]][3];
				int k=0;
				for(int j=startIndex; j<=endIndex; j++){
					if(sequence[i][j] != '-'){
						windowStructure[putInIndex][k][0] = structure_prob[i][j][0];
						windowStructure[putInIndex][k][1] = structure_prob[i][j][1];
						windowStructure[putInIndex][k][2] = structure_prob[i][j][2];
						k++;
					}
				}
				putInIndex++;
			}
		}
		
		Inputs in = new Inputs();
		char[][] bestSubAlignment = null;
		float bestScore = originalScore;
		String bestConfig = "";
		
		for(int i=0; numberOfNonBlankSequences > 1 && i<configList.length; i++){
			
			configList[i].useLeftTerminal = (startIndex <= 0);
			configList[i].useRightTerminal = (endIndex >= sequence[0].length-1);
			if(globalConfiguration.realignment_use_terminals == Configuration.REALIGNMENT_TERMINALS.NEVER) 
				configList[i].useRightTerminal = configList[i].useLeftTerminal = false;
			if(globalConfiguration.realignment_use_terminals == Configuration.REALIGNMENT_TERMINALS.ALWAYS) 
				configList[i].useRightTerminal = configList[i].useLeftTerminal = true;
			
			/** Workaround, because I reverted back to single terminal gap penalties **/
			//configList[i].setGammaTerm(configList[i].gamma);
			//configList[i].setLambdaTerm(configList[i].lambda);
			
			//System.err.println("Window: [" + startIndex + "," + endIndex + "] " + configList[i].useLeftTerminal + " " + configList[i].useRightTerminal);
			
			AlignmentMaker_SingleSequences am = new AlignmentMaker_SingleSequences();
			am.initialize(configList[i].sc.convertSeqsToInts(windowSequences), names, configList[i], in);
			int[][] subAlignmentInstance = am.buildAlignment();
			//am.printOutput(subAlignmentInstance, null);
			float score = Facet.defaultValue(new FacetAlignment(configList[i].sc.convertIntsToSeqs(subAlignmentInstance),windowStructure), globalConfiguration.useLegacyFacetFunction);
			if(globalConfiguration.useTCSforRealignment)
				score = (float)TCS.TCSValue(subAlignmentInstance, am, globalConfiguration);
			
			//globalConfiguration.realignmentLog += "Best Score:" + bestScore + "\tcurrent score:" + score + "\n";
			if(score>bestScore){
				bestScore = score;
				char[][] tempBestSubAlignment = configList[i].sc.convertIntsToSeqs(subAlignmentInstance);
				bestSubAlignment = new char[sequence.length][tempBestSubAlignment[0].length];
				putInIndex = 0;
				for(int j=0; j<sequence.length; j++){
					for(int k=0;k<tempBestSubAlignment[0].length;k++){
						if(numberOfNonGaps[j] > 0){
							bestSubAlignment[j][k] = tempBestSubAlignment[putInIndex][k];
						}else{
							bestSubAlignment[j][k] = '-';
						}
					}
					if(numberOfNonGaps[j] > 0) putInIndex++;
				}
				bestConfig = "used " + configList[i].toString() + " (" + numberOfNonBlankSequences + " of " + sequence.length + " sequences)";
				changed_one_regions = true;
			}
		}
		if(bestSubAlignment == null){
			bestSubAlignment = originalWindowSequences;
			bestConfig = "kept the same";
		}
		if(numberOfNonBlankSequences <= 1) bestConfig = "didn't try";
		globalConfiguration.realignmentLog += "Realigned window [" + startIndex + "," + endIndex + "] " + bestConfig + "\n"+realignmentIndexes+"\n";
		
		
		return bestSubAlignment;
	}
	
	public void simpleRealignment(int windowSize){
		simpleRealignment(windowSize, 0);
	}
	public void simpleRealignment(int windowSize, int itteration){
		if(globalConfiguration.realignment_itterations <= itteration) return;
		
		float[] scores = new float[sequence[0].length-windowSize+1];
		float scoreTotal = 0;
		for(int i=0;i<=sequence[0].length-windowSize;i++){
			scores[i] = evaluateWindow(i,i+windowSize-1);
			scoreTotal += scores[i];
		}
		
		String[] newAlignment = new String[sequence.length];
		for(int i=0;i<sequence.length;i++){
			newAlignment[i] = "";
		}
		boolean[] includeInRealignment = new boolean[sequence[0].length];
		
		if(globalConfiguration.realignment_threshold_type == Configuration.THRESHOLD_TYPE.AVERAGE_WINDOW ||
				globalConfiguration.realignment_threshold_type == Configuration.THRESHOLD_TYPE.WHOLE_ALIGNMENT ||
						globalConfiguration.realignment_threshold_type == Configuration.THRESHOLD_TYPE.VALUE){
			// TODO threshold scores
			float threshold = 0;
			//(float)0.75 * wholeAlignmentScore;
			
			//threshold = (float)0.75 * (scoreTotal/scores.length);
			if(globalConfiguration.realignment_save_threshold && itteration>0){
				threshold = fixed_threshold;
			}else{
				if(globalConfiguration.realignment_threshold_type == Configuration.THRESHOLD_TYPE.AVERAGE_WINDOW) 
					threshold = (float)globalConfiguration.realignment_threshold_value * (scoreTotal/scores.length);
				if(globalConfiguration.realignment_threshold_type == Configuration.THRESHOLD_TYPE.WHOLE_ALIGNMENT) 
					threshold = (float)globalConfiguration.realignment_threshold_value * wholeAlignmentScore;
				if(globalConfiguration.realignment_threshold_type == Configuration.THRESHOLD_TYPE.VALUE)
					threshold = globalConfiguration.realignment_threshold_value;
				
				fixed_threshold = threshold;
			}
			System.err.println(scoreTotal + "/" + scores.length + " = " + (scoreTotal/scores.length) + "\t" + wholeAlignmentScore + "\t" 
					+ globalConfiguration.realignment_threshold_type + "," + globalConfiguration.realignment_threshold_value + "(" + threshold + ")\t" 
					+ globalConfiguration.realignment_window_type + "," + globalConfiguration.realignment_window_value + "\t" 
					+ globalConfiguration.realignment_minimum_type + "," + globalConfiguration.realignment_minimum_window_value
					);
			
			
			int[] countGood = new int[sequence[0].length];
			for(int i=0;i<sequence[0].length-windowSize;i++){
				if(scores[i]>threshold){
					for(int j=0;j<windowSize;j++){
						countGood[i+j]++;
					}
				}
			}
				
			
			
			for(int i=0;i<sequence[0].length;i++){
				includeInRealignment[i] = (countGood[i]<windowSize);
			}
		}else if(globalConfiguration.realignment_threshold_type == Configuration.THRESHOLD_TYPE.TWO_AVERAGE ||
				globalConfiguration.realignment_threshold_type == Configuration.THRESHOLD_TYPE.TWO_VALUE ||
				globalConfiguration.realignment_threshold_type == Configuration.THRESHOLD_TYPE.TWO_WHOLE ||
				globalConfiguration.realignment_threshold_type == Configuration.THRESHOLD_TYPE.TWO_SD ||
				globalConfiguration.realignment_threshold_type == Configuration.THRESHOLD_TYPE.TWO_PERCENTAGE){
			
			float window_column_value[] = new float[windowSize];
			for(int j=0;j<=windowSize/2;j++){
				window_column_value[j] = (float)Math.pow(globalConfiguration.realignmentWindowWeightDecay, (windowSize/2)-j);
				window_column_value[windowSize-j-1] = (float)Math.pow(globalConfiguration.realignmentWindowWeightDecay, (windowSize/2)-j);
			}
			
			float[] column_scores = new float[sequence[0].length];
			float[] column_totals = new float[sequence[0].length];
			for(int i=0;i<=sequence[0].length-windowSize;i++){
				for(int j=0;j<windowSize;j++){
					//System.err.println("Window " + i + "\tColumn " + (i+j));
					column_scores[i+j]+=scores[i]*window_column_value[j];
					column_totals[i+j]+=window_column_value[j];
				}
			
			}
			
			for(int i=0;i<sequence[0].length;i++){
				//System.err.print("Column " + i + ": " + column_scores[i] + "/" + column_totals[i]);
				column_scores[i] /= column_totals[i];
				//System.err.println("\t" + column_scores[i]);
			}
			
			float good_threshold = 0;
			float bad_threshold = 0;
			
			if(globalConfiguration.realignment_save_threshold && itteration>0){
				good_threshold = fixed_good_threshold;
				bad_threshold = fixed_bad_threshold;
			}else{
				if(globalConfiguration.realignment_threshold_type == Configuration.THRESHOLD_TYPE.TWO_AVERAGE){ 
					good_threshold = (float)globalConfiguration.realignment_threshold_value * (sum(column_scores)/column_scores.length);
					bad_threshold = (float)globalConfiguration.realignment_threshold_value_lower * (sum(column_scores)/column_scores.length);
				}
				if(globalConfiguration.realignment_threshold_type == Configuration.THRESHOLD_TYPE.TWO_WHOLE) {
					good_threshold = (float)globalConfiguration.realignment_threshold_value * wholeAlignmentScore;
					bad_threshold = (float)globalConfiguration.realignment_threshold_value_lower * wholeAlignmentScore;
				}
				if(globalConfiguration.realignment_threshold_type == Configuration.THRESHOLD_TYPE.TWO_SD) {
					float standard_deviation = standardDeviation(column_scores);
					good_threshold = (float)globalConfiguration.realignment_threshold_value * standard_deviation + (sum(column_scores)/column_scores.length);
					bad_threshold = (float)globalConfiguration.realignment_threshold_value_lower * standard_deviation + (sum(column_scores)/column_scores.length);
				}
				if(globalConfiguration.realignment_threshold_type == Configuration.THRESHOLD_TYPE.TWO_VALUE){
					good_threshold = globalConfiguration.realignment_threshold_value;
					bad_threshold = globalConfiguration.realignment_threshold_value_lower;
				}
				if(globalConfiguration.realignment_threshold_type == Configuration.THRESHOLD_TYPE.TWO_PERCENTAGE){
					float[] sort_scores = column_scores.clone();
					Arrays.sort(sort_scores);
					good_threshold = sort_scores[scores.length - (int)(scores.length * globalConfiguration.realignment_threshold_value) - 1];
					bad_threshold = sort_scores[(int)(scores.length * globalConfiguration.realignment_threshold_value_lower) -1];
				}
				
				fixed_good_threshold = good_threshold;
				fixed_bad_threshold = bad_threshold;
			}
			int ct_good = 0;
			int ct_bad = 0;
			
			for(int i=0;i<sequence[0].length;i++){
				if(column_scores[i] >= good_threshold){
					includeInRealignment[i] = false;
					System.err.print("G");
					ct_good++;
				}
				else if(column_scores[i] <= bad_threshold){
					System.err.print("B");
					ct_bad++;
					includeInRealignment[i] = true;
					for(int j=i-1;j>=0 && column_scores[j] < good_threshold;j--){
						includeInRealignment[j] = true;
					}
					for(int j=i+1;j<sequence[0].length && column_scores[j] < good_threshold;j++){
						includeInRealignment[j] = true;
					}
				}else{
					System.err.print("-");
				}
			}

			System.err.println();
			System.err.println("good: " + ct_good + "\tbad: " + ct_bad + "\ttotal: " + sequence[0].length + "\tgood_threshold: " + good_threshold + "\tbad_threshold: " + bad_threshold);
		}
		//System.err.println();
		int ct_save = 0;
		for(int i=0;i<sequence[0].length;i++){
			System.err.print(includeInRealignment[i]?"X":" ");
			ct_save += (includeInRealignment[i]?0:1);
		}
		System.err.println();
		System.err.println("saved: " + ct_save + "\trealigned: " + (sequence[0].length - ct_save) + "\ttotal: " + sequence[0].length);
		
		int[] blockStart = new int[sequence[0].length];
		int blockAssignI = 1;
		for(int i=1;i<sequence[0].length;i++){
			if(includeInRealignment[i]!=includeInRealignment[blockStart[blockAssignI-1]]){
				blockStart[blockAssignI] = i;
				blockAssignI++;
			}
		}
		
		
		if(globalConfiguration.realignment_minimum_type != Configuration.WINDOW_SIZE_MINIMUM.NONE){
			// TODO merge low scoring windows
			int minimumBlockSize = 0;
			if(globalConfiguration.realignment_minimum_type == Configuration.WINDOW_SIZE_MINIMUM.VALUE)
				minimumBlockSize = (int)globalConfiguration.realignment_minimum_window_value;
			if(globalConfiguration.realignment_minimum_type == Configuration.WINDOW_SIZE_MINIMUM.WINDOW_MULTIPLIER)
				minimumBlockSize = (int)globalConfiguration.realignment_minimum_window_value * windowSize;
			
			//2*windowSize;
			
			//remove good blocks that are too short
			for(int i=0;i<blockAssignI-1;i++){
				if(blockStart[i+1] - blockStart[i] < minimumBlockSize && !includeInRealignment[blockStart[i]]){
					for(int j=blockStart[i];j<blockStart[i+1]; j++){
						includeInRealignment[j] = true;
					}
				}
			}
			
			//redo block definitions
			blockStart = new int[sequence[0].length];
			blockAssignI = 1;
			for(int i=1;i<sequence[0].length;i++){
				if(includeInRealignment[i]!=includeInRealignment[blockStart[blockAssignI-1]]){
					blockStart[blockAssignI] = i;
					blockAssignI++;
				}
			}
			//remove back blocks that are too short
			for(int i=0;i<blockAssignI-1;i++){
				if(blockStart[i+1] - blockStart[i] < minimumBlockSize && includeInRealignment[blockStart[i]]){
					for(int j=blockStart[i];j<blockStart[i+1]; j++){
						includeInRealignment[j] = false;
					}
				}
			}
		}
		
		

		changed_one_regions = false;
		globalConfiguration.realignmentLog += "Realignment Itteration " + itteration + " \n-------------------\n";
		
		// TODO realign using each configuration
		int lastInAlignmentAlready = -1;
		for(int i=0;i<sequence[0].length;i++){
			if(!includeInRealignment[i]){
				if(lastInAlignmentAlready+1 < i){
					char[][] realignedRegion = realignWindow(lastInAlignmentAlready+1, i-1);
					for(int j=0; j<sequence.length; j++){
						for(int k=0; k<realignedRegion[j].length;k++)
						newAlignment[j] += realignedRegion[j][k];
					}
				}
				
				for(int j=0; j<sequence.length; j++){
					newAlignment[j] += sequence[j][i];
				}
				lastInAlignmentAlready = i;
			}
		}
		//if the alignment ends in a realignment region, process it;
		if(includeInRealignment[sequence[0].length-1]){
			char[][] realignedRegion = realignWindow(lastInAlignmentAlready+1, sequence[0].length-1);
			for(int j=0; j<sequence.length; j++){
				for(int k=0; k<realignedRegion[j].length;k++)
				newAlignment[j] += realignedRegion[j][k];
			}
		}
		
		// TODO merge new alignments
		char[][] newAlignmentChar = new char[newAlignment.length][];
		for(int j=0; j<sequence.length; j++){
			newAlignmentChar[j] = newAlignment[j].toCharArray();
		}
		sequence = newAlignmentChar;
		derive_structure_array_from_sequence();
		alignmentInstance = globalConfiguration.sc.convertSeqsToInts(newAlignmentChar);
		
		if(changed_one_regions){
			simpleRealignment(windowSize, itteration + 1);
		}
	}
	
	public int[][] newAlignment(){
		return alignmentInstance;
	}
}


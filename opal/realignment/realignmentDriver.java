package opal.realignment;

import facet.Facet;
import facet.FacetAlignment;
import opal.IO.Configuration;
import opal.makers.AlignmentMaker_SingleSequences;
import opal.makers.AlignmentMaker;
import opal.IO.Inputs;

public class realignmentDriver {
	char[][] sequence;
	float[][][] structure_prob;
	Configuration[] configList;
	Configuration globalConfiguration;
	float wholeAlignmentScore;
	
	int[][] alignmentInstance;
	
	public realignmentDriver(char[][] se, float[][][] sp, Configuration[] cList, Configuration gConfig, float score) {
		sequence = se.clone();
		structure_prob = new float[se.length][se[0].length][3];
		wholeAlignmentScore = score;
		globalConfiguration = gConfig;
		
		for(int i=0;i<se.length;i++){
                int k=0;
                for(int j=0;j<se[0].length;j++){
                        if(sequence[i][j] == '-'){
                                structure_prob[i][j][0]=-1;
                                structure_prob[i][j][1]=-1;
                                structure_prob[i][j][2]=-1;
                        }
                        else{
                                structure_prob[i][j][0] = sp[i][k][0];
                                structure_prob[i][j][1] = sp[i][k][1];
                                structure_prob[i][j][2] = sp[i][k][2];
                                k++;
                        }
                }
        }
		
		configList = cList.clone();
		for(int i=0;i<configList.length;i++){
			configList[i].gammaTerm = configList[i].gamma;
			configList[i].lambdaTerm = configList[i].lambda;
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
		return Facet.defaultValue(new FacetAlignment(windowSequences, windowStructure));
	}

	private char[][] realignWindow(int startIndex, int endIndex){
		
		char[][] windowSequences = new char[sequence.length][endIndex-startIndex+1];
		float[][][] windowStructure = new float[sequence.length][][];
		String[] names = new String[sequence.length];
		for(int i=0; i<sequence.length; i++){
			int numberOfNonGaps = 0;
			names[i] = Integer.toString(i);
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
		
		Inputs in = new Inputs();
		char[][] bestSubAlignment = null;
		float bestScore = -1;
		
		for(int i=0; i<configList.length; i++){
			AlignmentMaker_SingleSequences am = new AlignmentMaker_SingleSequences();
			am.initialize(configList[i].sc.convertSeqsToInts(windowSequences), names, configList[i], in);
			int[][] subAlignmentInstance = am.buildAlignment();
			float score = Facet.defaultValue(new FacetAlignment(configList[i].sc.convertIntsToSeqs(subAlignmentInstance),windowStructure));
			if(bestScore == -1 || score>bestScore){
				bestScore = score;
				bestSubAlignment = configList[i].sc.convertIntsToSeqs(subAlignmentInstance);
			}
		}
		
		return bestSubAlignment;
	}
	
	public void simpleRealignment(int windowSize){
		float[] scores = new float[sequence[0].length-windowSize];
		for(int i=0;i<sequence[0].length-windowSize;i++){
			scores[i] = evaluateWindow(i,i+windowSize);
		}
		
		String[] newAlignment = new String[sequence.length];
		
		// TODO threshold scores
		float threshold = wholeAlignmentScore;
		int[] countGood = new int[sequence[0].length];
		for(int i=0;i<sequence[0].length-windowSize;i++){
			if(scores[i]>threshold){
				for(int j=0;j<windowSize;j++){
					countGood[i+j]++;
				}
			}
		}
			
		// TODO merge low scoring windows
		int minimumBlockSize = 2*windowSize;
		boolean[] includeInRealignment = new boolean[sequence[0].length];
		for(int i=0;i<sequence[0].length;i++){
			includeInRealignment[i] = (countGood[i]>=windowSize);
		}
		int[] blockStart = new int[sequence[0].length];
		int blockAssignI = 1;
		for(int i=1;i<sequence[0].length;i++){
			if(includeInRealignment[i]!=includeInRealignment[blockStart[blockAssignI-1]]){
				blockStart[blockAssignI] = i;
				blockAssignI++;
			}
		}
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
		
		
		
		// TODO realign using each configuration
		int lastInAlignmentAlready = -1;
		for(int i=0;i<sequence[0].length;i++){
			if(!includeInRealignment[i]){
				char[][] realignedRegion = realignWindow(lastInAlignmentAlready+1, i-1);
				for(int j=0; j<sequence.length; j++){
					for(int k=0; k<realignedRegion[j].length;k++)
					newAlignment[j] += realignedRegion[j][k];
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
		alignmentInstance = globalConfiguration.sc.convertSeqsToInts(newAlignmentChar);
	}
	
	public int[][] newAlignment(){
		return alignmentInstance;
	}
}


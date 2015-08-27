package opal.IO;

import java.io.*;
import java.util.regex.Pattern;

import opal.align.StructureAlignment;
import opal.exceptions.GenericOpalException;
import com.traviswheeler.libs.LogWriter;

public class StructureFileReader {

	public static String[] names;
	public static float[][] helices;
	public static float[][] sheets;
	public static float[][] loops;
	public static float[][] structureLevels;
	public static float[][] structureNeighborLevels;
	public static float[][][] structure;

	public static int window = 0;
	
	private static String getContentFromFname(String filename) throws GenericOpalException{
		InputStream is = null;
		
		try {
			is = new FileInputStream(filename);
		} catch (FileNotFoundException e) {
			LogWriter.stdErrLogln("The file '" + filename + "' cannot be found.  Qutting");
			throw new GenericOpalException(e.getMessage());
		}
	
		byte b[] = null;
		try {
			int x = is.available();
			b = new byte[x];
			is.read(b);
			is.close();
		} catch ( IOException e) {
			LogWriter.stdErrLogln("Error reading file '" + filename + "'");
			throw new GenericOpalException(e.getMessage());
		}
		
		return new String(b);
		
	}
	
	public static void initialize (String filename){
		initialize(filename, null);
	}
	
	public static void initialize (String filename, String filename2) {
		
		String content = getContentFromFname(filename);
		if(filename2 != null) content += "\n" + getContentFromFname(filename2);
		
		String[] sequences = content.split(">");
		
		names = new String[sequences.length-1];
		helices = new float[sequences.length-1][];
		sheets = new float[sequences.length-1][];
		loops = new float[sequences.length-1][];;
		structure = new float[sequences.length-1][][];
		structureLevels = new float[sequences.length-1][];
		structureNeighborLevels  = new float[sequences.length-1][];
		
		float[][] helicesTmp = new float[sequences.length-1][];
		float[][] sheetsTmp = new float[sequences.length-1][];
		float[][] loopsTmp = new float[sequences.length-1][];

		int preWindow = (window-1)/2;
		int postWindow = window/2;
		
		for (int i=0; i<sequences.length-1; i++ ) {

			int eol =	sequences[i+1].indexOf("\n");
			names[i] = sequences[i+1].substring(0,eol).trim();

			String[] triples = sequences[i+1].substring(eol+1).split("\n");
			int len = triples.length;
		
			structure[i] = new float[len][3];
			structureLevels[i] = new float[len];
			structureNeighborLevels[i]  = new float[len];
			
			loops[i] = new float[len];
			helices[i] = new float[len];
			sheets[i] = new float[len]; 
			loopsTmp[i] = new float[len];
			helicesTmp[i] = new float[len];
			sheetsTmp[i] = new float[len]; 
			
			for (int j=0; j<len; j++) {
		        Pattern p = Pattern.compile("[\\s]+");
		        String[] vals = p.split(triples[j]);
		        
		        
		        loopsTmp[i][j] = Float.parseFloat(vals[0]);
		        helicesTmp[i][j] = Float.parseFloat(vals[1]);
		        sheetsTmp[i][j] = Float.parseFloat(vals[2]);
		        
		        float tot = loopsTmp[i][j] + helicesTmp[i][j] + sheetsTmp[i][j];
		     
		        loopsTmp[i][j] /=  tot;
		        helicesTmp[i][j] /= tot;
		        sheetsTmp[i][j] /= tot;
		        
			}

			for (int j=0; j<len; j++) {
				int cnt = 1 + (j - Math.max(0, j-preWindow)) + (Math.min(len-1, j+postWindow) - j ) ;
				float loopSum=loopsTmp[i][j];
				float helixSum=helicesTmp[i][j];
				float sheetSum=sheetsTmp[i][j];
				for (int k=Math.max(0, j-preWindow); k<j; k++) {
					loopSum += loopsTmp[i][k];
					helixSum += helicesTmp[i][k];
					sheetSum += sheetsTmp[i][k];					
				}
				for (int k=j+1; k<=Math.min(len-1, j+postWindow); k++) {
					loopSum += loopsTmp[i][k];
					helixSum += helicesTmp[i][k];
					sheetSum += sheetsTmp[i][k];					
				}
				
//		        loops[i][j] = loopSum/cnt; 
//		        helices[i][j] = helixSum/cnt; 
//		        sheets[i][j] = sheetSum/cnt; 

				//3 digitd of precision
		        loops[i][j] =  (float)Math.round(1000*loopSum/cnt) / 1000;
		        helices[i][j] = (float)Math.round(1000*helixSum/cnt) / 1000;
		        sheets[i][j] = (float)Math.round(1000*sheetSum/cnt) / 1000;


		        structure[i][j][0] = (float)Math.round(1000*loopSum/cnt) / 1000;;
		        structure[i][j][1] = (float)Math.round(1000*helixSum/cnt) / 1000;;
		        structure[i][j][2] = (float)Math.round(1000*sheetSum/cnt) / 1000;;

			}

			for (int j=0; j<len-1; j++) {
				float structElemProb = helices[i][j] + sheets[i][j];
				//int lvl = (int)Math.floor( structElemProb * conf.gapLevelCnt);
				//if (lvl == conf.gapLevelCnt) lvl--;
				structureLevels[i][j] = structElemProb;

				float structElemProbNeighbor = (structElemProb + helices[i][j] + sheets[i][j])/2;
				//lvl = (int)Math.floor( structElemProbNeighbor * conf.gapLevelCnt);
				//if (lvl == conf.gapLevelCnt) lvl--;
				structureNeighborLevels[i][j] = structElemProbNeighbor;
			}
		
		}
	}

}

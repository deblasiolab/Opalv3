/**
 * 
 */
package opal.IO;

import opal.makers.AlignmentMaker;

import java.io.IOException;
import java.util.*;
/**
 * @author deblasio
 *
 */
public class TCS {
	private static String execCmd(String cmd) throws java.io.IOException {
        java.util.Scanner s = new java.util.Scanner(Runtime.getRuntime().exec(cmd).getInputStream());
        
        String rtn = "";
        while(s.hasNext()){
        	String line = s.nextLine();
        	if(line.startsWith("SCORE=")){
        		rtn = line.substring(6);
        		break;
        	}
        }
        return rtn;
    }
	
	public static double TCSValue(int[][] alignmentInstance, AlignmentMaker am, Configuration config){
		String filename = config.temporaryFileDirectory + "tmp" + alignmentInstance.hashCode() + "_" + am.hashCode() + "_" + config.hashCode()  + "_" + am.in.hashCode();
		
		int oldVerbosity = am.in.verbosity;
		am.in.verbosity = -1;
		am.printOutput(alignmentInstance, filename);
		am.in.verbosity = oldVerbosity;
		
		String command = config.tcoffeeDirectory + "t_coffee -infile " + filename + " -evaluate -method proba_pair -output score_ascii -outfile=stdout";
		String result = "-1";
		try {
			result = execCmd(command);
			execCmd("rm " + filename);
		} catch (NumberFormatException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		//System.err.println("Command: \""+command+"\"\nPre Parse Line: \""+result+"\"");
		return Double.parseDouble(result)/100.0;
	}
}
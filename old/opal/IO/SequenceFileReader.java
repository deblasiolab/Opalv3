package opal.IO;

import java.io.*;

import opal.exceptions.GenericOpalException;
import com.traviswheeler.libs.LogWriter;

public class SequenceFileReader {
	char[][] seqs;	
	String[] names;

	public SequenceFileReader (String filename, boolean separateSequences) {
		
		InputStream is = null;
		
		try {
			if (null == filename) {
				
				is = System.in;
				int x = 0;
				try {
					x= is.available();					
				} catch ( IOException e) {
					LogWriter.stdErrLogln("Error reading from STDIN");
					throw new GenericOpalException(e.getMessage());			
				}

				if (x==0) {  //STILL??  what am I supposed to align?
					LogWriter.stdErrLogln("No sequences to align. Quitting");
					//logger.printUsage();
					throw new GenericOpalException("No sequences to align. Quitting");
				}

			} else {
				is = new FileInputStream(filename);
			}
		} catch (FileNotFoundException e) {
			LogWriter.stdErrLogln("The file '" + filename + "' cannot be found.  Qutting");
			throw new GenericOpalException("The file '" + filename + "' cannot be found.  Qutting");
		}
	
		byte b[] = null;
		try {
			int x= is.available();
			b= new byte[x];
			is.read(b);
		} catch ( IOException e) {
			LogWriter.stdErrLogln("Error reading file '" + filename + "'");
			throw new GenericOpalException(e.getMessage());			
		}
		String content = new String(b);
//		content = content.replaceAll("-","");	
		
		String[] tmp = content.split(">");
		char[][] tmpChars = new char[tmp.length-1][];
		names = new String[tmp.length-1];
		int length = 0;
		boolean cleanDashCols = true;
		for (int i=1; i<tmp.length; i++ ) {
			int eol =	tmp[i].indexOf("\n");
			names[i-1] = tmp[i].substring(0,eol).trim();

			
			if (separateSequences) {
				tmpChars[i-1] = tmp[i].substring(eol+1).replaceAll(" |\r|\n|>|-","").toCharArray();
			} else { 
				tmpChars[i-1] = tmp[i].substring(eol+1).replaceAll(" |\r|\n|>","").toCharArray();
				if (i==1)
					length = tmpChars[i-1].length;
				else if  (length != tmpChars[i-1].length ) {
					LogWriter.stdErrLogln( "Invalid input. Input sequences are not all of the same length. (" + filename + ")" );
					throw new GenericOpalException("Invalid input. Input sequences are not all of the same length. (" + filename + ")");
				}
			}			
		}
				
		
		
		//clean out columns with only dashes
		if (separateSequences) {
			seqs = tmpChars;
		} else {
			int k;
			if (cleanDashCols) {
				k=0;
				for (int i=0; i<tmpChars[0].length; i++){ //columns
					boolean good = false;
					for (int j=0; j<tmpChars.length; j++) {
						if (tmpChars[j][i] != '-') {
							good = true;
							break;
						}
					}
					if (good) { 
						for (int j=0; j<tmpChars.length; j++) {
							tmpChars[j][k] = tmpChars[j][i]; 
						}
						k++;
					}			
				}
			
				seqs = new char[tmp.length-1][k];
				for (int i=0; i<k; i++)  //columns
					for (int j=0; j<tmpChars.length; j++) 
						seqs[j][i] = tmpChars[j][i];
			}
		}
	}
	
	public char[][] getSeqs () {
		return seqs;
	}

	public String[] getNames () {
		return names;
	}
}

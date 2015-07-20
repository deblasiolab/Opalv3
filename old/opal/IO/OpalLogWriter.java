package opal.IO;

import com.traviswheeler.libs.LogWriter;
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;

import opal.Opal;
import opal.exceptions.GenericOpalException;


public class OpalLogWriter extends LogWriter {

	private static String version = "";
	
	
	public OpalLogWriter () {
		Package p = this.getClass().getPackage();	
		version = p.getImplementationVersion();
	}

	public static void printVersion ( ) {
		OpalLogWriter l = new OpalLogWriter(); // to get version info
		logger.logln( "\nOpal v" + version + " by Travis Wheeler and John Kececioglu ");
		
	}
	
	public static void printUsage () {
		printVersion();
		
		try {
			InputStream is = Opal.class.getResourceAsStream("IO/usage.txt");
			BufferedReader br = new BufferedReader(new InputStreamReader(is));
			String line;
		    while (null != (line = br.readLine())) {
		    	logger.logln(line);
		     }
		} catch (FileNotFoundException e) {		
			logger.logln("The usage file cannot be found.  Qutting");
			throw new GenericOpalException("The usage file cannot be found.  Qutting");
		} catch (IOException e) {		
			logger.logln("Problem reading usage file cannot be found.  Qutting");
			throw new GenericOpalException("Problem reading usage file cannot be found.  Qutting");
		}
		
		logger.logln("");
	}

	
}

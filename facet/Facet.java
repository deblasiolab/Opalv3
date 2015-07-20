package facet;

import java.io.FileNotFoundException;
import java.io.File;
import java.util.*;
import java.lang.StringBuffer;
import java.io.PrintWriter;
/**
 * 
 */

/**
 * @author deblasio
 *
 */
public class Facet {

	public static float defaultValue(FacetAlignment a){
		double total = 0.171730397;
		Configuration c = new Configuration();
		c.normalizeBigN = true;
		total += 0.105434529 * PercentIdentity.replacement_score(a, c);
		total += 0.172122922 * GapDensity.open(a, c);
		total += 0.174107269 * Blockiness.evaluate(a, c);
		total += 0.176015402 * PercentIdentity.structure(a, c);
		total += 0.200589481 * Support.probability(a, c);
		 return (float) total;
	}
	
	public static void outputDefaultFeatures(String fname, FacetAlignment a){
		try{
			File f = new File(fname);
			if(f.getParentFile()!=null) f.getParentFile().mkdirs();
			PrintWriter writer = new PrintWriter(f, "UTF-8");
			Configuration c = new Configuration();
			writer.println(PercentIdentity.replacement_score(c) + "\t" + PercentIdentity.replacement_score(a, c));
			writer.println(GapDensity.open(c) +"\t"+ GapDensity.open(a, c));
			writer.println(Blockiness.evaluate(c) +"\t"+ Blockiness.evaluate(a, c));
			writer.println(PercentIdentity.structure(c) +"\t"+ PercentIdentity.structure(a, c));
			writer.println(Support.probability(c) +"\t"+ Support.probability(a, c));
			writer.close();
		}catch(Exception e){
			System.err.println("Error while writing feature file " + fname);
			System.exit(11);
		}
	}

}

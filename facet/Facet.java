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
					
			Configuration c6 = new Configuration();
			Configuration c10 = new Configuration();
			Configuration c20 = new Configuration();
			c6.normalizeBigN = true;
			c10.normalizeBigN = true;
			c20.normalizeBigN = true;
			c6.equivelanceClassSize = 6;
			c10.equivelanceClassSize = 10;
			c20.equivelanceClassSize = 20;
			if(a.type == FacetAlignment.AlignmentType.DNA){
				c6.matrix = Configuration.ReplacementMatrix.DNA;
				c10.matrix = Configuration.ReplacementMatrix.DNA;
				c20.matrix = Configuration.ReplacementMatrix.DNA;
			}
			if(a.type == FacetAlignment.AlignmentType.RNA){
				c6.matrix = Configuration.ReplacementMatrix.RNA;
				c10.matrix = Configuration.ReplacementMatrix.RNA;
				c20.matrix = Configuration.ReplacementMatrix.RNA;
			}
			
			writer.println(PercentIdentity.replacement_score(c20) + "\t" + PercentIdentity.replacement_score(a, c20));
			writer.println(GapDensity.open(c6) + "\t" + GapDensity.open(a, c6));
			writer.println(GapDensity.extension(c6) + "\t" + GapDensity.extension(a, c6));
			writer.println(GapDensity.consensus(c6) + "\t" + GapDensity.consensus(a, c6));
			writer.println(PercentIdentity.sequence(c10) + "\t" + PercentIdentity.sequence(a, c10));
			writer.println(CoreColumn.percentage(c10) + "\t" + CoreColumn.percentage(a, c10));
			writer.println(CoreColumn.consensus(c20) + "\t" + CoreColumn.consensus(a, c20));
			writer.println(InformationContent.evaluate(c20) + "\t" + InformationContent.evaluate(a, c20));
			
			if(a.type == FacetAlignment.AlignmentType.Protein){
				writer.println(GapCoil.percentage(c6) + "\t" + GapCoil.percentage(a, c6));
				writer.println(Blockiness.evaluate(c6) + "\t" + Blockiness.evaluate(a, c6));
				writer.println(PercentIdentity.structure(c10) + "\t" + PercentIdentity.structure(a, c10));
				writer.println(Support.probability(c10) + "\t" + Support.probability(a, c10));
			}

			writer.close();

		}catch(Exception e){
			System.err.println("Error while writing feature file " + fname);
			System.err.println(e.toString());
			System.exit(11);
		}
	}

}

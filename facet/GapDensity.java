package facet;

public class GapDensity {
	private static int number_gaps;
	private static int number_opens;
	private static int total_width;
	
	public static String open(Configuration c){
		String rtn = "gapOpenDensity";
		if(c.normalizeBigN) rtn += "-bigN";
		return rtn;
	}
	public static float open(FacetAlignment a, Configuration c){
		run_loops(a,c);
		if(c.normalizeBigN){
			return (float) 1.0-(number_opens/a.bigN);
		}else{
			return (float) 1.0-(number_opens/number_gaps);
		}
	}

	public static String extension(Configuration c){
		String rtn = "gapExtensionDensity";
		if(c.normalizeBigN) rtn += "-bigN";
		return rtn;
	}
	public static float extension(FacetAlignment a, Configuration c){
		run_loops(a,c);
		if(c.normalizeBigN){
			return (float) 1.0-(number_gaps/a.bigN);
		}else{
			return (float) 1.0-(number_gaps/total_width);
		}
	}

	private static void run_loops(FacetAlignment a, Configuration c){
		number_gaps = 0 ; number_opens = 0 ; total_width = 0;
		for(int i=0;i<a.numberOfSequences;i++){
			boolean in_gap = false;
			for(int j=0;j<a.width;j++){
				total_width++;
				if(a.sequence[i].charAt(j) == '-'){
					if(!in_gap) number_opens++;
					number_gaps++;
					in_gap = true;
				}else{
					in_gap = false;
				}
			}
		}
	}
}

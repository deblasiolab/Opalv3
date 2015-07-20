package facet;

public class PercentIdentity {
	public static String sequence(Configuration c){
		String rtn = "percentIdentity";
		if(c.equivelanceClassSize!=20) rtn += "-E" + c.equivelanceClassSize;
		return rtn;
	}
	public static float sequence(FacetAlignment a, Configuration c){
		int pairs=0;
		int pairsMatch=0;
		for(int i1=0;i1<a.numberOfSequences;i1++){
			for(int i2=i1+1;i2<a.numberOfSequences;i2++){
				for(int j=0;j<a.width;j++){
					if(!(a.sequence[i1].charAt(j) == '-' && a.sequence[i2].charAt(j)=='-')){
						pairs++;
						if(c.basesEquivelant(a.sequence[i1].charAt(j), a.sequence[i2].charAt(j))){
							pairsMatch++;
						}
					}
				}
			}
		}
		return (float) pairsMatch/pairs;
	}
	public static String structure(Configuration c){
		return "structurePercentIdentity";
	}
	public static float structure(FacetAlignment a, Configuration c){
		int pairs=0;
		int pairsMatch=0;
		for(int i1=0;i1<a.numberOfSequences;i1++){
			for(int i2=i1+1;i2<a.numberOfSequences;i2++){
				for(int j=0;j<a.width;j++){
					if(!(a.structure[i1].charAt(j) == '-' && a.structure[i2].charAt(j)=='-')){
						pairs++;
						if(a.structure[i1].charAt(j)==a.structure[i2].charAt(j)){
							pairsMatch++;
						}
					}
				}
			}
		}
		return (float) pairsMatch/pairs;
	}
	
	public static String structure_prob(Configuration c){
		return "structurePercentIdentity-structureProbability";
	}
	public static float structure_prob(FacetAlignment a, Configuration c){
		int pairs=0;
		float pairsMatch=0;
		for(int i1=0;i1<a.numberOfSequences;i1++){
			for(int i2=i1+1;i2<a.numberOfSequences;i2++){
				for(int j=0;j<a.width;j++){
					if(!(a.structure_prob[i1][j][0] == -1 && a.structure_prob[i2][j][0]==-1)){
						pairs += a.structure_prob[i1][j][0] * a.structure_prob[i2][j][0];
						pairs += a.structure_prob[i1][j][1] * a.structure_prob[i2][j][1];
						pairs += a.structure_prob[i1][j][2] * a.structure_prob[i2][j][2];						
					}
				}
			}
		}
		return (float) pairsMatch/pairs;
	}
	
	public static String replacement_score(Configuration c){
		String rtn = "avgReplacementScore";
		if(c.normalizeBigN) rtn += "-bigN";
		return rtn;
	}
	
	public static float replacement_score(FacetAlignment a, Configuration c){
		int pairs=0;
		float pairsMatch=0;
		for(int i1=0;i1<a.numberOfSequences;i1++){
			for(int i2=i1+1;i2<a.numberOfSequences;i2++){
				for(int j=0;j<a.width;j++){
					if(!(a.sequence[i1].charAt(j) == '-' || a.sequence[i2].charAt(j)=='-')){
						pairs++;
						pairsMatch += c.replacementValue(a.sequence[i1].substring(j, j+1), a.sequence[i2].substring(j, j+1));
					}
				}
			}
		}
		if(c.normalizeBigN){
			return (float) pairsMatch/(a.bigN*100);
		}else{
			return (float) pairsMatch/(pairs*100);
		}
	}
}

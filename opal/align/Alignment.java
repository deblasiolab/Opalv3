package opal.align;

import java.util.Hashtable;

import com.traviswheeler.libs.LogWriter;
import opal.IO.SequenceConverter;
import opal.exceptions.GenericOpalException;
import opal.IO.Configuration;

public class Alignment {

	static int alphabetLength;

	//just a container for the sequences and all these profile frequency tables, 
	// so they can be passed around easily e.g. into Tree and Shape objects).
	// No getters/setters, since encapsulation isn't a concern but speed is.
	public int[] seqIds;
	public int[][] seqs;
	public int M, K;

	public int[][] chars;
	public int[][] freqs;
	public int[] f0, f1;
	public int[] firstLetterCount;
	public int[] lastLetterCount;
	public int[] gapsBeforeFirst;
	public int[] gapsAfterLast;
		
	public int[] firstLetterLoc;
	public int[] lastLetterLoc;

	public int[] positions;
	public int[][] posInUgappedString;
	
	public boolean isReverse;
	public int sequencesStartAtInStructure;

	public Hashtable<String,Integer> kmers = new Hashtable<String,Integer>();
	public int kmerK = 4;
	
	//extra profile arrays for profile alignment
	public int[] f00, f01, f10, f11;
			
	public Configuration conf;
		
	public Alignment (int[] A, int id, Configuration c){
		this(A,id,c,0);
	}
	
	public Alignment (int[] A, int id, Configuration c, int startIndex) { // single sequence		
		K = 1;
		seqs = new int[1][];
		seqs[0] = A;
		seqIds = new int[1];
		seqIds[0] = id + startIndex;
		M = A.length;
		conf = c;
		alphabetLength = conf.cost.getChars().length;
		sequencesStartAtInStructure = startIndex;
		
		
		if (M==0) {
			LogWriter.stdErrLogln("Opal can't align zero-length sequences ( seq #" + (id+1) + "). Quitting" );
			throw new GenericOpalException("Opal can't align zero-length sequences ( seq #" + (id+1) + "). Quitting" );
		}
			
		
		makeProfile();
	}

	public Alignment (int[][] A, int[] ids, boolean isReverse, Configuration c) {
		this(A,ids, isReverse, c, 0);
	}
	
	public Alignment (int[][] A, int[] ids, boolean isReverse, Configuration c, int startIndex) {
			
	
		this.isReverse = isReverse;
		conf = c;

		alphabetLength = conf.cost.getChars().length;
		
		K = A.length;
		M = A[0].length;
		seqs = A;
		seqIds = ids;
		for(int i=0; i<seqIds.length; i++){
			seqIds[i] += startIndex;
		}
		sequencesStartAtInStructure = startIndex;
		makeProfile();
	}

	public Alignment (int[][] A, int[] ids, Configuration c) {
		this(A, ids, false, c, 0);
	}
	
	public Alignment (int[][] A, int[] ids, Configuration c, int startIndex) {
		this(A, ids, false, c, startIndex);
	}

	
	public static void setAlphabetLength (int len) {
		alphabetLength = len;
	}
	
	protected void makeProfile () {
		
		positions = new int[K];
		posInUgappedString = new int[K][M+1];//one-based
		for (int i=0; i<K; i++)
			positions[i] = 0;
		
		
		
		//the +1 thing here shifts over the characters to a 1-based counting method
		chars = new int[M+1][];
		freqs = new int[M+1][];
		
		
		f0 = new int[M+1];
		f1 = new int[M+1];
		
		// values required for boundary conditions and terminal gap cost calculations
		firstLetterCount = new int[M+1];
		lastLetterCount = new int[M+1];
		gapsBeforeFirst = new int[M+1];
		gapsAfterLast = new int[M+1];
		
		firstLetterLoc = new int[K];
		lastLetterLoc = new int[K];

		
		f0[0] = 0;
		f1[0] = K;
		
		int[] counts = new int[alphabetLength];  	

		int[] cs;
		int[] fs;
		int k;
		int cnt1;
		int cnt00, cnt01, cnt10, cnt11;
				
		f00 = new int[M+1];
		f01 = new int[M+1];
		f10 = new int[M+1];
		f11 = new int[M+1];
		

		for (int i=0; i<M; i++) { // each column
			//initialize for new column
			int letters = 0;
			cnt1 = cnt00 = cnt01 = cnt10 = cnt11 = 0;
			
			for (int j=0; j<alphabetLength; j++)
				counts[j]=0;
			
			for (int j=0; j<K; j++) { //each row
				
				
				int c = seqs[j][i];
				int prev = i==0 ? 1 : (seqs[j][i-1]);  //just used char = 1 as a proxy, to note that some character was there
				
				if (c != SequenceConverter.GAP_VAL){ 
					positions[j]++;

					if ( counts[c] == 0 ) 
						letters++;
					counts[c]++;
					cnt1++;
					
					if (prev == SequenceConverter.GAP_VAL)
						cnt01++;
					else
						cnt11++;
				} else {
					if (prev == SequenceConverter.GAP_VAL)
						cnt00++;
					else
						cnt10++;
				}
				posInUgappedString[j][i+1] = positions[j];
			}
			
			cs = new int[letters];
			fs = new int[letters];
			k=0;
			for (int j=0; j<alphabetLength; j++){
				if (counts[j]>0) {
					cs[k] = j;
					fs[k] = counts[j];
					k++;
				}
			}
			chars[i+1] = cs;
			freqs[i+1] = fs;
			f1[i+1] = cnt1;
			f0[i+1] = K-cnt1;

			f00[i+1] = cnt00;
			f01[i+1] = cnt01;
			f10[i+1] = cnt10;
			f11[i+1] = cnt11;			

		}
		// values required for boundary conditions and terminal gap cost calculations
		firstLetterCount[0] = lastLetterCount[0] =
			gapsBeforeFirst[0] = gapsAfterLast[0] = 0;
		
		for (int i=0; i<K; i++) {
			int j=0;
			while ( seqs[i][j] == SequenceConverter.GAP_VAL ) {
				gapsBeforeFirst[++j]++; 
			}
			firstLetterCount[j+1]++;
			firstLetterLoc[i] = j+1;
			j=seqs[0].length-1;
			while ( seqs[i][j] == SequenceConverter.GAP_VAL ) {
				gapsAfterLast[j+1]++;
				j--;
			}
			lastLetterCount[j+1]++;
			lastLetterLoc[i] = j+1;
		}	
	}
	
	public static int[][] getDegappedCopy (int[][] A) {

		int K = A.length;
		int M = A[0].length;
		//count gap columns
		boolean[] goodCols = new boolean[M];
		int nonGapColCount = 0;
		for (int i=0; i<A[0].length; i++){ //columns
			goodCols[i] = false;
			for (int j=0; j<K; j++) {
				if (A[j][i] != SequenceConverter.GAP_VAL) {
					nonGapColCount++;
					goodCols[i] = true;
					break;
				}
			}
		}
		int[][] seqs = new int[K][nonGapColCount];
		M=0;
		for (int i=0; i<A[0].length; i++){ //columns
			if ( goodCols[i] ) { 
				for (int j=0; j<K; j++) {
					seqs[j][M] = A[j][i]; 
				}
				M++;
			}			
		}
		return seqs;
	}

	public static int[][] getDegappedSeqs (int[][] A) {

		int K = A.length;
		int M = A[0].length;
		int[][] seqs = new int[K][];
		
		for (int i=0; i<K; i++){ //for each seq
			int nonGapCnt = 0;
			for (int j=0; j<M; j++) { //for each column
				if (A[i][j] != SequenceConverter.GAP_VAL) {
					nonGapCnt++;
				}
			}
			int[] newSeq = new int[nonGapCnt];
			int x = 0;
			for (int j=0; j<M; j++) { //for each column
				if (A[i][j] != SequenceConverter.GAP_VAL) {
					newSeq[x++] = A[i][j];
				}
			}
			seqs[i] = newSeq;
		}
		return seqs;
	}

	
	public void countKmers () {
		if (K>1) {
			LogWriter.stdErrLogln("countKmers called for an alignment of more than one sequence. That's not right.");
			throw new GenericOpalException("countKmers called for an alignment of more than one sequence. That's not right.");
		}		

		for (int i=0; i<M-kmerK+1; i++) {
			StringBuffer sb = new StringBuffer();
			for (int j=0; j<kmerK; j++ )
				sb.append( seqs[0][i+j] + ",");
			//	sb.append( SequenceConverter.compressedVals[seqs[0][i+j]] + ","); 			
			String str = sb.toString();
			Integer x = kmers.get(str);
			kmers.put(str, 1 + (x == null ? 0 : x.intValue()) );
		}
	}

	public static Alignment buildNewAlignment (int[] seq, int seqId, Configuration c){
		return buildNewAlignment(seq, seqId, c, 0);
	}
	
	public static Alignment buildNewAlignment (int[] seq, int seqId, Configuration c, int startIndex) {
		if (c.useStructure)
			return new StructureAlignment(seq, seqId, c, startIndex);
		else
			return  new Alignment(seq, seqId, c, startIndex);
	}	


	public static Alignment buildNewAlignment (int[][] seqs, int[] seqIds, boolean isReverse, Configuration c) {
		return buildNewAlignment(seqs, seqIds, isReverse, c, 0);
	}
	
	public static Alignment buildNewAlignment (int[][] seqs, int[] seqIds, boolean isReverse, Configuration c, int startIndex) {
		if (c.useStructure)
			return new StructureAlignment(seqs, seqIds, isReverse,c, startIndex );
		else
			return  new Alignment(seqs, seqIds, isReverse, c, startIndex);
	}
	
	public static Alignment buildNewAlignment (int[][] seqs, int[] seqIds, Configuration c, int startIndex) {
		return buildNewAlignment(seqs, seqIds, false, c, startIndex);
	}	
	
	public static Alignment buildNewAlignment (int[][] seqs, int[] seqIds, Configuration c) {
		return buildNewAlignment(seqs, seqIds, false, c, 0);
	}
}

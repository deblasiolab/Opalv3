package opal.IO;

import java.util.ArrayList;
import opal.align.Aligner;
import opal.align.Alignment;
import opal.exceptions.GenericOpalException;

import com.traviswheeler.libs.LogWriter;

public class SequenceConverter {

	public static char GAP_CHAR = '-';
	public static int GAP_VAL = -2;
	public static int UNKNOWN_VAL = -1;

	char[] alphabet;
	int [] reverseLookup = new int[256];
	public int[] compressedVals = new int[256];

	public SequenceConverter (char[] chars ) {
		alphabet = chars; 
		buildReverseLookup();
	}

	public char[] getAlphabet() {
		return alphabet;
	}
	
	private void buildReverseLookup () {		
		for (int i=0; i<reverseLookup.length; i++) 
			reverseLookup[i] = UNKNOWN_VAL;

		for (int i=0; i<alphabet.length; i++) { 
			reverseLookup[Character.toLowerCase(alphabet[i])] = i;
			reverseLookup[Character.toUpperCase(alphabet[i])] = i;			
		}
	}
	
	public int[][] convertSeqsToInts (char[][] seqs) {
		int[][] converted = new int[seqs.length][]; 
		for (int i=0; i<seqs.length; i++)  {//columns
			converted[i] = new int[seqs[i].length];
			for (int j=0; j<seqs[i].length; j++) { 
				converted[i][j] = seqs[i][j] == GAP_CHAR ? GAP_VAL : reverseLookup[seqs[i][j]];
				if (converted[i][j] == UNKNOWN_VAL) {
					LogWriter.stdErrLogln("Unknown character found in string #" + i + " position " + j +" : '" + seqs[i][j] +"'");
					throw new GenericOpalException("Unknown character found in string #" + i);
				}
			}
		}
		return converted;
	}

	public char[][] convertIntsToSeqs (int[][] seqs) {
		char[][] converted = new char[seqs.length][];
		for (int i=0; i<seqs.length; i++)  {//columns
			converted[i] = new char[seqs[i].length];
			for (int j=0; j<seqs[i].length; j++) { 
				converted[i][j] = seqs[i][j] == GAP_VAL ? GAP_CHAR :  alphabet[seqs[i][j]];
				if (seqs[i][j] == UNKNOWN_VAL) {
					LogWriter.stdErrLogln("Unknown character found in string #" + i);
					throw new GenericOpalException("Unknown character found in string #" + i);
				}
			}
		}
		return converted;
	}
	
	public char[][] convertIntArrayToCharAlignment (int[][] alignment, char[][] strs /* these are required in order to recover capitalization */) {
		int len=0;
		int K=0;
		for (int i=0; i<alignment.length; i++) { 
			if (alignment[i] != null) {
				if (alignment[i].length > len) len = alignment[i].length ;
				K++;
			}
		} 
		
//		char[][] converted = new char[alignment.length][alignment[0].length];
		char[][] converted = new char[K][len];
		int k=0;
		for (int i=0; i<alignment.length; i++) { // for each sequence
			if (alignment[i] == null || alignment[i].length ==0) {
				
			} else {
			
				int x = 0;
				for (int j=0; j<alignment[i].length; j++)  {//columns
					if (alignment[i][j] == GAP_VAL) { 
						converted[k][j] =  GAP_CHAR;
					} else { 
						if (strs != null)
							converted[k][j] = strs[i][x++];
						else 
							converted[k][j] = alphabet[alignment[i][j]];
					}
				}
				k++;
			}
		}
		
		return converted;
			
	}

	public char[][] convertPathToCharAlignment (ArrayList<Aligner.Direction> path, char[][] A, char[][] B) {
		int K = A.length;
		int L = B.length;
		char[][] C = new char[K+L][path.size()];
		int i = 0, j = 0, k = 0;
		
//		iterate over vector, adding columns/dashes to C as appropriate 
		//for ( Iterator iter = path.iterator();  iter.hasNext() ; )   {
		for (Aligner.Direction direc: path) {
//		   Integer direc = (Integer)iter.next();
		   if ( Aligner.Direction.vert == direc || Aligner.Direction.diag == direc ) {
			   // copy col i of A  into first rows
			   for (int x=0; x<K; x++) 
				   C[x][k] = A[x][i];
			   i++;
		   } else {
			   // dashes into first rows
			   for (int x=0; x<K; x++) 
				   C[x][k] = GAP_CHAR;
		   }
		   if ( Aligner.Direction.horiz == direc || Aligner.Direction.diag == direc ) {
			   // copy col j of B into last rows
			   for (int x=0; x<L; x++) 
				   C[K+x][k] = B[x][j];
			   j++;
		   } else {
			   // dashes into last rows
			   for (int x=0; x<L; x++) 
				   C[K+x][k] = GAP_CHAR;
		   }
		   k++;		   
		}		
		return C;
	}

	public int[][] convertPathToIntAlignment (ArrayList<Aligner.Direction> path, Alignment A, Alignment B) {
		int K = A.K;
		int L = B.K;
		int[][] C = new int[K+L][path.size()];
		int i = 0, j = 0, k = 0;
		
//		iterate over vector, adding columns/dashes to C as appropriate 
//		for ( Iterator iter = path.iterator();  iter.hasNext() ; )   {
		for (Aligner.Direction direc : path) {
//		   Integer direc = (Integer)iter.next();
			
		   if ( Aligner.Direction.vert == direc || Aligner.Direction.diag == direc ) {
			   // copy col i of A  into first rows
			   for (int x=0; x<K; x++) 
				   C[x][k] = A.seqs[x][i];
			   i++;
		   } else {
			   // dashes into first rows
			   for (int x=0; x<K; x++) 
				   C[x][k] = GAP_VAL;
		   }
		   if ( Aligner.Direction.horiz == direc || Aligner.Direction.diag == direc ) {
			   // copy col j of B into last rows
			   for (int x=0; x<L; x++) {
				   C[K+x][k] = B.seqs[x][j];
			   }
			   j++;
		   } else {
			   // dashes into last rows
			   for (int x=0; x<L; x++) 
				   C[K+x][k] = GAP_VAL;
		   }
		   k++;		   
		}		
		return C;
	}

	public int[][] buildReverseAlignment (int[][] A){
		int K = A.length;
		int N = A[0].length;
		int[][] rev = new int[K][N];
		for (int i=0; i<K; i++){
			for (int j=0; j<N; j++){
				rev[i][j] = A[i][N-j-1];
			}
		}
		return rev;
	}
	
	public void fillCompressedAlph() {
		/*Dayhoff - Amino acid groups according to MAFFT (sextet5)
		# 0 =  A G P S T
		# 1 =  I L M V
		# 2 =  N D Q E B Z
		# 3 =  R H K
		# 4 =  F W Y
		# 5 =  C
		# 6 =  X U
		*/
		
		compressedVals[reverseLookup['a']] = compressedVals[reverseLookup['A']] = 
		compressedVals[reverseLookup['g']] = compressedVals[reverseLookup['G']] =
		compressedVals[reverseLookup['p']] = compressedVals[reverseLookup['P']] =
		compressedVals[reverseLookup['s']] = compressedVals[reverseLookup['S']] =
		compressedVals[reverseLookup['t']] = compressedVals[reverseLookup['T']] = 0;
		
		compressedVals[reverseLookup['i']] = compressedVals[reverseLookup['I']] = 
		compressedVals[reverseLookup['l']] = compressedVals[reverseLookup['L']] =
		compressedVals[reverseLookup['m']] = compressedVals[reverseLookup['M']] =
		compressedVals[reverseLookup['v']] = compressedVals[reverseLookup['V']] = 1;

		compressedVals[reverseLookup['n']] = compressedVals[reverseLookup['N']] = 
		compressedVals[reverseLookup['d']] = compressedVals[reverseLookup['D']] =
		compressedVals[reverseLookup['q']] = compressedVals[reverseLookup['Q']] =
		compressedVals[reverseLookup['e']] = compressedVals[reverseLookup['E']] =
		compressedVals[reverseLookup['b']] = compressedVals[reverseLookup['B']] =
		compressedVals[reverseLookup['z']] = compressedVals[reverseLookup['Z']] = 2;

		compressedVals[reverseLookup['r']] = compressedVals[reverseLookup['R']] = 
		compressedVals[reverseLookup['h']] = compressedVals[reverseLookup['H']] =
		compressedVals[reverseLookup['k']] = compressedVals[reverseLookup['K']] = 3;

		compressedVals[reverseLookup['f']] = compressedVals[reverseLookup['F']] = 
		compressedVals[reverseLookup['y']] = compressedVals[reverseLookup['Y']] =
		compressedVals[reverseLookup['v']] = compressedVals[reverseLookup['V']] = 4;

		compressedVals[reverseLookup['c']] = compressedVals[reverseLookup['C']] = 5;
	
		compressedVals[reverseLookup['x']] = compressedVals[reverseLookup['X']] =
		compressedVals[reverseLookup['u']] = compressedVals[reverseLookup['U']] = 6;
	}
}





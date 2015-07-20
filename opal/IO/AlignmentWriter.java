package opal.IO;

import java.util.ArrayList;

import com.traviswheeler.libs.LogWriter;

import opal.align.Aligner;

public abstract class AlignmentWriter {

	public static enum OutputType {fasta, clustal}; 
	public static OutputType outFormat = OutputType.fasta;

	
	String[] namesA;
	String[] namesB;
	char[][] alignment;
	int K, L;
	ArrayList<Aligner.Direction> path;
	boolean toUpper;
	
	protected int width;
	
	public AlignmentWriter(String[] namesA, String[] namesB, char[][] alignment, int K, int L, boolean toUpper) {	
		this.namesA = namesA;
		this.namesB = namesB;
		this.alignment = alignment;
		this.K = K;
		this.L = L;	
		this.toUpper = toUpper;
		setDefaultOutputWidth();
	}

	public AlignmentWriter(String[] names, char[][] alignment, int K, boolean toUpper) {	
		this (names, null, alignment, K, 0, toUpper);
	}

	public void setPath(ArrayList<Aligner.Direction> path) {
		this.path = path;
	}

	public void write(int width) {
		if (width > 0) 
			this.width = width;
		write();
	}
	
	abstract public void write();
	abstract protected void setDefaultOutputWidth();
}


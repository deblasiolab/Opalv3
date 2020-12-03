package opal.IO;

import java.util.ArrayList;

import java.io.PrintStream;

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
	PrintStream out;
	
	protected int width;
	
	public AlignmentWriter(String[] namesA, String[] namesB, char[][] alignment, int K, int L, boolean toUpper, PrintStream outStream) {	
		this.namesA = namesA;
		this.namesB = namesB;
		this.alignment = alignment;
		this.K = K;
		this.L = L;	
		this.toUpper = toUpper;
		setDefaultOutputWidth();
		out = outStream; 
	}

	public AlignmentWriter(String[] names, char[][] alignment, int K, boolean toUpper, PrintStream outStream) {	
		this (names, null, alignment, K, 0, toUpper, outStream);
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


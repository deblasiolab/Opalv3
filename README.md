# Opalv3

**Pleaee note, this is the development version of this application. Any bug reports are very much appreciated.**

Contents
--------
* Introduction
* Quick-start
* Alignment based on psipred-predicted protein structure
* Bypassing the shell script
* Alignment using parameter advising
* Alignment using adaptive local realignment
* Default parameters


Introduction
------------

Java >= 5 (aka v1.5) is required.

http://opal.cs.arizona.edu/

Opal is software for aligning multiple biological sequences. It can align 
both protein and DNA sequences, and expects inputs to be in fasta format.

In benchmark tests, Opal aligns sequences with good accuracy. On DNA, 
accuracy is similar to that of MAFTT and Muscle; on protein, accuracy is 
similar to MAFFT, and better than Muscle. 

Details of the algorithms used in Opal are available in the original ISMB 
paper, which should be cited in the event Opal is used:
   Wheeler, T.J. and Kececioglu, J.D.
   Multiple alignment by aligning alignments,
   Proceedings of the 15th ISCB Conference on Intelligent Systems  
   for Molecular Biology (ISMB), Bioinformatics 23, i559-i568, 2007.
   
Significant acceleration with negligible loss in accuracy was achieved
after publication of that paper, through heuristic guide-tree construction 
(similar to the approach used in MAFFT and Muscle).

Quick-start
-----------

(1) After cloning this repository, simply run `make` from the command line, 
        this will create a `bin` folder and a new `Opal.jar` file. 

(2) To build a alignment, use one of the following commands:

       ./opal unaligned_seqs.fasta > alignment.fasta
       (or)
       ./opal --in unaligned_seqs.fasta --out alignment.fasta


    To align two fixed alignments, run:
       ./opal --in alignment1.fasta --in2 alignment2.fasta


    For a full list of Opal arguments, use the --help flag:
       ./opal --help


Opal has been tested on unix and linux operating systems (including Mac OS X),
and should work on cygwin. It should also work in Windows by bypassing the 
shell script and directly calling Java.

If you receive an "out of memory" error message, you can increase the 
memory allocated to the Java VM by calling opal with the --mem flag, 
for example:
./opal --mem 2G unaligned_seqs.fasta > alignment.fasta

A more permanent increase to memory allocation can be acheived by editing 
the "mem" line in the opal shell script, for example, changing
   mem="1G"
to
   mem="2G"   



Alignment based on psipred-predicted protein structure 
------------------------------------------------------

Additional accuracy on protein sequence alignment is achieved by using a 
modified scoring scheme that incorporates protein secondary structure, as 
predicted by psipred. Details of this approach are available in the following 
paper, which should be cited if secondary-structure-based alignment is 
performed: 

   Kim, E., Wheeler, T.J., and Kececioglu, J.D.  Learning models for aligning
   protein sequence with predicted secondary structure, Proceedings of the 
   13th Conference on Research in Computational Molecular Biology (RECOMB), 
   Springer-Verlag Lecture Notes in Bioinformatics 5541: 586-605, 2009.

If you choose to use this approach follow these steps:

(1) Get psipred, e.g.
$ mkdir psipred32
$ cd psipred32
$ wget http://bioinfadmin.cs.ucl.ac.uk/downloads/psipred/psipred32.tar.gz
$ tar -xzf psipred32.tar.gz
$ cd src
$ make

(2) Open the psipred script runpsipred_single (in the psipred installation
directory), and modify the two lines:
set execdir
set datadir
to match to the absolute path of those directories

(3) Ensure that the the Opal script predict_structure.pl can find psipred.
Either:
  (a) place runpsipred_single in your PATH
  or
  (b) edit predict_structure.pl to set the path in which runpsipred_single
      is found, by editing the line:
        my $psipred_dir ; # if empty, use $PATH

(4) Opal may then be run with one additional argument: "--use_struct", as in:
  ./opal --use_struct unaligned_seqs.fasta > alignment.fasta          
  or
  ./opal --use_struct --in unaligned_seqs.fasta --out alignment.fasta



Facet Scores
------------
Opal can output an estimated accuracy value using Facet (short for the feature 
based accuracy estimator, see facet.cs.arizona.edu). To enable this functionality 
you first need to have a working version of psipred as described in the first 
3 steps in the section above (entitled "Alignment based on psipred-predicted 
protein structure"). You can then inform Opal to produce a Facet score 
using the "--facet" argument.

  ./opal --facet unaligned_seqs.fasta > alignment.fasta          
  or
  ./opal --facet --in unaligned_seqs.fasta --out alignment.fasta
  
in both cases the facet score will be printed to standard error.


Global parameter advising
-------------------------
Opal can natively preform global parameter advising using Facet 
(see facet.cs.arizona.edu). To enable this functionality you first need to have 
a working version of psipred as described in the first steps in the section above 
(entitled "Alignment based on psipred-predicted protein structure"). Enable global 
parameter advising by giving Opal a local realignment advising set using the argument: "--advising_configuration_file" 

  ./opal  --advising_configuration_file configuration.set \
      --in unaligned_seqs.fasta --out alignment.fasta

Advising sets can be found on the facet website http://facet.cs.arizona.edu.

Adaptive local realignment
--------------------------
Enable adaptive local realignment by giving Opal a local realignment advising set
using the argument: "--realignment_configuration_file" 

For example to use adaptive local realignment using the default parameter choices 
to produce the initial alignment use the following command: 

  ./opal --realignment_configuration_file configuration.set \
      --in unaligned_seqs.fasta --out alignment.fasta

To combine adaptive local realignment with global parameter advising either use the combined "--configureation_file" argument: 

  ./opal --configuration_file configuration.set \
      --in unaligned_seqs.fasta --out alignment.fasta
  or 
  ./opal --realignment_configuration_file configuration.set \
      --advising_configuration_file configuration.set \
      --in unaligned_seqs.fasta --out alignment.fasta
      
to specify seperate global and local advising set.

Advising sets can be found on the facet website http://facet.cs.arizona.edu.


Bypassing the shell script 
--------------------------

If you wish to bypass the shell script, you can run Opal.jar directly with
a command like:
 java -server -Xmx1G -jar Opal.jar unaligned_seqs.fasta > alignment.fasta

If you are using structure predictions while calling Opal in this way, you 
must provide a structure file of the sort produced by predict_structure.pl.
Suppose such a file is named file.ss; then call Opal like:

 java -server -Xmx1G -jar Opal.jar --structure_file file.ss file.fasta > output

If you take this approach because of some failing in the opal shell script, 
please e-mail me so I can improve the shell script

Source code is included in the jar. Extract the jar 
('mv opal runopal; jar xf Opal.jar'), then enter the directory opal.


Default parameters
------------------

Default model parameters (substituion and gap scores) were selected by 
the process described in Kececioglu and DeBlasio 2012.

Default DNA model parameters (substituion and gap scores) were selected by a 
process described in the 2007 ISMB paper. Those details specifically relate to 
parameters for protein sequences. A similar approach was used when choosing
parameters for DNA alignment, training on BRaliBase III 
(http://projects.binf.ku.dk/pgardner/bralibase/bralibase3/).


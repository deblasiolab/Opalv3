#!/usr/bin/perl


# Synopsis:  Given an input multi-sequence fasta-format file, run psipred_single
#            on each individual sequence, then merge the results into a file
#            with the format required by Opal 
# Incept:    TJW, Mon Jan 10 11:33:23 EST 2011 [Janelia]

use strict;
use Getopt::Long;
use Pod::Usage;
use Cwd;
use File::Temp qw/ tempfile tempdir /; # for mafft polishing routine

#=======================================================
# The next line may be edit
#=======================================================

my $psipred_dir ; # if empty, use $PATH
#my $psipred_dir = "/Volumes/Portable_2TB/ParamAdvising/psipred35" ; 

#=======================================================
# You shouldn't need to edit the script below this point
#=======================================================
my $man = 0;
my $help = 0;
my $seq_file;
my $struct_file;
my $cwd = getcwd;
my $tmp_dir = $cwd;
my $verbose = 0;
my $psipred_cmd="runpsipred_single";

my $result = GetOptions (
        'help!'    => \$help, 
        "man!"     => \$man,
        "in=s"     => \$seq_file,
        "out=s"    => \$struct_file,
        "tmpdir=s" => \$tmp_dir,
        "verbose"  => \$verbose
        ) 
        or pod2usage(1);
      
pod2usage(-verbose  => 2) if $man;
pod2usage(0) if ($help  || !$seq_file);

unless ($struct_file) {
	$struct_file = $seq_file;
    if ($struct_file =~ /(fa|fasta)$/) {
    	$struct_file =~ s/$1/ss/;
    } else {
    	$struct_file .= ".ss";
    } 
}

if ($psipred_dir) {
    $psipred_dir .= "/";	
}


my $tmpdir = tempdir( DIR => $tmp_dir , CLEANUP => 1 );


open (OUT, ">$struct_file");

$/ = ">";
open (IN, "<$seq_file");
<IN>; # throw out the first empty entry

chdir($tmpdir);
print "using tmpdir: $tmpdir\n" if $verbose;
while (my $entry = <IN>) {
    my ($name, $seq) = ($entry =~ /^(\S+?)\n(.+?)(>|$)/s);
    open (TMPOUT, ">$name.fasta");
    print TMPOUT ">$name\n$seq";
    close TMPOUT;
    do_cmd( "$psipred_dir$psipred_cmd $name.fasta 2&>1");
    #secondary structure is computed; now extract data needed by Opal
    print OUT ">$name\n";
    $/ = "\n"; # read one line at a time
    open (SS, "<$name.ss2");
    while (my $line = <SS>) {
    	if ($line =~ /(\d+)\s+(\S)\s+(\S)\s+(\d\.\d+)\s+(\d\.\d+)\s+(\d\.\d+)/) {
    	   	print OUT "$4 $5 $6\n";
    	} 
    } 
    $/ = ">"; # going back to the input sequence file, read one block at a time
}
chdir($cwd); # required to allow tmpdir to be removed

 sub do_cmd {
    my $cmd = shift;
    print "$cmd\n" if $verbose;
    return `$cmd`;
}

__END__

=head1 predict_structure.pl

Using predict_structure.pl

=head1 SYNOPSIS

predict_structure.pl --in <sequences.fasta> (others)

  Options:
    -help          brief help message
    -man           full documentation
    -in            input fasta file
    -out           output structure file (Opal input format)
    -tmpdir        directory in which temporary directory is formed    


=head1 OPTIONS

=over 8


=item B<--out <structure_file>>

File in which the Opal-formatted secondary structure is stored. If none 
is entered, the default is to use a modification of the name of the 
fasta file containing the input sequences. The new file will have the 
extension .ss. If the input file contains the extension  "fasta" or "fa",
that extension will be removed.  
(e.g. filename.fa --> filename.ss;  filename.seqs --> filename.seqs.ss) 


=item B<--tmpdir <temporary_directory>>

Directory (x) in which a temporary directory (y) will be placed. 
Temporary psipred files will be stored in directory x/y, which 
is removed at the completion of processing. If none is chosen, a 
random one will be generated in the current working directory 


=back

=head1 DESCRIPTION

B<predict_structure.pl> takes as input a sequence file in fasta format. 
For each sequence in the file, psipred_single is run. The results of these
psipred runs are merged into a single file with the format required by Opal 

=cut

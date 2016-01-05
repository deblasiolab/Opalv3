#!/usr/bin/perl
use strict;

my $folder = shift;
chomp $folder;
my @benchmarks = `ls $folder | grep -v BLOSUM | grep -v VTML`;
chomp @benchmarks;

foreach my $benchmark(@benchmarks){
#	print STDERR "wc -l $folder/$benchmark*\n";
	my @files_pre = `wc -l $folder/$benchmark* | grep " 0 "`;
	
	if(scalar(@files_pre) > 0){
		print "rm $folder/$benchmark*\n";
		print STDERR "$benchmark A\n";
		next;
	}

	foreach my $file(`ls $folder/$benchmark*`){
		chomp $file;
		my $qscore = `qscore -ref /Volumes/Portable\\ 2TB/ParamAdvising/1028_paramadvisor_data_transfer/benchmark/combined/ref/$benchmark -test $file 2>/dev/null`;
		if($qscore eq ""){
			print "rm $folder/$benchmark*\n";
			print STDERR "$benchmark B\n";
			last;
		}
	}

}

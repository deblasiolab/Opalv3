#!/usr/bin/perl
use strict;

my $folder = shift;
chomp $folder;
my @benchmarks = `ls $folder | grep -v BLOSUM | grep -v VTML`;
chomp @benchmarks;

foreach my $benchmark(@benchmarks){
	#print STDERR "ls $folder/$benchmark*pre\n";
	my @files_pre = `ls $folder/$benchmark*pre`;
	my $best_pre_value = -1;
	my $best_pre = "";
	foreach my $file_pre(@files_pre){
		my $facet_value = $file_pre;
		$facet_value = s/.*facetScore(.*).pre/\1/;
		if($facet_value > $best_pre_value){
			$best_pre_value = $facet_value;
			$best_pre = $file_pre;
		}
	}
	chomp $best_pre;

	my $file_post_prefix = $best_pre;
	$file_post_prefix =~ s/facetScore.*//;
	my $file_post = `ls $file_post_prefix* | grep -v pre`;
	chomp $file_post;

	my $qscore_default = `qscore -ref /Volumes/Portable\\ 2TB/ParamAdvising/1028_paramadvisor_data_transfer/benchmark/combined/ref/$benchmark -test realignment_removeBlankLines_window2_1.0belowWhole/$benchmark.pre`;
	my $qscore_pre = `qscore -ref /Volumes/Portable\\ 2TB/ParamAdvising/1028_paramadvisor_data_transfer/benchmark/combined/ref/$benchmark -test $best_pre`;
	my $qscore_post = `qscore -ref /Volumes/Portable\\ 2TB/ParamAdvising/1028_paramadvisor_data_transfer/benchmark/combined/ref/$benchmark -test $file_post`;
	my $qscore_best = `qscore -ref /Volumes/Portable\\ 2TB/ParamAdvising/1028_paramadvisor_data_transfer/benchmark/combined/ref/$benchmark -test $folder/$benchmark`;

	next if $qscore_default eq "";
	next if $qscore_pre eq "";
	next if $qscore_post eq "";
	next if $qscore_best eq "";
	
	$qscore_default =~ s/.*Q=(.*);TC=.*/\1/;
	$qscore_pre =~ s/.*Q=(.*);TC=.*/\1/;
	$qscore_post =~ s/.*Q=(.*);TC=.*/\1/;
	$qscore_best =~ s/.*Q=(.*);TC=.*/\1/;

	chomp $qscore_default;
	chomp $qscore_pre;
	chomp $qscore_post;
	chomp $qscore_best;

	print "$benchmark\t$qscore_default\t$qscore_pre\t$qscore_post\t$qscore_best\n";
}

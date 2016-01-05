#!/usr/bin/perl
use strict;

my $fname = shift;

open FILE, "$fname" or die("$fname: $!\n");

my @histogram;

my $countInFile;
while(<FILE>){
	if($_ =~ /output file.*pre/){
		$histogram[$countInFile]++;
		$countInFile = 0;
	}elsif($_ =~ /Realigned window/){
		$countInFile++;
	}
}
$histogram[$countInFile]++;


foreach my $i(0...scalar(@histogram)){
	print "$i\t $histogram[$i]\n";
}

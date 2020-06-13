#!/usr/bin/perl
use strict;

my $newMax = shift;
my $fname = shift;
open FILE, $fname or die("$fname: $!\n");

sub round{
	my $a = shift;
	my $f = int($a);
	if(($a-$f)<0.5){ return $f; }
	return $f+1;
}

sub min{
 my $a = shift;
 my $b = shift;
 if($a>$b){ return $b; }
 return $a;
}

sub max{
 my $a = shift;
 my $b = shift;
 if($a<$b){ return $b; }
 return $a;
}

my @array;
while(<FILE>){
	chomp $_;
	if(substr($_,0,1) ne "#"){
		my @spl = split(/\s+/,$_);
		if($spl[1] ne "A"){
			push(@array,[@spl[1...scalar(@spl)]]);
		}else{
		}
	}
}

my $min = -1* $array[0][0];
my $max = -1 *$array[0][0];
for(my $i=0;$i<scalar(@array);$i++){
	for(my $j=0;$j<scalar(@array);$j++){
		$array[$i][$j] = -1*$array[$i][$j];
		print $array[$i][$j]."\t";
		$min = min($min,$array[$i][$j]);
		$max = max($max,$array[$i][$j]);
	}
	print "\n";
}

print "\n\n";
print "min: $min\tmax: $max\n\n";

my $max = $array[0][0];
for(my $i=0;$i<scalar(@array);$i++){
        for(my $j=0;$j<scalar(@array);$j++){
		$array[$i][$j] = $array[$i][$j]-$min;
			
		$max = max($max,$array[$i][$j]);
		print $array[$i][$j]."\t";
	}
	print "\n";
}

print "\n\n";
my $scaleBy = $newMax/$max;
print "min: $min\tmax: $max\tscaleBy: $scaleBy\n\n";

my $max = $array[0][0];
for(my $i=0;$i<scalar(@array);$i++){
        print "{";
	for(my $j=0;$j<=$i;$j++){
                $array[$i][$j] *= $scaleBy;
                print round($array[$i][$j]);
		if($j!=$i){ print ","; }
        }
        print "},\n";
}


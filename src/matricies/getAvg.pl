#!/usr/bin/perl
use strict;


sub aaFrequencyInd{
        my $A = shift;
        if($A =~ /[Aa]/){return 0.0760;}
        if($A =~ /[Cc]/){return 0.0189;}
        if($A =~ /[Dd]/){return 0.0521;}
        if($A =~ /[Ee]/){return 0.0632;}
        if($A =~ /[Ff]/){return 0.0397;}
        if($A =~ /[Gg]/){return 0.0719;}
        if($A =~ /[Hh]/){return 0.0228;}
        if($A =~ /[Ii]/){return 0.0529;}
        if($A =~ /[Kk]/){return 0.0581;}
        if($A =~ /[Ll]/){return 0.0917;}
        if($A =~ /[Mm]/){return 0.0229;}
        if($A =~ /[Nn]/){return 0.0436;}
        if($A =~ /[Pp]/){return 0.0520;}
        if($A =~ /[Qq]/){return 0.0417;}
        if($A =~ /[Rr]/){return 0.0523;}
        if($A =~ /[Ss]/){return 0.0715;}
        if($A =~ /[Tt]/){return 0.0587;}
        if($A =~ /[Vv]/){return 0.0649;}
        if($A =~ /[Ww]/){return 0.0131;}
        if($A =~ /[Yy]/){return 0.0321;}

        return 0;

}

my $fname = shift;
my $order = "ARNDCQEGHILKMFPSTWYVBZX";
my @order = split(//,$order);
open FILE, "$fname" or die("$fname: $!\n");


my $sum = 0;
for(my $i=0;$i<scalar(@order);$i++){
for(my $j=0;$j<=$i;$j++){
my $a = $order[$i];
my $b = $order[$j];
$sum += (aaFrequencyInd($a) * aaFrequencyInd($b));
}}
#print "sum: $sum\n";
$sum = 0;
my $ct = 0;
foreach my $a (@order){
	my $line = <FILE>;
	chomp $line;
	my @line = split(/[{},]+/,$line);
	for(my $i=1;$i<scalar(@line);$i++){
		$sum += (aaFrequencyInd($a) * aaFrequencyInd($order[$i-1]) * $line[$i]);
		$ct++;
	}
}
print "$fname: $sum\n";#/$ct = ".($sum/$ct)."\n";

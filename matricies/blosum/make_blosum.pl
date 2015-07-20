#!/usr/bin/perl

#WRITTEN BY Travis Wheeler
#EDITED BY Dan DeBlasio 
#        13 Feb 2012

$multiplier1 = $ARGV[0] ;
$multiplier2 = $ARGV[1] ;

open FH, "<$ARGV[2]";

for (1..3) {<FH>}

$char_str = <FH>;
chomp $char_str;
@chars = split (/\s+/, $char_str);
shift @chars;

$qsum = 0;
for $i (0..$#chars) {
	$line = <FH>;
	@vals = split(/\s+/, $line);
	
	for $j (0..$i) {
		$qij[$i][$j] = $qij[$j][$i] = $vals[$j];
#		$qsum += $vals[$j];
#		$qsum += $vals[$j] if ($i != $j);
#		print "$chars[$i] $chars[$j] : $qij[$i][$j]\n";
	}
}

print "qsum = $qsum\n";

$tot_aas = 0;
for $i (0..$#chars){
	$pi[$i] = $qij[$i][$i];
	for $j (0..$#chars) {
		#$pi[$i] += $qij[$i][$j]/2 unless $i == $j;  #this is what the paper says, but it's wrong ... and also isn't what their code does (see line 377 of blosum.c)
		$pi[$i] += $qij[$i][$j] unless $i == $j;
	}
	print "pi[i] : $chars[$i] $pi[$i]\n";
	$tot_aas += $pi[$i];
}


#print "\ntotal: $tot_aas \n\n";


for $i (0..$#chars){
	for $j (0..$#chars) {
		if ( $i == $j ){
			$eij[$i][$i] = $pi[$i] * $pi[$i];
		} else {
#			$eij[$i][$j] = 2 * $pi[$i] * $pi[$j];
			$eij[$i][$j] = $pi[$i] * $pi[$j]; #this is what the paper says, but it's wrong ... and also isn't what their code does (see line 379 of blosum.c)
		}
	}
}


for $i (0..$#chars){
	for $j (0..$i) {
		$oij[$i][$j] = $oij[$j][$i] = $qij[$i][$j]/$eij[$i][$j];
	}
}

for $i (0..$#chars){
	for $j (0..$i) {
		$sij[$i][$j] = $sij[$j][$i] =  log($oij[$i][$j])/log(2);
	}
}

#if (0==1) {
# Compute the B, Z and X columns 
for $i (0..$#chars){
	#B (20) is the weighted average of N (2) and D (3)
	$tmp = ( $pi[2]*$sij[$i][2] + $pi[3]*$sij[$i][3] ) / ($pi[2] + $pi[3]);
	$sij[$i][20] = $sij[20][$i] = $tmp;
	#Z (21) is the weighted average of Q (5) and E (6)
	$tmp = ( $pi[5]*$sij[$i][5] + $pi[6]*$sij[$i][6] ) / ($pi[5] + $pi[6]);
	$sij[$i][21] = $sij[21][$i] = $tmp;
}

# value for BB is weighted average of all N x D combos
$tmp = $pi[2] * $pi[2] * $sij[2][2];
$tmp += $pi[2] * $pi[3] * $sij[2][3] * 2;
$tmp += $pi[3] * $pi[3] * $sij[3][3];
$sij[20][20] = $tmp / ($pi[2] + $pi[3]) ** 2;
# value for ZZ is weighted average of all Q x E combos 
$tmp = $pi[5] * $pi[5] * $sij[5][5];
$tmp += $pi[5] * $pi[6] * $sij[5][6] * 2;
$tmp += $pi[6] * $pi[6] * $sij[6][6];
$sij[21][21] = $tmp / ($pi[5] + $pi[6]) ** 2;

#score for ZB is weighted average of all (N,D)x(Q,E) combos
$tmp = $pi[2] * $pi[5] * $sij[2][5];
$tmp += $pi[2] * $pi[6] * $sij[2][6];
$tmp += $pi[3] * $pi[5] * $sij[3][5];
$tmp += $pi[3] * $pi[6] * $sij[3][6];
$tmp /= ( ($pi[2] + $pi[3]) * ($pi[5] + $pi[6])) ;
$tmp = $tmp ** 2;
$sij[21][20] = $sij[20][21] = $tmp;

# Compute values for unknown entry X (22) as the average
#    of row scores weighted by frequency



$xx = $tmp = 0;

for $i (0..$#chars){
	$x = 0;
	for $j (0..$#chars){
		$x += $pi[$j] * $sij[$i][$j];
		$xx += $pi[$i] * $pi[$j] * $sij[$i][$j];
		$tmp += $pi[$i] * $pi[$j];
	}
	$x /= $tot_aas;
	$sij[$j][22] = $sij[22][$i] = $x;
}


$tmp1 = $xx/$tmp;
$sij[22][22] = $tmp1;

#Now fill in (X) x (B,Z) ------------------------------*/
$x = $xx = 0;
for $i (0..$#chars){
	$x += $pi[$i] * $pi[2] * $sij[$i][2] +
		$pi[$i] * $pi[3] * $sij[$i][3];
	$xx += $pi[$i] * $pi[5] * $sij[$i][5] +
		$pi[$i] * $pi[6] * $sij[$i][6];
}
$x /= ($pi[2] + $pi[3]) * $tot_aas;
$xx /= ($pi[5] + $pi[6]) * $tot_aas;
$sij[20][22] = $sij[22][20] = $x;
$sij[21][22] = $sij[22][21] = $xx;

push @chars, qw(B Z X);
#}


print "score matrix\n-----------\n";
print join (" ", @chars);
print "\n";
$max = 0;
for $i (0..$#chars){
	for $j (0..$i) {
		$max = $sij[$i][$j] if $max < $sij[$i][$j];
	}
}
for $i (0..$#chars){
	for $j (0..$i) {
#		$sij2[$i][$j] = $sij2[$j][$i] = $sij[$i][$j] ;
		$sij2[$i][$j] = $sij2[$j][$i] = $sij[$i][$j]/$max ;
		if ($multiplier1) {
			$sij2[$i][$j] =  $multiplier1 * $sij2[$i][$j]  + .5 * ($sij2[$i][$j] <=> 0);
			$sij2[$i][$j] = $sij2[$j][$i] = int ($sij2[$i][$j]);
		}
	}
}


for $i (0..$#chars){
	for $j (0..$i) {
		if ($multiplier1) {
			printf ("%4d", $sij2[$i][$j] );
		} else {
			printf ("%0.4f ", $sij2[$i][$j] );
		}
	}
	print "\n";
}



#now turn it into a cost matrix
$max = 0;
for $i (0..$#chars){
	for $j (0..$i) {
		$max = $sij[$i][$j] if $max < $sij[$i][$j];
	}
}

print "max1: $max\n";

for $i (0..$#chars){
	for $j (0..$i) {
		$sij[$i][$j] = $sij[$j][$i] = -1 *  ($sij[$i][$j] - $max);
	}
}

#scale it
$max = 0;
for $i (0..$#chars){
	for $j (0..$i) {
		$max = $sij[$i][$j] if $max < $sij[$i][$j];
	}
}
print "max2: $max\n";

for $i (0..$#chars){
	for $j (0..$i) {
		$sij[$i][$j] = $sij[$j][$i] = $sij[$i][$j]/$max ;
		if ($multiplier2) {
			$sij[$i][$j] = $sij[$j][$i] = int ( $multiplier2 * $sij[$i][$j] + 0.5) ;
		}
	}
}


print "cost matrix\n-----------\n";
print join (" ", @chars);
print "\n";

for $i (0..$#chars){
	print "{";
	for $j (0..$i) {
		if ($multiplier2) {
			printf ("%d",  $sij[$i][$j] );
			print "," if $j != $i;
		} else {
			printf ("%0.4f ", $sij[$i][$j] );
		}
	}
	print "}";
	print "," if $i <=> $#chars;
	print "\n";
}

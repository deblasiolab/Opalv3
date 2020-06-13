#!/usr/bin/perl


my $fname = shift;
if($fname eq ""){
	open FILE, "<-";
}else{
	open FILE, $fname or die("$fname: $!\n");
}

my @sij;
while(<FILE>){
	chomp $_;
	if(substr($_,0,1) ne "#"){
		my @spl = split(/\s+/,$_);
		if($spl[1] ne "A"){
			push(@sij,[@spl[1...scalar(@spl)]]);
		}else{
		}
	}
}

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
#!/usr/bin/perl
use strict;
use warnings;

my $folder = shift;
chomp $folder;
open FILE,">$folder/META-INF/MANIFEST.MF"  or die("$folder/META-INF/MANIFEST.MF: $!");
print FILE "Manifest-Version: 1.0\n";
my $time = localtime();
print FILE "Implementation-Version: 3.2.x-dev (compiled ".$time.")\n";
print FILE "Main-Class: opal.Opal\n";

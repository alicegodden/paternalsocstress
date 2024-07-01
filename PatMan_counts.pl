# Title: miRNA/piRNA counts from PatMan
# Author: Dr. Leighton Folkes


#!/usr/bin/perl -w
use strict;
my $input = shift;
my %features;
if (!$input){
print STDERR "Sums feature counts from patman output\n";
print STDERR "$0: <patman file>\n";
exit();
}
open (IN, $input) or die "Can't open file $input\n";
while (<IN>){
if ($_=~m/^(\S+).*\t+(\S+)\((\d+)\).+/){
my $feature = $1;
my $sequence = $2;
my $count = $3;
$features{$feature}+=$count;
}
}
foreach my $element (keys %features){
my $count = $features{$element};
print "$element\t$count\n";
}

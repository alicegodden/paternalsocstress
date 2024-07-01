# Title: Perl script- Remove adaptors
# Author: Dr. Leighton Folkes

#!/usr/bin/perl
use strict;
use Bio::SeqIO;
use lib "/usr/local/Bioperl/lib/perl5/site_perl/5.8.3";
print STDERR "Removes adaptors - also removes short seqs (< 16nt)
and low complexity sequences by default\n";
my $input = shift;
my $adaptor = shift;
if (!$input){
print STDERR "No input file specified\n";
exit();
}
if (!$adaptor){
print STDERR "No adaptor sequence specified\n";
exit();
}
my %seqs;
#print STDERR "Please enter adaptor sequence to match:\n";
#my $adaptor = <STDIN>;
chomp $adaptor;
my $total_count = 0;
my $sequences = Bio::SeqIO->new(-format => 'fasta', -file =>
"$input");
while (my $seq = $sequences->next_seq ) {
$total_count++;
my $sequence = $seq->seq();
my $id = $seq->id();
 if ($sequence =~m/^(\w+)($adaptor\w*)$/){
 my $processed_seq = $1;
 my $adaptor = $2;
 my %bases;
 if (length $processed_seq >= 16){
 foreach (split('',$processed_seq)) {
$bases{$_} = 1 ;
}
my $bcount = keys %bases;
if ($bcount <3) { # low complexity sequence
next;
}
else{
$seqs{$processed_seq}++;
 }
 }
 }
}
foreach my $element (keys %seqs){
my $count = $seqs{$element};
print ">$element\($count\)\n$element\n";}

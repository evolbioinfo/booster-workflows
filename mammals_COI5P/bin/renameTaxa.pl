#!/usr/bin/env perl

use strict;
use warnings;
use Fatal;


my ($inputAlign,$outputAlign,$outputNames) = @ARGV;

my %seqs;
my %sizes;
my $head = "";
my $seq = "";
open(IN,$inputAlign);
while(<IN>){
    chomp;
    if(/^>(.*)$/){
	$head=$1;
    }else{
	$seqs{$head}  .= $_;
	$sizes{$head} += length($_);
    }
}
close(IN);

# We add a digit for the identical shortnames
my %spNames;
for my $h (keys(%seqs)){
    my $short = $h;
    $short =~ s/[: ]//g;
    $short =~ s/_//g;
    $short=substr($short,0,14);
    push @{$spNames{$short}},$h;
}

my %short;
for my $s (keys(%spNames)){
    my $size = scalar(@{$spNames{$s}});
    if($size > 1){
	my $i=1;
	for my $name (@{$spNames{$s}}){
	    $short{$name} = $s.$i;
	    $i++;
	}
    }else{
	for my $name (@{$spNames{$s}}){
	    $short{$name} = $s."0";
	}
    }
}

open(OUTPUT,">".$outputAlign);
open(OUTNAMES,">".$outputNames);
my @heads = keys(%seqs);
for my $head (@heads){
    my $h = $head;
    print OUTNAMES "$h\t".$short{$h}."\n";
    print OUTPUT ">".$short{$h}."\n";
    my $s = $seqs{$head};
    $s=~s/\*/-/g;
    print OUTPUT $s."\n";
}
close(OUTNAMES);
close(OUTPUT);

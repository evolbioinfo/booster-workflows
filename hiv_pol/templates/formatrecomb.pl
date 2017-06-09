#!/usr/bin/env perl

use strict;

my $curseq="";
my @curtype;
open(OUT,">formated.txt");
open(IN,"!{recomb}");
while(<IN>){
    chomp;
    if(/^\>([^\s]*)/){
	if($curseq ne ""){
	    print OUT $curseq." ".join(",",@curtype)."\n";
	}
	$curseq=$1;
	undef @curtype;
    }elsif($_ ne ""){
	my @cols = split(/\t/);
	push @curtype,$cols[2];
    }
}
if($curseq ne ""){
    print OUT $curseq." ".join(",",@curtype)."\n";
}
close(IN);
close(OUT);

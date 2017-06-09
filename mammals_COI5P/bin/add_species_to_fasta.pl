#!/usr/bin/env perl

use strict;
use warnings;
use Fatal;

my ($specimenOriginalFile, 
    $fastaFile) 
    = @ARGV;

my %species;
open(IN,$specimenOriginalFile);
<IN>;
while(<IN>){
    chomp;
    my @tab = split(/\t/);
    if(! ($tab[20] =~ /^\s*$/) ){
	$species{$tab[0]} = $tab[20];
    }
}
close(IN);

my $not_defined = 0;
open(IN,$fastaFile);
my $toprint=0;
while(<IN>){
    chomp;
    if(/^>COI-5P__(.*)/){
	my $processid = $1;
	
	if(defined $species{$processid}){
	    my $name = $species{$processid};
	    $name=~s/\s/_/g;
	    $_=">".$name."|".$processid;
	    $toprint = 1;
	}else{
	    $not_defined++;
	    $toprint = 0;
	}
    }
    if($toprint){
	print $_."\n";
    }
}
close(IN);

print STDERR "Not assigned Bold IDs: ".$not_defined."\n";

#!/usr/bin/perl -w

use strict;
use Fatal;

my $annotfile = shift(@ARGV);
my @subtypes = @ARGV;

my %subtypes;

open(IN,$annotfile);
while(<IN>){
    chomp;
    my @tab = split(/\s/);
    my $tax = $tab[0];
    my $type = $tab[1];

    $tax=~s/_/\./g;

    if($type=~/,/) {
	$type="recomb";
    }
    if ($type eq "A1" || $type eq "A2"){
	$type = "A";
    }
    if ($type eq "F1" || $type eq "F2"){
	$type = "F";
    }

    if ($type eq "01_AE"){
	$type = "recomb";
    }
    
    push @{$subtypes{$type}},$tax;
}
close(IN);

print "(";
my $i=0;
for my $t (keys(%subtypes)){
    if($i>0){
	print ",";
    }
    print "(";
    print join(",",@{$subtypes{$t}});
    print ")".$t."_subtype";
    $i++;
}
print ");";

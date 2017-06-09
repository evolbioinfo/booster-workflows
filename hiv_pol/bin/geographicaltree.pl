#!/usr/bin/perl -w
# Creates a tree having one branch separating hiv sequences
# from a subtype and a geographical area from the others.
use strict;
use Fatal;

my $annotfile = shift(@ARGV);
my $geocodes_str = shift(@ARGV);
my $subtypes_str = shift(@ARGV);

my @subtypeOK;
my @subtypeNotOK;

my @geocodes = split(/,/,$geocodes_str);
my @subtypes = split(/,/,$subtypes_str);

my %geocodes_h;
my %subtypes_h;

for my $v (@geocodes){
    $geocodes_h{$v}=1;
}
for my $v (@subtypes){
    $subtypes_h{$v}=1;
}

open(IN,$annotfile);
while(<IN>){
    chomp;
    my @tab = split(/\s/);
    my $tax = $tab[0];
    my $type = $tab[1];

    my $geo = "none";
    if($tax=~/_(.*)_/){
	$geo = $1;
    }
    
    $tax=~s/_/\./g;
    
    my $ok = 0;
    
    if(defined $subtypes_h{$type} && defined $geocodes_h{$geo}){
	$ok =1;
	push @subtypeOK,$tax;
    }else{
	    push @subtypeNotOK,$tax;
    }
}
close(IN);

if($#subtypeNotOK<$#subtypeOK){
    print "(";
    print join(",",@subtypeOK);
    print ",(";
    print join(",",@subtypeNotOK);
    print ")".join(".",@subtypes)."_".join(".",@geocodes)."_subtype);";
}else{
    print "(";
    print join(",",@subtypeNotOK);
    print ",(";
    print join(",",@subtypeOK);
    print ")".join(".",@subtypes)."_".join(".",@geocodes)."_subtype);";
}

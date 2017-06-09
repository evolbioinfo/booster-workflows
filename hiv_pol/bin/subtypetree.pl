#!/usr/bin/perl -w
use strict;
use Fatal;

my $annotfile = shift(@ARGV);
my @subtypes = @ARGV;

my @subtypeOK;
my @subtypeNotOK;

open(IN,$annotfile);
while(<IN>){
    chomp;
    my @tab = split(/\s/);
    my $tax = $tab[0];
    my $type = $tab[1];

    $tax=~s/_/\./g;
    
    my $ok = 0;
    for my $s (@subtypes){
	if($type eq $s){
	    $ok =1;
	    push @subtypeOK,$tax;
	}
    }
    if($ok == 0){
	    push @subtypeNotOK,$tax;
    }
}
close(IN);

if($#subtypeNotOK<$#subtypeOK){
    print "(";
    print join(",",@subtypeOK);
    print ",(";
    print join(",",@subtypeNotOK);
    print ")".join(".",@subtypes)."_subtype);";
}else{
    print "(";
    print join(",",@subtypeNotOK);
    print ",(";
    print join(",",@subtypeOK);
    print ")".join(".",@subtypes)."_subtype);";
}

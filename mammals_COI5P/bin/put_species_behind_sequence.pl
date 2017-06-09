#!/usr/bin/perl -w

use strict;
use warnings;
use Fatal;

my $species = "";
my $id = "";
while(<>){
    chomp;
    if(/^>(.*)\|(.*)/){
	($species,$id) = ($1,$2);
    }else{
	print $species."\t".$id."\t".$_."\n";
    }
}

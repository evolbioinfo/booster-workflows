#!/usr/bin/perl -w
# ARG1 : NCBI TAXONOMY  nodes.dmp
# ARG2 : specimen_Mammalia.tsv : For having a correspondance between taxid and species names
# ARG3 : Un fichier contenant les noms à garder (provenant de l'arbre de l'analyse par exemple)
# ARG4 : NCBI taxonomy: names.dmp

use strict;
use Fatal;
use warnings;

my $ncbinodeFile = shift(@ARGV);
my $specimenFile = shift(@ARGV);
my $speciesFile  = shift(@ARGV);
my $ncbinameFile = shift(@ARGV);

my %parent;
my %child;

my $nbSp = 0;

#Anoura_caudifera
#Cephalophus_maxwelli
#Galeopterus_variegates
#Hypsugo_arabicus
#Mazama_gouazoubira

my %conversion=(
    "Anoura_caudifera" => "Anoura_cultrata", 
    "Aotus_azarae_azarai"=> "Aotus_azarai_azarai",
    "Bos_taurus_indicus"=> "Bos_indicus",
    "Callosciurus_sp._1_MG2013"=>"Callosciurus_sp._1_MG-2013",
    "Cephalophus_maxwelli"=>"Cephalophus_sp._PG-2014", 
    "Cercopithecus_aethiops"=>"Cercopithecus_sp._WMS-2012",
    "Choeroniscus_sp."=>"Choeroniscus_sp._BOLD:AAC3416",
    "Cynopterus_JLE_sp._A"=>"Cynopterus_sp._A_JLE-2012",
    "Erethizon_dorsata"=>"Erethizon_dorsatum",
    "Galeopterus_variegates"=>"Galeopterus variegatus peninsulae",
    "Gerbillus_sp._1_TCB2013"=>"Gerbillus_sp._1_TCB-2013",
    "Hipposideros_cf._bicolor"=>"Hipposideros_sp._bicolor31",
    "Hipposideros_cf._larvatus"=>"Hipposideros_cf._larvatus_CMF-2010",
    "Homo_sapiens_ssp_Denisova"=>"Homo_sapiens_ssp._Denisova",
    "Hylomyscus_sp._1"=>"Hylomyscus_sp._1_TCD-2013",
    "Hylomyscus_sp._2"=>"Hylomyscus_sp._2_VN-2012",
    "Hylomyscus_sp._6"=>"Hylomyscus_sp._6_VN-2012",
    "Hypsugo_arabicus"=>"Hypsugo sp. R1", 
    "Kerivoula_cf._hardwickii"=>"Kerivoula_cf._hardwickii_CMF-2010",
    "Kerivoula_cf._lenis"=>"Kerivoula_cf._lenis_CMF-2010",
    "Mazama_gouazoubira"=>"Mazama_gouazoupira",
    "Murina_annamitica"=>"Murina_sp._B_CMF-2010",
    "Murina_fionae"=>"Murina_sp._A_MR-2012", 
    "Murina_lorelieae"=>"Murina_sp._a_JLE-2011",
    "Murina_sp."=>"Murina_sp._b_JLE-2011", 
    "Murina_walstoni"=>"Murina_sp._HZ-2008",
    "Myotis_annatessae"=>"Myotis_yanbarensis",
    "Myotis_cf._laniger"=>"Myotis_cf._laniger_CMF-2010",
    "Necromys_urichi"=>"Bolomys_urichi",
    "Oecomys_cf._rex"=>"Oecomys_cf._rex_BOLD:AAD7215",
    "Oecomys_sp._CMV2014"=>"Oecomys_sp._CMV-2014",
    "Phataginus_tetradactyla"=>"Manis_tricuspis",
    "Pipistrellus_cf._coromandra"=>"Pipistrellus_cf._coromandra_CMF-2010",
    "Plecotus_cf._strelkovi"=>"Plecotus_strelkovi",
    "Plecotus_macrobullaris_alpinus"=>"Plecotus_macrobullaris_macrobullaris",
    "Praomys_sp._A"=>"Praomys_sp._A_VN-2012",
    "Rattus_sp._abtc_43216"=>"Rattus_sp._ABTC_43216",
    "Rattus_sp._abtc_45409"=>"Rattus_sp._ABTC_45409",
    "Rhinolophus_cf._pusillus"=>"Rhinolophus_cf._pusillus_BOLD:AAB9127",
    "Rhinolophus_cf._thomasi"=>"Rhinolophus_cf._thomasi_CMF-2010",
    "Rhinolophus_yunanensis"=>"Rhinolophus_sp._1_KS-2008",
    "Rhinopithecus_bieti_2_RL2012"=>"Rhinopithecus_bieti_2_RL-2012",
    "Smutsia_gigantea"=>"Manis_sp."
    );
my %inversion;

for my $k ( keys(%conversion)){
    $inversion{$conversion{$k}} = $k;
} 

#
# On charge l'arbre taxonomique du ncbi
print STDERR "Loading NCBI taxonomy\n";
open(IN,$ncbinodeFile);
while(<IN>){
    chomp;
    my @temp = split(/\t\|\t/);
    my $tax    = $temp[0];
    my $parent = $temp[1];
    
    #print $parent." ".$tax."\n";
    if($parent ne $tax){
	$parent{$tax} = $parent;
	$child{$parent}{$tax}=1;
    }
}
close(IN);

#
# On charge la correspondance entre les id de taxons et les noms d'espèces
my %ncbinames;
my %ncbitaxid;
print STDERR "Loading NCBI names\n";
open(IN,$ncbinameFile);
while(<IN>){
    chomp;
    my @temp = split(/\t*\|\t*/);
    my $taxid = $temp[0];
    my $name  = $temp[1];
    my $type  = $temp[3];
    $name =~ s/\s/_/g;
    if($type eq "scientific name" || $type eq "synonym"){
        $ncbinames{$taxid}{$name} = 1;
	$ncbitaxid{$name} = $taxid;
    }
}
close(IN);

#
# On load les coorespondances entre espèce et tax id de bold
# les process seront utilisés après pour faire correspondre les séquences 
# de l'alignement avec les espèces
print STDERR "Loading SPECIMEN File\n";
my %species;
my %process;
my %boldSpeciesToId;
open(IN,$specimenFile);
while(<IN>){
    chomp;
    my @tab=split(/\t/);
    
    #if(! ($tab[20] =~ /^\s+$/) ){
    $tab[20]=~s/\s/_/g;
    $species{$tab[19]} = $tab[20];
    $process{$tab[0]} = $tab[20];
    $boldSpeciesToId{$tab[20]} = $tab[19];
    if($tab[19] eq "12341"){
	print STDERR $tab[19]."\t".$tab[20]."\n";
    }
    #}
}
close(IN);

#
# On load les noms de séquences
# On récupère les séquences que l'on va garder.
print STDERR "Getting Species of the Species File\n";
my $not_defined = 0;
open(IN,$speciesFile);
my $toprint=0;
my %keepSpecies;
while(<IN>){
    chomp;
    $keepSpecies{$_} = 1;
}
close(IN);

#print STDERR "Getting Taxa nodes\n";
#print join(",",@{$child{1}})."\n";
#my %taxa;
#getTaxa(1);

print STDERR "Converting Ncbi TaxId to Bold Id\n";
my %ncbiTaxIdToBoldId;
convertTaxIdToBoldId(1);

print STDERR "98661\t".$ncbiTaxIdToBoldId{98661}."\n";
print STDERR "69077\t".$ncbiTaxIdToBoldId{69077}."\n";

print STDERR "Keeping only fasta species && Cleaning Tree\n";
cleanTree(1);

print STDERR "Printing clade\n";
printTree(1);
print ";\n";

# Ecris l'arbre au format Newick sur la sortie standard
# Attention: Si noeud interne a une correspondance avec une séquence 
# (feuille donc) de bold, on cré une feuille
sub printTree{
    my ($node) = @_;
    #no children
    if(! defined $child{$node} || scalar(keys(%{$child{$node}}))==0){ #&& defined $taxa{$node}
	my $boldids = $ncbiTaxIdToBoldId{$node};
	if(defined $boldids){
	    my $f = 1;
	    # May be several boldids for on ncbi id,
	    # We test them all
	    for my $bid (keys(%{$boldids})){
		if(defined $species{$bid} && defined $keepSpecies{$species{$bid}}){
		    print "," if(!$f);
		    my $sp=$species{$bid};
		    #$sp=~s/\s/_/g;
		    print $sp;
		    $f = 0;
		}
	    }
	}else{
	    print $node;
	}
    }else{
	my $boldids = $ncbiTaxIdToBoldId{$node};
	my $f = 1;
	if(defined $boldids){
	    for my $bid (keys(%{$boldids})){
		if(defined $species{$bid} && defined $keepSpecies{$species{$bid}}){
		    print "," if(!$f);
		    my $sp=$species{$bid};
		    #$sp=~s/\s/_/g;
		    print $sp;
		    $f=0;
		}
	    }
	}
	print "," if (!$f);

	my @children = keys(%{$child{$node}});
	if($#children>0){
	    print "(";
	}
	my $first = 1;
	for my $cnode (@children){
	    if(!$first){
		print ",";
	    }
	    printTree($cnode);
	    $first = 0;
	}
	if($#children>0){
	    print ")".join(".",keys(%{$ncbinames{$node}}));
	}
    }
}

# Récupère les feuilles de l'arbre (taxons)
#sub getTaxa{
#    my ($node) = @_;
#
#    if(! defined $child{$node}){
#        $taxa{$node}=1;
#    }else{
#        my @children = keys(%{$child{$node}});
#        for my $cnode (@children){
#            getTaxa($cnode);
#        }
#    }
#    
#}

# Utilise les infos de nom d'espèce de Bold et de NCBI taxonomy (plus synonymes de NCBI)
# pour faire correspondre les identifiants de taxons de bold et de NCBI
sub convertTaxIdToBoldId{
    my ($node) = @_;
    
    if(defined $ncbinames{$node}){
	for my $name (keys(%{$ncbinames{$node}})){
	    if(defined $boldSpeciesToId{$name}){
		$ncbiTaxIdToBoldId{$node}{$boldSpeciesToId{$name}} = 1;
		#print $node."\t".$ncbiTaxIdToBoldId{$node}."\n";
	    }elsif(defined $inversion{$name} && defined $boldSpeciesToId{$inversion{$name}}){
		$ncbiTaxIdToBoldId{$node}{$boldSpeciesToId{$inversion{$name}}} = 1;
	    }
	}
    }

    if(defined $child{$node}){
        my @children = keys(%{$child{$node}});
        for my $cnode (@children){
            convertTaxIdToBoldId($cnode);
        }
    }
}

# Supprime de l'arbre les taxons qui n'ont pas d'équivalent 
# dans notre alignement
# Supprime également les noeuds internes qui restent tout seuls
sub cleanTree{
    my ($node) = @_;
   
    if(! defined $child{$node} || scalar(keys(%{$child{$node}}))==0){
	#if($node == 9606){
	#    print "No child: ".$species{}
	#}
	#my $boldid = $ncbiTaxIdToBoldId{$node};
	#if(defined $boldid){
	#    print STDERR $node."\t".$boldid."\n";
	#}
	
	#if(!defined $taxa{$node}){
	#    delete $child{$parent{$node}}{$node};
	#    delete $parent{$node};	    
	#}else{
	    if(!existSpecies($node)){
		delete $child{$parent{$node}}{$node};
		delete $parent{$node};
	    }
	#}
    }else{
	for my $c (keys(%{$child{$node}})){
	    cleanTree($c);
	}
	if(scalar(keys(%{$child{$node}})) == 0 && !existSpecies($node)){
	    delete $child{$parent{$node}}{$node};
	    delete $parent{$node};
	}
    }
}

sub existSpecies{
    my ($ncbiid) = @_;

    my $boldids = $ncbiTaxIdToBoldId{$ncbiid};
    if(!defined $boldids){
	return 0;
    }

    for  my $bid (keys(%{$boldids})){
	if(defined $species{$bid} && defined $keepSpecies{$species{$bid}}){
	    return 1;
	}
    }
    return 0;
}

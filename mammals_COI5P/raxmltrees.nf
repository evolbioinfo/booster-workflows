#!/usr/bin/env nextflow

params.datadir  = "$baseDir/../data/mammals"
params.specimen = "specimen_Mammalia.tsv.gz"
params.alignment= "Mammalia_COI5P_final_AA.aln.bz2"
params.resultDir= 'result'

specimen  = file([params.datadir, params.specimen].join(File.separator))
alignment = file([params.datadir, params.alignment].join(File.separator))
resultDir = file(params.resultDir)
alignmentDir = file([resultDir,"aligns"].join(File.separator))
treeDir = file([resultDir,"trees"].join(File.separator))

resultDir.with {
    mkdirs()
}

alignmentDir.with {
    mkdirs()
}

treeDir.with {
	mkdirs()
}

divisionChan = Channel.from([1,8,8,8,8,8,8,8,8,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64])

seedChan = Channel.from([31144,12009,13569,18804,29987,15882,8344,28523,1713,19400,6912,19245,10822,8025,15690,17179,12941,12163,6326,8070,31300,20913,23301,22351,16340,11283,14802,17187,23822,24089,30558,20564,30420,13486,26648,1397,31624,26165,17647,13495,7222,4676,22104,30726,11393,291,17438,21219,5454,24323,5613,5632,17970,28414,12444,24733,20892,29138,23614,23461,8130,15453,7911,2695,25677,14808,30782,19926,17091,21971,18887,29599,19896])

/**
  Convert the original alignment to FASTA
  And add the species from which the sequences
  come from
*/
process alignmentToFasta {
	input :
	  file specimen
	  file alignment
	output:
	  file alignFasta

	shell:
	'''
	#!/usr/bin/env bash
	add_species_to_fasta.pl <(zcat !{specimen}) <(bzcat !{alignment}) > alignFasta
	'''
}

/**
  Removes sequences corresponding to species 
  that have more than one sequence
*/
process selectOneSequenceSpecies {
	input : 
	file alignFasta
	
	output:
	set file("original.fa.gz"),file("name_map.txt.gz") into originalAlignment, originalAlignmentCopy

	shell:
	'''
	#!/usr/bin/env bash
	put_species_behind_sequence.pl !{alignFasta} | cut -f 1,3 > outalign.tmp
	select_sequences_with_one_species.Rscript outalign.tmp original_tmp.fa
	goalign trim name -i original_tmp.fa -m name_map.txt -o original.fa -n 20
	gzip original.fa
	gzip name_map.txt
	rm out_align.tmp
	rm original_tmp.fa
	'''
}

originalAlignmentCopy.subscribe{
	file,file2 -> file.copyTo(alignmentDir.resolve(file.name)); file2.copyTo(alignmentDir.resolve(file2.name))
}

/**
  Samples the original alignment using different divisions:
  1, 2, 4, 8, 16, 32, 64, 128
*/
process divideAlignment{
	tag "${align} : div ${div} - seed ${seed}"

	input:
	set file(align),file(map) from originalAlignment.first()
	val(div) from divisionChan
	val(seed) from seedChan

	output:
	set val(div), val(seed), file("align_${div}_${seed}.fa.gz") into dividedAlignment, divideAlignmentCopy

	shell:
	'''
	#!/usr/bin/env bash

	NBTAX=`goalign stats nseq -i !{align}`
	NB=$(($NBTAX/!{div}))

	if [ "!{div}" -eq "1" ]
	then
	    ln -s !{align} align_!{div}_!{seed}.fa.gz
	else
	    goalign sample seqs -i !{align} -n $NB -s !{seed} -o align_!{div}_!{seed}.fa
	    gzip align_!{div}_!{seed}.fa
	fi
	'''
}

divideAlignmentCopy.subscribe{
	div,seed, file -> file.copyTo(alignmentDir.resolve(file.name))
}

/**
  The Process that will reconstruct the Reference trees and Bootstrap Trees
*/
process runRAxML {
	tag "${refAlign} : div ${div} - seed ${seed}"

	input:
	set val(div), val(seed), file(refAlign) from dividedAlignment

	output:
	set val(div), val(seed), file("refTree.nw.gz"),file("bootTrees.nw.gz") into trees
	
	shell:
	'''
	#!/usr/bin/env bash
	goalign reformat phylip -i !{refAlign} > al.phy
	goalign reformat fasta  -i !{refAlign} > al.fa
	raxmlHPC-PTHREADS -f a -m PROTGAMMAWAG -c 6 -s al.phy -n TEST -T !{task.cpus} -p $RANDOM -x $RANDOM -# 1000
	gzip -c RAxML_bestTree.TEST > refTree.nw.gz
	gzip -c RAxML_bootstrap.TEST > bootTrees.nw.gz

	rm -f al.fa al.phy* tmp.nw RAxML_*
	'''
}

trees.subscribe{
	div, seed, ref,boot -> ref.copyTo(treeDir.resolve("ref_"+div+"_"+seed+".nw.gz")); boot.copyTo(treeDir.resolve("boot_"+div+"_"+seed+".nw.gz"));
}

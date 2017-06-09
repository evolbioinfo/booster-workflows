#!/usr/bin/env nextflow

/* "True" tree dir (tree that will serve as reference for simulations) */
params.datadir       = "$baseDir/../mammals_COI5P/result/trees/"
params.resultdir     = 'result'
params.inittree      = 'ref_phyml_1_31144.nw.gz'
params.seqlen        = 200
params.minbrlen      = 0
/* Percentage of alignment positions we will shuffle */
params.rateshuffle   = 0.5
/* Percentage of alignment we will recombine */
params.raterecombi   = 0
params.lengthrecombi = 0
params.raterogue     = 0
/*params.lengthrogue   = 0.5*/
params.simulator     = "seqgen"

params.genmodel  = 'WAG'
params.treemodel = 'WAG'

originalTree = file([params.datadir,params.inittree].join(File.separator))
resultDir    = file(params.resultdir)
treeDir      = file([params.resultdir,"trees"].join(File.separator))
alignmentDir = file([params.resultdir,"aligns"].join(File.separator))

seqlen        = params.seqlen
minbrlen      = params.minbrlen
rateShuffle   = params.rateshuffle
rateRecombi   = params.raterecombi
lengthRecombi = params.lengthrecombi
rateRogue     = params.raterogue
/*lengthRogue = params.lengthrogue*/
simulator     = params.simulator

genmodel  = params.genmodel
treemodel = params.treemodel

fasttreemodel = ''
if(treemodel.equals("WAG")){
    fasttreemodel = '-wag'
}

resultDir.with {mkdirs()}
alignmentDir.with {mkdirs()}
treeDir.with {mkdirs()}

divisionChan = Channel.from([1, 8, 8, 8, 8, 8, 8, 8, 8, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64])

seedChan = Channel.from([65512081, 67314813, 29193783, 22809492, 29601577, 64499408, 5813842, 38806899, 53759907, 2689187, 63333043, 33698375, 33513082, 90261325, 70327223, 43834756, 86136850, 22949422, 5532617, 58303852, 96496114, 92158152, 31698872, 54929374, 23454316, 71561917, 34702350, 32591600, 28339787, 41201882, 27190056, 46180751, 37998346, 27504745, 38335643, 35759018, 30618122, 40363485, 93118766, 56295590, 42010332, 20498030, 26641928, 47847232, 12192448, 60565314, 25997758, 22058783, 74692809, 62882754, 42022524, 48494987, 94042845, 29281382, 55553992, 64329656, 98284480, 53380616, 43086462, 95835734, 90288732, 34669284, 46635252, 44865659, 16560924, 52470104, 72511404, 64150685, 86344560, 29054110, 14254813, 60960293, 98508141])

process getOriginalTree{

	cpus 1
	memory '500M'
	time '1m'

	input:
	file originalTree

	output:
	file("tree.nh.gz") into firsttree

	shell:
	'''
	#!/usr/bin/env bash
	gotree minbrlen -l !{minbrlen} -i !{originalTree} \
	       | gzip -c                                  \
	       > tree.nh.gz
	'''
}

process simulateFasta{
	cpus 1

	scratch true

	cpus 1
	memory '1G'
	time '10m'

	input :
	file(tree) from firsttree

	output:
	file("original.fa.gz") into alignOriginalFasta

	shell:
	if (simulator == "seqgen")
		template 'seqgen.sh'
	else if (simulator == "indelible")
		template 'indelible.sh'
	else
		error "Invalid simulation mode: ${simulator} (must be indelible or seqgen)"
}

process shuffleOriginalFasta{
	scratch true

	cpus 1
	memory '500M'
	time '10m'

	input:
	file(original) from alignOriginalFasta

	output:
	file("shuffled.fa.gz") into alignFasta
	file("rogues.txt") into rogue

	/* Works also if rate shuffle and rate recombi are 0 */
	shell:
	'''
	#!/usr/bin/env bash
	goalign shuffle recomb -i !{original}          \
		               -n !{rateRecombi}       \
			       -l !{lengthRecombi}     \
			       -s ${RANDOM}            \
		| goalign shuffle sites                \
			       -r !{rateShuffle}       \
			       -s ${RANDOM}            \
			       --rogue !{rateRogue}    \
			       --rogue-file rogues.txt \
		| gzip -c - >  shuffled.fa.gz
	'''
}

rogue.subscribe{
	file -> file.copyTo(alignmentDir.resolve(file.name))
}

alignFasta.into{originalAlignment; originalAlignmentCopy}

originalAlignmentCopy.subscribe{
	file -> file.copyTo(alignmentDir.resolve(file.name))
}

/**
  Divides the original alignment
*/
process divideAlignment{
	tag "${align} : div ${div} - seed ${seed}"

	scratch true

	cpus 1
	memory '1G'
	time '10m'

	input:
	file(align) from originalAlignment.first()
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
	    gunzip -c !{align} > align_!{div}_!{seed}.fa
	else
	    goalign sample -i !{align} -n $NB -s !{seed} -o align_!{div}_!{seed}.fa
	fi
	gzip align_!{div}_!{seed}.fa
	'''
}

divideAlignmentCopy.subscribe{
	div,seed,file -> file.copyTo(alignmentDir.resolve(file.name))
}


/**
  The Process that will reconstruct the Reference trees and Bootstrap Trees
*/
process runRAxML {
	tag "${refAlign} : div ${div} - seed ${seed}"

	memory '20G'

	cpus 12
	scratch true

	input:
	set val(div), val(seed), file(refAlign) from dividedAlignment

	output:
	set val(div), val(seed), file("refTree.nw.gz"),file("bootTrees.nw.gz") into trees
	
	shell:
	'''
	#!/usr/bin/env bash
	goalign reformat phylip -i !{refAlign} > al.phy
	raxmlHPC-PTHREADS -f a -m PROTGAMMAWAG -c 4 -s al.phy -n TEST -T 12 -p $RANDOM -x $RANDOM -# 1000
	gzip -c RAxML_bestTree.TEST > refTree.nw.gz
	gzip -c RAxML_bootstrap.TEST > bootTrees.nw.gz

	rm -f al.phy* tmp.nw RAxML_*
	'''
}

trees.subscribe{
	div, seed, ref,boot -> ref.copyTo(treeDir.resolve("ref_"+div+"_"+seed+".nw.gz")); boot.copyTo(treeDir.resolve("boot_"+div+"_"+seed+".nw.gz"));
}

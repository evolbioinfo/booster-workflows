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

divisionChan = Channel.from([8, 8, 8, 8, 8, 8, 8, 8])

seedChan = Channel.from([67314813, 29193783, 22809492, 29601577, 64499408, 5813842, 38806899, 53759907])

process getOriginalTree{

	input:
	file originalTree

	output:
	file "tree.nh.gz" into firsttree1, firsttree2

	shell:
	'''
	#!/usr/bin/env bash
	gotree minbrlen -l !{minbrlen} -i !{originalTree} \
	       | gzip -c                                  \
	       > tree.nh.gz
	'''
}

process simulateFasta{

	input :
	file(tree) from firsttree1

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

	input:
	file(align) from originalAlignment.first()
	val(div) from divisionChan
	val(seed) from seedChan

	output:
	set val(div), val(seed), stdout, file("align_${div}_${seed}.fa.gz") into dividedAlignment, divideAlignmentCopy

	shell:
	'''
	#!/usr/bin/env bash

	NBTAX=`goalign stats nseq -i !{align}`
	NB=$(($NBTAX/!{div}))

	if [ "!{div}" -eq "1" ]
	then
	    gunzip -c !{align} > align_!{div}_!{seed}.fa
	else
	    goalign sample seqs -i !{align} -n $NB -s !{seed} -o align_!{div}_!{seed}.fa
	fi
	goalign stats length -i align_!{div}_!{seed}.fa
	gzip align_!{div}_!{seed}.fa
	'''
}

divideAlignmentCopy.subscribe{
	div,seed,length, file -> file.copyTo(alignmentDir.resolve(file.name))
}

/**
We divide the channel to the channels :
   1) The Channel that will bootstrap the alignment then construct the bootstrap trees;
   2) The channel that will reconstruct the reference tree
*/
dividedAlignment.into{dividedAlignmentToBoot; refAlignmentFastTree; dividedAlignmentToBootSimu}

/**
   The Process that will reconstruct the Reference trees
*/
process runRefFastTree {
	tag "${refAlign} : div ${div} - seed ${seed}"

	input:
	set val(div), val(seed), val(length), file(refAlign) from refAlignmentFastTree

	output:
	set val(div), val(seed), val(length), file("refTree.nw.gz") into refTreeOutput
	
	shell:
	'''
	#!/usr/bin/env bash
	goalign reformat fasta -i !{refAlign} -o tmp.fa
	goalign reformat phylip -i !{refAlign} -o tmp.phylip
	FastTree -nopr -nosupport !{fasttreemodel} -gamma tmp.fa > tmp.nw
        phyml -i tmp.phylip -b 0 -m !{treemodel} -a e -t e -o lr -u tmp.nw -d aa --quiet 1>&2
        mv tmp.phylip_phyml_tree.txt refTree.nw
	gzip refTree.nw
	rm -f tmp.fa tmp.phylip* tmp.nw
	'''
}

refTreeOutput.subscribe{
	div, seed, length, fasttree -> fasttree.copyTo(treeDir.resolve("ref_"+div+"_"+seed+".nw.gz"))
}


/**
We build the 1000 bootstraps alignment for a given reference alignment in a set of gz files
*/
process bootstrapAlignments {
	tag "${refAlign} : div ${div} - seed ${seed}"

	input:
	set val(div), val(seed), val(length), file(refAlign) from dividedAlignmentToBoot

	output:
	set val(div), val(seed), val(length), val("boot"), file("boot.tar") into bootAlignment

	shell:
	'''
	#!/usr/bin/env bash
	goalign build seqboot -n 1000 -o boot_ -s !{seed} -S -t !{task.cpus} -i !{refAlign} --gz
	tar -cf boot.tar boot_* > /dev/null
	rm -f boot_*
	'''
}


/**
We build the 1000 bootstraps alignment for a given reference alignment in a set of gz files
*/
process bootstrapSimuAlignments {
	tag "${refAlign} : div ${div} - seed ${seed}"

	input:
	set val(div), val(seed), val(length), file(refAlign) from dividedAlignmentToBootSimu
	file tree from firsttree2

	output:
	set val(div), val(seed), val(length), val("simu"), file("bootsimu.tar") into bootSimuAlignment

	shell:
	template 'indelible_boot.sh'
}


/**
	We first divide the bootstrap trees then For each bootstrap alignment, we run FastTree
*/
process divideBootAlign{
	tag "${bootFile} : div ${div} - seed ${seed}"

	input:
	set val(div), val(seed), val(length), val(type), file(bootFile) from bootAlignment.mix(bootSimuAlignment)

	output:
	set val(div), val(seed), val(type), file("boot_*.fa.gz") into bootAlignToTree mode flatten
	
	shell:
	'''
	#!/usr/bin/env bash
	tar -xf !{bootFile}
	'''
}

process runBootFastTree {
	tag "${bootFile} : div ${div} - seed ${seed}"

	module 'FastTree/2.1.8'
	module 'perl/5.22.0'

	cpus 1
	memory '5G'
	time '15h'

	scratch true

	input:
	set val(div), val(seed), val(type), file(bootFile) from dividedBootAlignToTree

	output:
	set val(div), val(seed), val(type), file("${bootFile.baseName}.nw.gz") into bootTreeOutput
	
	shell:
	'''
	#!/usr/bin/env bash
        goalign reformat fasta -i !{bootFile}                       \
		| FastTree -nopr -nosupport !{fasttreemodel} -gamma \
		| gzip -c -                                         \
		> !{bootFile.baseName}.nw.gz
	'''
}

/**
 We group the bootstrap trees per ref alignment and alignment length
 => and send that into the Channel 
 joinedBootTreeOutput
*/
joinedBootTreeOutput = bootTreeOutput.groupTuple(by : [0, 1, 2])

/**
This process reads the Channel joinedBootTreeOutput
and Concatenate all bootstrap trees into a single file
*/
process concatBootTreeFiles {

	input:
	set val(div), val(seed), val(type), file(bootstrapFileList) from joinedBootTreeOutput
	
	output:
	set val(div), val(seed), val(type), file("bootTrees.nw.gz") into concatBootTreeOutput

	shell:
	'''
	gunzip -c  !{bootstrapFileList} | gzip -c - > bootTrees.nw.gz
	'''
}


concatBootTreeOutput.subscribe{
	div, seed, length, bootTree -> bootTree.copyTo(treeDir.resolve("boot_"+div+"_"+seed+".nw.gz"));
}

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
	executor 'local'
	cpus 1
	memory '500M'
	time '1m'

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
	module 'R/3.2.3'
	executor 'local'

	cpus 1
	memory '1G'
	time '10m'

	input : 
	file alignFasta
	
	output:
	set file("original.fa.gz"),file("name_map.txt.gz") into originalAlignment, originalAlignmentCopy

	shell:
	'''
	#!/usr/bin/env bash
	put_species_behind_sequence.pl !{alignFasta} | cut -f 1,3 > outalign.tmp
	select_sequences_with_one_species.Rscript outalign.tmp original_tmp.fa
	renameTaxa.pl original_tmp.fa original.fa name_map.txt
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
	executor 'local'

	cpus 1
	memory '1G'
	time '10m'

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
	    goalign sample -i !{align} -n $NB -s !{seed} -o align_!{div}_!{seed}.fa
	    gzip align_!{div}_!{seed}.fa
	fi
	'''
}

divideAlignmentCopy.subscribe{
	div,seed, file -> file.copyTo(alignmentDir.resolve(file.name))
}

/**
  We divide the channel to the channels :
    1) The Channel that will bootstrap the alignment then construct the bootstrap trees;
    2) The channel that will reconstruct the reference tree
*/
dividedAlignment.into{dividedAlignmentToBoot; refAlignmentFastTree; refAlignmentPhyML; refAlignmentRaxML}

/**
  The Process that will reconstruct the Reference trees
*/
process runRefFastTree {
	tag "${refAlign} : div ${div} - seed ${seed}"

	module 'FastTree/2.1.8'
        module 'phyml'

	cpus 1
	memory { div == 1 ? '2G' : '1G' }

	time '15h'

	cpus 1 
	scratch true

	input:
	set val(div), val(seed), file(refAlign) from refAlignmentFastTree

	output:
	set val(div), val(seed), file("refTree.nw.gz") into refTreeOutput
	
	shell:
	'''
	#!/usr/bin/env bash
	goalign reformat phylip -i !{refAlign} > al.phy
	goalign reformat fasta  -i !{refAlign} > al.fa
	FastTree -nopr -nosupport -wag -gamma al.fa > tmp.nw
        phyml -i al.phy -b 0 -m WAG -a e -t e -o lr -u tmp.nw -d aa --quiet 1>&2
        mv al.phy_phyml_tree.txt refTree.nw
	gzip refTree.nw
	rm -f al.fa al.phy* tmp.nw
	'''
}

refTreeOutput.subscribe{
	div, seed, fasttree -> fasttree.copyTo(treeDir.resolve("ref_"+div+"_"+seed+".nw.gz"))
}


/**
  The Process that will reconstruct the Reference trees Fully with PhyML (No FastTree)
*/
process runRefPhyML {
	tag "${refAlign} : div ${div} - seed ${seed}"

	module 'phyml'

	cpus 1
        time { div <= 2 ? '10d' : ( div <= 4 ? '48h' : '10h')}
	memory { div == 1 ? '2G' : '1G' }
 

	cpus 1 
	scratch true

	input:
	set val(div), val(seed), file(refAlign) from refAlignmentPhyML

	output:
	set val(div), val(seed), file("tree.nw.gz") into refTreePhyMLOutput
	
	shell:
	'''
	#!/usr/bin/env bash
	goalign reformat phylip -i !{refAlign} > al.phy
	phyml -i al.phy -b 0 -m WAG -a e -t e -o tlr -d aa --quiet 1>&2
	mv al.phy_phyml_tree.txt tree.nw
	gzip tree.nw
	rm -f al.phy_* al.phy
	'''
}

refTreePhyMLOutput.subscribe{
	div, seed, phymltree -> phymltree.copyTo(treeDir.resolve("ref_phyml_"+div+"_"+seed+".nw.gz"))
}

/**
  The Process that will reconstruct the Reference trees with RaxML
*/
process runRefRaxML {
	tag "${refAlign} : div ${div} - seed ${seed}"

	cpus 2
	memory '1G'
        time { div <= 2 ? '10d' : ( div <= 4 ? '48h' : '10h')}
        /*clusterOptions { div <= 4 ? '--qos=long' : '--qos=normal' }*/

	cpus 1 
	scratch true

	input:
	set val(div), val(seed), file(refAlign) from refAlignmentRaxML

	output:
	set val(div), val(seed), file("tree.nw.gz") into refTreeRaxMLOutput
	
	shell:
	'''
	#!/usr/bin/env bash
	goalign reformat phylip -i !{refAlign} > al.phy
	raxmlHPC-PTHREADS -f o -p ${RANDOM} -m PROTGAMMAWAG -c 6 -s al.phy -n TEST -T 2
	mv RAxML_result.TEST tree.nw
	gzip tree.nw
	rm -f al.phy_* al.phy
	'''
}

refTreeRaxMLOutput.subscribe{
	div, seed, raxmltree -> raxmltree.copyTo(treeDir.resolve("ref_raxml_"+div+"_"+seed+".nw.gz"))
}

/**
  We build the 1000 bootstraps alignment for a given reference alignment in a set of gz files
*/
process bootstrapAlignments {
	tag "${refAlign} : div ${div} - seed ${seed}"

	scratch true

	cpus 1
	memory '1G'
	time '10m'

	input :
	set val(div), val(seed), file(refAlign) from dividedAlignmentToBoot

	output:
	set val(div), val(seed), stdout, file("boot.tar.gz") into bootAlignment

	shell:
	'''
	#!/usr/bin/env bash
	goalign build seqboot -i !{refAlign} -n 1000 -o boot -S -s ${RANDOM} --tar --gz
	printf `goalign stats nseq -i !{refAlign}`
	'''
}

dividedBootAlignment = Channel.create()
groupedBootAlignment = Channel.create()

/* If the alignment is too large, we infer bootstrap trees separately */
bootAlignment.choice( dividedBootAlignment, groupedBootAlignment ) { item -> item[3].toInteger() > 180 ? 0 : 1 }

/* Bootstrap trees are inferred in the same process */
process runGroupedBootFastTree {
	tag "${bootFile} : div ${div} - seed ${seed}"

	module 'FastTree/2.1.8'

	cpus 1
	memory '1G'
	time '24h'

	scratch true

	input:
	set val(div), val(seed), val(nbtaxa), file(bootFile) from groupedBootAlignment

	output:
	set val(div), val(seed), file("${bootFile.baseName}.nw.gz") into groupedBootTreeOutput
	
	shell:
	template 'bootstrapfastree.py'
}

/**
  Bootstrap trees will be inferred independently
  We first divide the bootstrap trees then for each bootstrap alignment, we run FastTree
*/
process divideBootAlign{
	tag "${bootFile} : div ${div} - seed ${seed}"

	module 'FastTree/2.1.8'

	cpus 1
	memory '500M'
	time '10m'

	scratch true

	input:
	set val(div), val(seed), val(nbtaxa), file(bootFile) from dividedBootAlignment

	output:
	set val(div), val(seed), file("boot*.fa.gz") into dividedBootAlignToTree mode flatten
	
	shell:
	'''
	#!/usr/bin/env bash
	tar -xzvf !{bootFile}
	gzip boot*.fa
	'''
}

/* On FastTree process per bootstrap tree */
process runBootFastTree {
	tag "${bootFile} : div ${div} - seed ${seed}"

	module 'FastTree/2.1.8'

	cpus 1
	memory '1G'
	time '15h'

	scratch true

	input:
	set val(div), val(seed), file(bootFile) from dividedBootAlignToTree

	output:
	set val(div), val(seed), file("${bootFile.baseName}.nw.gz") into bootTreeOutput
	
	shell:
	'''
	#!/usr/bin/env bash
	goalign reformat fasta -i !{bootFile} > al.fa
        FastTree -nopr -nosupport -wag -gamma al.fa | gzip -c - > !{bootFile.baseName}.nw.gz
	rm -f al.fa
	'''
}

/**
 We group the bootstrap trees per ref alignment
 => and send that into the Channel
 joinedBootTreeOutput
*/
joinedBootTreeOutput = bootTreeOutput.groupTuple(by : [0, 1])

/**
This process reads the Channel joinedBootTreeOutput
and Concatenates all bootstrap trees into a single file
*/
process concatBootTreeFiles {

	scratch true

	cpus 1
	memory '500M'
	time '10m'

	input:
	set val(div), val(seed), file(bootstrapFileList) from joinedBootTreeOutput
	
	output:
	set val(div), val(seed), file("bootTrees.nw.gz") into concatBootTreeOutput

	shell:
	'''
	gunzip -c  !{bootstrapFileList} | gzip -c - > bootTrees.nw.gz
	'''
}

allBootTreeOutput = concatBootTreeOutput.mix(groupedBootTreeOutput)

allBootTreeOutput.subscribe{
	div, seed, bootTree -> bootTree.copyTo(treeDir.resolve("boot_"+div+"_"+seed+".nw.gz"));
}

#!/usr/bin/env nextflow

/* "True" tree dir (tree that will serve as reference for simulations) */
params.datadir       = "$baseDir/../mammals_COI5P/results/fasttree/trees"
params.resultdir     = 'result'
params.inittree      = 'ref_phyml_1_31144.nw.gz'
params.seqlen        = 200
params.nboot         = 1000

params.genmodel  = 'WAG'
params.raterogue   = 0.0
params.rateshuffle = 0.0

originalTree = file([params.datadir,params.inittree].join(File.separator))
resultDir    = file(params.resultdir)
treeDir      = file([params.resultdir,"trees"].join(File.separator))

nboot     = params.nboot
seqlen    = params.seqlen
genmodel  = params.genmodel

raterogue = params.raterogue
rateshuffle = params.rateshuffle

resultDir.with {mkdirs()}
treeDir.with {mkdirs()}

divisions = Channel.from([1])
seeds = Channel.from([65512081])

// Generates subsamples of the initial tree 
// by sub sampling tips and removing the others
process gentruetrees{
	input:
	file tree from originalTree
	val seed from seeds
	val div from divisions

	output:
	set val(div), val(seed), file("sampledtree.nh.gz") into truetrees

	shell:
	'''
        NBTAX=`gotree stats -i !{tree} | cut -f 3 | tail -n 1`
        NB=$(($NBTAX/!{div}))
	gotree prune -i !{tree} --random $NB -r -s !{seed} | gzip -c - > sampledtree.nh.gz
	'''
}

// One simulated alignment per true tree
process simulatefasta{

	input :
	val genmodel
	val raterogue
	val rateshuffle
	set val(div), val(seed), file(tree) from truetrees

	output:
	set val(div), val(seed), file(tree), file("original.fa.gz") into simfasta, simfastacopy
	set val(div), val(seed), file("rogues.txt.gz") into refrogues, refrogues2

	shell:
	template 'indelible.sh'
}

simfastacopy.subscribe{
	div, seed, tree, fasta -> 
	     fasta.copyTo(treeDir.resolve("align_"+div+"_"+seed+".fa.gz")); 
	     tree.copyTo(treeDir.resolve("true_"+div+"_"+seed+".nw.gz"));
}

refrogues.subscribe{
	div, seed, rogues -> rogues.copyTo(treeDir.resolve("true_"+div+"_"+seed+"_rogues.txt.gz"));
}

/**
We divide the channel to the channels :
   1) The Channel that will bootstrap the alignment then construct the bootstrap trees;
   2) The channel that will reconstruct the reference tree
*/
simfasta.into{simfastaboot; simfastabootsimu; simfastareftree}

/**
   The Process that will reconstruct the Reference trees
*/
process inferreftree {
	tag "${refAlign} : div ${div} - seed ${seed}"

	input:
	set val(div), val(seed), file(truetree), file(align) from simfastareftree

	output:
	set val(div), val(seed), file("reftree.nw.gz") into reftrees, reftreesupport
	
	shell:
	'''
	#!/usr/bin/env bash
	goalign reformat phylip -i !{align} -o al.phy
	raxmlHPC-PTHREADS -f d -p ${RANDOM} -m PROTGAMMAWAG -c 6 -s al.phy -n TEST -T !{task.cpus}
	mv RAxML_result.TEST reftree.nw
	gzip reftree.nw
	rm -f al.phy_* al.phy
	'''
}

reftrees.subscribe{
	div, seed, reftree -> reftree.copyTo(treeDir.resolve("ref_"+div+"_"+seed+".nw.gz"))
}


/**
We build the nboot bootstraps alignment for a given reference alignment in a set of gz files
*/
process bootalignments {
	tag "${refAlign} : div ${div} - seed ${seed}"

	input:
	set val(div), val(seed), file(truetree), file(align) from simfastaboot
	val nboot

	output:
	set val(div), val(seed), val("boot"), file("boot.tar") into bootaligns, bootalignscopy

	shell:
	'''
	#!/usr/bin/env bash
	goalign build seqboot -n !{nboot} -o boot_ -s !{seed} -S -t !{task.cpus} -i !{align} --gz
	tar -cf boot.tar boot_* > /dev/null
	rm -f boot_*
	'''
}


/**
We build the nboot bootstraps alignment for a given reference alignment in a set of gz files
*/
process bootsimualignments {
	tag "${refAlign} : div ${div} - seed ${seed}"

	input:
	set val(div), val(seed), file(truetree), file(align) from simfastabootsimu
	val raterogue
	val rateshuffle
	val nboot

	output:
	set val(div), val(seed), val("simu"), file("bootsimu.tar") into bootsimualigns, bootsimualignscopy
	set val(div), val(seed), file("rogues.txt.gz") into bootrogues

	shell:
	template 'indelible_boot.sh'
}


bootalignscopy.mix(bootsimualignscopy).subscribe{
        div, seed, type, aligns -> 
	     aligns.copyTo(treeDir.resolve("boot_"+type+"_"+div+"_"+seed+".tar"));
}

bootrogues.subscribe{
	     div, seed, rogues -> rogues.copyTo(treeDir.resolve("boot_simu_"+div+"_"+seed+"_rogues.txt.gz"));
}

/**
We first divide the bootstrap trees then For each bootstrap alignment, we run FastTree
*/
process divideboots{
	tag "${bootFile} : div ${div} - seed ${seed}"

	input:
	set val(div), val(seed), val(type), file(bootFile) from bootaligns.mix(bootsimualigns)

	output:
	set val(div), val(seed), val(type), file("boot*.fa.gz") into divbootaligns mode flatten
	
	shell:
	'''
	#!/usr/bin/env bash
	tar -xf !{bootFile}
	'''
}

process inferboottrees {
	tag "${bootFile} : div ${div} - seed ${seed}"

	input:
	set val(div), val(seed), val(type), file(bootFile) from divbootaligns

	output:
	set val(div), val(seed), val(type), file("${bootFile.baseName}.nw.gz") into bootTreeOutput, bootTreeOutputCopy
	
	shell:
	'''
	#!/usr/bin/env bash
	goalign reformat phylip -i !{bootFile} -o al.phy
	raxmlHPC-PTHREADS -f d -p ${RANDOM} -m PROTGAMMAWAG -c 6 -s al.phy -n TEST -T !{task.cpus}
	mv RAxML_result.TEST !{bootFile.baseName}.nw
	gzip !{bootFile.baseName}.nw
	rm -f al.phy_* al.phy
	'''
}

/**
 We group the bootstrap trees per ref alignment
 => and send that into the Channel 
 joinedBootTreeOutput
*/
joinedBootTreeOutput = bootTreeOutput.groupTuple(by : [0, 1, 2])

/**
This process reads the Channel joinedBootTreeOutput
and Concatenate all bootstrap trees into a single file
*/
process concatboottrees {

	input:
	set val(div), val(seed), val(type), file(bootstrapFileList) from joinedBootTreeOutput
	
	output:
	set val(div), val(seed), val(type), file("bootTrees.nw.gz") into concatBootTreeOutput, concatBootTreeOutputCopy

	shell:
	'''
	gunzip -c  !{bootstrapFileList} | gzip -c - > bootTrees.nw.gz
	'''
}

concatBootTreeOutputCopy.subscribe{
	div, seed, type, bootTree -> bootTree.copyTo(treeDir.resolve("boot_"+div+"_"+type+"_"+seed+".nw.gz"));
}

concatBootTreeOutput.combine(reftreesupport, by: [0,1]).into{tosupport; tosampsupport}

process support {
	input:
	set val(div), val(seed), val(type), file(boottrees), file(reftree) from tosupport

	output:
	set val(div), val(seed), val(type), file("fbp.nw"), file("tbe.nw"), file("tbe.log") into support, support2

	shell:
	'''
	gotree compute support classical -i !{reftree} -b !{boottrees} -o fbp.nw -t !{task.cpus}
	gotree compute support   booster -i !{reftree} -b !{boottrees} -o tbe.nw -t !{task.cpus} -l tbe.log --moved-taxa --dist-cutoff 0.3
	'''
}

support.subscribe{
	div, seed, type, fbp, tbe,log -> fbp.copyTo(treeDir.resolve("boot_"+type+"_"+div+"_"+seed+"_fbp.nw"));
	     	   	    	     	    tbe.copyTo(treeDir.resolve("boot_"+type+"_"+div+"_"+seed+"_tbe.nw"));
	     	   	    	     	    log.copyTo(treeDir.resolve("boot_"+type+"_"+div+"_"+seed+"_tbe.log"));
}

process supportsample {
	input:
	set val(div), val(seed), val(type), file(boottrees), file(reftree) from tosampsupport
	val nboot
	each s from 10,11,12,13

	output:
	set val(div), val(seed), val(type), val(s), file("fbp_samp.nw"), file("tbe_samp.nw") into supportsamp


	shell:
	'''
	gotree sample --replace -n !{nboot} -i !{boottrees} -s !{s} | gzip -c > boot_samp.nw.gz
	gotree compute support classical -i !{reftree} -b boot_samp.nw.gz -o fbp_samp.nw -t !{task.cpus}
	gotree compute support   booster -i !{reftree} -b boot_samp.nw.gz -o tbe_samp.nw -t !{task.cpus}
	'''
}

supportsamp.subscribe{
	div, seed, type, s, fbp, tbe -> fbp.copyTo(treeDir.resolve("boot_"+type+"_"+div+"_"+seed+"_samp"+s+"_fbp.nw"));
	     	   	    	     	tbe.copyTo(treeDir.resolve("boot_"+type+"_"+div+"_"+seed+"_samp"+s+"_tbe.nw"));
}


simu = Channel.create()
boot = Channel.create()
support2.choice( simu, boot ) { item -> item[2] == "simu" ? 0 : 1 }
supportmerge = simu.combine(boot,by:[0,1]).combine(refrogues2, by: [0,1])

process toSupportData{

	input:
	set val(div), val(seed), val(typesimu), file(fbpsimu:'fbpsimu.nw'), file(tbesimu:'tbesimu.nw'), file(tbelogsimu:'tbesimu.log'),val(typeboot), file(fbpboot:'fbpboot.nw'), file(tbeboot:'tbeboot.nw'), file(tbelogboot:'tbeboot.log'), file(refrogues) from supportmerge

	output:
	set val(div), val(seed), file("fbp.dat"), file("tbe.dat"), file(tbelogsimu), file(tbelogboot), file(refrogues) into supportdata


	shell:
	'''
	paste <(gotree stats edges -i !{tbeboot} | cut -f 1,2,3,4,7) <(gotree stats edges -i !{tbesimu} | cut -f 4) > tbe.dat
	paste <(gotree stats edges -i !{fbpboot} | cut -f 1,2,3,4,7) <(gotree stats edges -i !{fbpsimu} | cut -f 4) > fbp.dat
	'''		
}

process plotCompareFBPTBE {

	input:
	set val(div), val(seed), file(fbp), file(tbe), file(tbelogsimu), file(tbelogboot), file(refrogues) from supportdata

	output:
	file "*.svg" into compareplots mode flatten
	file "*.txt" into comparetexts mode flatten


	shell:
	template 'plot.R'
}

compareplots.subscribe{
	f -> f.copyTo(treeDir.resolve(f.name));
}
comparetexts.subscribe{
	f -> f.copyTo(treeDir.resolve(f.name));
}

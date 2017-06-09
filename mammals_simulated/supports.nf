#!/usr/bin/env nextflow

params.datadir  = "$baseDir/result/"
params.resultdir= 'result'
params.truetree = "$baseDir/../mammals_COI5P/result/trees/ref_phyml_1_31144.nw.gz"
params.roguefile = "$baseDir/../mammals_COI5P/result/align/rogues.txt"

/* Input files & directories */
alignmentDir = file([params.datadir,"aligns"].join(File.separator))
treeDir      = file([params.datadir,"trees"].join(File.separator))
trueTree     = file(params.truetree)
roguefile    = file(params.roguefile)

/* Output directories */
resultDir    = file(params.resultdir)
outtreeDir   = file([params.resultdir, "supports"].join(File.separator))
resultDir.with {mkdirs()}
outtreeDir.with {mkdirs()}

divisionChan = Channel.from([1, 8, 8, 8, 8, 8, 8, 8, 8, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64])

seedChan = Channel.from([65512081, 67314813, 29193783, 22809492, 29601577, 64499408, 5813842, 38806899, 53759907, 2689187, 63333043, 33698375, 33513082, 90261325, 70327223, 43834756, 86136850, 22949422, 5532617, 58303852, 96496114, 92158152, 31698872, 54929374, 23454316, 71561917, 34702350, 32591600, 28339787, 41201882, 27190056, 46180751, 37998346, 27504745, 38335643, 35759018, 30618122, 40363485, 93118766, 56295590, 42010332, 20498030, 26641928, 47847232, 12192448, 60565314, 25997758, 22058783, 74692809, 62882754, 42022524, 48494987, 94042845, 29281382, 55553992, 64329656, 98284480, 53380616, 43086462, 95835734, 90288732, 34669284, 46635252, 44865659, 16560924, 52470104, 72511404, 64150685, 86344560, 29054110, 14254813, 60960293, 98508141])

process getFiles{
	cpus 1
	memory '500M'
	time '30s'

	input:
	val(div) from divisionChan
	val(seed) from seedChan

	output:
	set val(div), val(seed), file(refAlign), file(refTree), file(bootTrees) into originalFiles

	shell:
	'''
	cp !{alignmentDir}/align_!{div}_!{seed}.fa.gz refAlign
	cp !{treeDir}/ref_!{div}_!{seed}.nw.gz refTree
	cp !{treeDir}/boot_!{div}_!{seed}.nw.gz bootTrees
	'''
}

originalFiles.into{originalFilesSupport; originalFilesSupportTRUE}

/************************/
/*                      */     
/*  Analyze Supports    */
/*                      */
/************************/
process computeSupports {
	tag "${refAlign.name} - ${div} ${seed}"
	scratch true

	cpus 10
	memory '1G'
	time '24h'

	input:
	set val(div), val(seed), file(refAlign), file(refTree), file(bootTrees) from originalFilesSupport
	file(roguefile)

	output:
	set val(div), val(seed), file("fbp.nw"), file("tbe.nw"), file("tbe.log"), file("fbp_norogue.nw") into supporttrees
	/* We output only the collapsed version of the trees with supports */

	shell:
	'''
	gunzip -c !{bootTrees} > boot.nw
	gunzip -c !{refTree}   > ref.nw
	gunzip -c !{refAlign}  > ali.fa	

	gotree compute support classical -i ref.nw -b boot.nw -o fbp.nw -t 5
	gotree compute support booster -i ref.nw -o boot.nw -t 10 -l tbe.log --moving-taxa --dist-cutoff 0.3 -o tbe.nw

	# We recompute FBP using trees without rogues
	gotree prune -i ref.nw -f !{roguefile} -o ref_norogue.nw
	gotree prune -i boot.nw -f !{roguefile} -o boot_norogue.nw
	gotree compute support classical -i ref_norogue.nw -b boot_norogue.nw -t 4 -o fbp_norogue.nw

	rm -f boot.nw ref.nw ali.fa
	'''
}

supporttrees.into{ supporttreescopy; supporttreesnext}
supporttreescopy.subscribe{
	div,seed,fbp, tbe, tbelog, fbpnorogue -> 
		fbp.copyTo(outtreeDir.resolve(div+"_"+seed+"_"+fbp.name));
		tbe.copyTo(outtreeDir.resolve(div+"_"+seed+"_"+tbe.name));
		tbelog.copyTo(outtreeDir.resolve(div+"_"+seed+"_"+tbelog.name));
		fbpnorogue.copyTo(outtreeDir.resolve(div+"_"+seed+"_"+fbpnorogue.name));
}

process analyzeSupports {
	tag "${div} ${seed}"

	scratch true

	cpus 1
	memory '1G'
	time '10m'

	input:
	set val(div), val(seed), file(fbp), file(tbe), file(tbelog), file(fbpnorogue) from supporttreesnext

	output:
	set val(div), val(seed), file(fbp), file(tbe), file("supports.txt") into supportstats
	
	shell:
	'''
	gotree stats edges -i !{tbe}      | awk '{if($4 != "N/A" && NR>1){print "!{div}\t!{seed}\tnot_collapsed\tTBE\t" $4 "\t" $7}}' > supports.txt
	gotree stats edges -i !{fbp} | awk '{if($4 != "N/A" && NR>1){print "!{div}\t!{seed}\tnot_collapsed\tFBP\t" $4 "\t" $7}}' >> supports.txt
	'''
}

supportstats.map{
	div, seed, fbp, tbe, supports -> supports
}.collectFile(name: 'allsupports.txt').subscribe{
	file -> file.copyTo(resultDir.resolve(file.name))
}

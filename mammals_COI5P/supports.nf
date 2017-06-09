#!/usr/bin/env nextflow

params.datadir  = "$baseDir/results/fasttree"
params.resultdir= 'results/fasttree'

/* Input files & directories */
alignmentDir     = file([params.datadir,"aligns"].join(File.separator))
treeDir          = file([params.datadir,"trees"].join(File.separator))
name_map         = file([alignmentDir, "name_map.txt.gz"].join(File.separator))

/* Output directories */
resultDir    = file(params.resultdir)
outtreeDir   = file([params.resultdir, "support"].join(File.separator))

resultDir.with {
    mkdirs()
}

outtreeDir.with {
    mkdirs()
}

divisionChan = Channel.from([1,8,8,8,8,8,8,8,8,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64])

seedChan = Channel.from([31144,12009,13569,18804,29987,15882,8344,28523,1713,19400,6912,19245,10822,8025,15690,17179,12941,12163,6326,8070,31300,20913,23301,22351,16340,11283,14802,17187,23822,24089,30558,20564,30420,13486,26648,1397,31624,26165,17647,13495,7222,4676,22104,30726,11393,291,17438,21219,5454,24323,5613,5632,17970,28414,12444,24733,20892,29138,23614,23461,8130,15453,7911,2695,25677,14808,30782,19926,17091,21971,18887,29599,19896])

process getFiles{
	cpus 1
	memory '500 MB'
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

/************************/
/*                      */     
/*  Analyze Supports    */
/*                      */
/************************/
process computeSupports {
	tag "${refAlign.name} - ${div} ${seed}"
	scratch true

	cpus 4
	memory '1 GB'
	time '24h'

	input:
	set val(div), val(seed), file(refAlign), file(refTree), file(bootTrees) from originalFiles

	output:
	set val(div), val(seed), file("fbp.nw"), file("tbe.nw"),file("fbp_norogue.nw") into supporttrees
	file("rogue_*.txt") into rogues

	shell:
	'''
	gunzip -c !{bootTrees} > boot.nw
	gunzip -c !{refTree}   > ref.nw
	gunzip -c !{refAlign}  > ali.fa	

	gotree compute support classical -i ref.nw -b boot.nw -o fbp.nw 
	booster  -i ref.nw -b boot.nw -o tbe.nw -s $RANDOM -@ 4 -n auto

	# We compute Rogue Taxa with RAxML
	gotree resolve -i boot.nw > boot_resolve.nw
	raxmlHPC-PTHREADS -J MR_DROP -z boot_resolve.nw -m GTRCAT -n boot  -m PROTGAMMAWAG -c 6 -T 4
	grep "dropping" RAxML_info.boot | cut -f 2 -d ":" | tr '\n' ' ' | tr ',' '\n' | sed 's/ //g' > rogue_!{div}_!{seed}.txt

	# We recompute FBP using trees without rogues
	gotree prune -i ref.nw -f rogue_!{div}_!{seed}.txt -o ref_norogue.nw
	gotree prune -i boot.nw -f rogue_!{div}_!{seed}.txt -o boot_norogue.nw
	gotree compute support classical -i ref_norogue.nw -b boot_norogue.nw -t 4 -o fbp_norogue.nw

	rm -f boot.nw ref.nw ali.fa boot_resolve.nw ref_norogue.nw boot_norogue.nw 
	'''
}

rogues.subscribe{
	f ->  f.copyTo(outtreeDir.resolve(f.name));
}

supporttrees.into{ supporttreescopy; supporttreesnext}
supporttreescopy.subscribe{
	div,seed,fbp,tbe,fbpnorogue -> 
		fbp.copyTo(outtreeDir.resolve(div+"_"+seed+"_"+fbp.name));
		tbe.copyTo(outtreeDir.resolve(div+"_"+seed+"_"+tbe.name));
		fbpnorogue.copyTo(outtreeDir.resolve(div+"_"+seed+"_"+fbpnorogue.name));
}

process analyzeSupports {
	tag "${div} ${seed}"

	scratch true

	cpus 1
	memory '1 GB'
	time '10m'

	input:
	set val(div), val(seed), file(fbp), file(tbe), file(fbpnorogue) from supporttreesnext

	output:
	file("supports.txt") into supportstats
	
	shell:
	'''
	# We print $4=support | $7=depth
	gotree stats edges -i !{tbe} | awk '{if($4 != "N/A" && NR>1){print "!{div}\t!{seed}\tnot_collapsed\tTBE\t" $4 "\t" $7}}' > supports.txt
	gotree stats edges -i !{fbp} | awk '{if($4 != "N/A" && NR>1){print "!{div}\t!{seed}\tnot_collapsed\tFBP\t" $4 "\t" $7}}' >> supports.txt
	'''
}

supportstats.collectFile(name: 'allsupports.txt').subscribe{
	file -> file.copyTo(resultDir.resolve(file.name))
}

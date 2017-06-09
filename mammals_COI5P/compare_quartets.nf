#!/usr/bin/env nextflow

params.div1tree  = "$baseDir/results/fasttree/trees/ref_1_31144.nw.gz"
params.truetree  = "$baseDir/../data/ncbitax/ncbi.nw"
params.boottrees = "$baseDir/results/fasttree/trees/boot_1_31144.nw.gz"
params.mapfile   = "$baseDir/results/fasttree/aligns/name_map.txt.gz"
params.resultdir= 'results/fasttree/quartets'

/* Output directories */
resultDir    = file(params.resultdir)

resultDir.with {mkdirs()}

truetree = Channel.from(file(params.truetree))
div1tree = Channel.from(file(params.div1tree))
boottrees= Channel.from(file(params.boottrees))
mapfile  = Channel.from(file(params.mapfile))

/*
This will rename taxa of ncbi into renamed taxa of the analysis (using mapfile)
And will prune ncbi+inferred+bootstrap trees to get the same set of taxa
*/
process preprocessTree {
	cpus 1
	memory '1G'
	time '10m'

	input:
	file div1tree 
	file truetree
	file boottrees
	file mapfile

	output:
	file("div1.nw") into unzipeddiv1tree
	file("true.nw") into unzipedtruetree
	file("boot.nw.gz") into preprocessedboottrees

	shell:
	'''
	#!/bin/bash
	if [ !{mapfile} != "null" ]
	then
		gotree rename -i !{truetree} -m !{mapfile} -o true_full.nw
	else
		gunzip -c !{truetree} > true_full.nw
	fi
	gotree prune -i !{div1tree}  -c true_full.nw > div1.nw
	gotree prune -i true_full.nw -c !{div1tree}  > true.nw
	gotree prune -i !{boottrees} -c true_full.nw | gzip -c > boot.nw.gz
	'''
}

unzipeddiv1tree.into{div1tree1; div1tree2; div1tree3; div1tree4; div1tree5}
unzipedtruetree.into{truetree1; truetree2; truetree3}
preprocessedboottrees.into{boottrees1; boottrees2}

/* We divide the inferred tree in as many trees as there are edges in the tree */
/* Each tree will have only one resolved branch */
process divideTreeIntoEdges {
	cpus 5
	memory '5G'
	time '10m'

	input:
	file(divtree) from div1tree1

	output:
	file("edge*") into staredges mode flatten

	shell:
	'''
	gotree compute edgetrees -i !{divtree} -o edge -t 10
	'''
}

/* For each edge tree, we compare the quartets with the NCBI tree */
process analyzeStarEdges {
	tag {staredge.baseName}
	cpus 1
	memory '1G'
	/*time '2m'*/

	input:
	file(staredge) from staredges
	file(truetree) from truetree1.first()

	output:
	file("quartet") into quartets

	shell:
	template 'quartet.sh'
}

quartets.collectFile(name: 'quartets.txt').into{quartetscopy; quartetsanalyze}
quartetscopy.subscribe{
	file -> file.copyTo(resultDir.resolve(file.name))
}

/* 
  Now we compare the transfer support of the inferred tree 
  compared to the true tree 
*/
process transferDistTrue {
	cpus 5
	memory '5G'
	time '100m'
	/*cache false*/

	input:
	file(divtree) from div1tree2
	file(truetree) from truetree2.first()

	output:
	file("transfer_to_true.txt")   into truetrans
	file("div_1_transfer_true.nw") into truetranstree

	shell:
	'''
	#!/bin/bash
	booster -i !{divtree} -b !{truetree} -n empirical -o div_1_transfer_true.nw -@ 5
	gotree stats edges -i div_1_transfer_true.nw > transfer_to_true.txt
	'''
}

truetranstree.subscribe{
	file -> file.copyTo(resultDir.resolve(file.name))
}

/* Which bipartitions are really true? */
process classicalDistTrue {
	cpus 1
	memory '5G'
	time '100m'

	input:
	file(divtree) from div1tree4
	file(truetree) from truetree3.first()

	output:
	file("class_to_true.txt")   into trueclass
	file("div_1_class_true.nw") into trueclasstree

	shell:
	'''
	#!/bin/bash
	gotree compute support classical -i !{divtree} -b !{truetree} -o div_1_class_true.nw
	gotree stats edges -i div_1_class_true.nw > class_to_true.txt
	'''
}

trueclasstree.subscribe{
	file -> file.copyTo(resultDir.resolve(file.name))
}

/* Compute transfer support of branches of inferred tree */
process transferSupport {
	cpus 5
	memory '1G'
	time '3h'

	input:
	file(divtree) from div1tree3.first()
	file(boottrees1)

	output:
	file("transfer_values.txt") into transvalues
	file("div_1_transfer.nw") into transtree

	shell:
	'''
	#!/bin/bash
	gunzip -c !{boottrees1} > boot.nw
	booster -i !{divtree} -b boot.nw -@ 5 -o div_1_transfer.nw
	gotree stats edges -i div_1_transfer.nw > transfer_values.txt
	'''
}

transtree.subscribe{
	file -> file.copyTo(resultDir.resolve(file.name))
}

/* Classical support to bootstrap trees */
process classicalSupport {
	cpus 1
	memory '5G'
	time '1h'

	input:
	file(divtree) from div1tree5
	file(boottrees) from boottrees2

	output:
	file("class_values.txt") into classvalues
	file("div_1_class.nw") into classtree

	shell:
	'''
	#!/bin/bash
	gotree compute support classical -i !{divtree} -b !{boottrees} -t 1 > div_1_class.nw
	gotree stats edges -i div_1_class.nw > class_values.txt
	'''
}

classtree.subscribe{
	file -> file.copyTo(resultDir.resolve(file.name))
}

/* We group all these information into one single file */
process groupInfos{
	input:
	file(quartetsanalyze)
	file(truetrans)
	file(trueclass)
	file(transvalues)
	file(classvalues)

	output:
	file("groupedinfos.txt") into allinfos

	shell:
	'''
	#!/bin/bash
	# Not numeric sort: the first column is of the form edge_000000.nw (0 padded)
	sort -k 1 !{quartetsanalyze} | cut -d ' ' -f 2,3,4,5,6 > quarts
	# Other files are already sorted by edge id (corresponding to the same edge than the sorted file)
	awk '{if($5=="false"){print $2 " " $3 " " $7 " " $4}}' !{truetrans} > mtrue
	awk '{if($5=="false"){print $4}}' !{trueclass}> ctrue
	awk '{if($5=="false"){print $4}}' !{classvalues}> cval
	awk '{if($5=="false"){print $4}}' !{transvalues}> mval

	paste mtrue quarts ctrue mval cval > groupedinfos.txt
	'''
}

allinfos.subscribe{
	file -> file.copyTo(resultDir.resolve(file.name))
}

#!/usr/bin/env nextflow

params.div1tree  = "$baseDir/result/trees/ref_1_65512081.nw.gz"
params.truetree  = "$baseDir/../mammals_COI5P/result/trees/ref_phyml_1_31144.nw.gz"
params.boottrees = "$baseDir/result/trees/boot_1_65512081.nw.gz"
params.roguefile = "$baseDir/result/aligns/rogues.txt"

params.mapfile   = "null"
params.resultdir= 'result/quartets'

/* Output directories */
resultDir    = file(params.resultdir)
resultDir.with {mkdirs()}

truetree = Channel.from(file(params.truetree))
div1tree = Channel.from(file(params.div1tree))
boottrees= Channel.from(file(params.boottrees))
mapfile  = Channel.from(file(params.mapfile))
roguefile= Channel.from(file(params.roguefile))

process preprocessTree {
	input:
	file div1tree 
	file truetree
	file boottrees
	file mapfile
	file roguefile

	output:
	file "div1.nw" into unzipeddiv1tree
	file "true.nw" into unzipedtruetree
	file "boot.nw.gz" into preprocessedboottrees

	shell:
	'''
	#!/bin/bash
	if [ "!{mapfile}" != "null" ]
	then
		gotree rename -i !{truetree} -m !{mapfile} -o true_full.nw
	else
		gunzip -c !{truetree} > true_full.nw
	fi
 
	# PRUNE SIMULATED ADDED ROGUE TAXA (or not if no file)
	gotree prune -i !{div1tree}  -c true_full.nw | gotree prune -f !{roguefile} > div1.nw
	gotree prune -i true_full.nw -c !{div1tree}  | gotree prune -f !{roguefile} > true.nw
	gotree prune -i !{boottrees} -c true_full.nw | gotree prune -f !{roguefile} | gzip -c > boot.nw.gz
	'''
}

unzipeddiv1tree.into{div1tree1; div1tree2; div1tree3; div1tree4; div1tree5}
unzipedtruetree.into{truetree1; truetree2; truetree3}
preprocessedboottrees.into{boottrees1; boottrees2}

/* We divide the inferred tree in as many trees as there are edges in the tree */
/* Each tree will have only one resolved branch */
process divideTreeIntoEdges {
	input:
	file divtree from div1tree1

	output:
	file "edge*" into staredges mode flatten

	shell:
	'''
	gotree compute edgetrees -i !{divtree} -o edge -t !{task.cpus}
	'''
}

/* For each edge tree, we compare the quartets with the true tree */
process analyzeStarEdges {
	tag {staredge.baseName}
	input:
	file staredge from staredges
	file truetree from truetree1.first()

	output:
	file "quartet" into quartets

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
	input:
	file divtree from div1tree2
	file truetree from truetree2.first()

	output:
	file "transfer_to_true.txt"   into truetrans
	file "div_1_transfer_true.nw" into truetranstree

	shell:
	'''
	#!/bin/bash
	booster -i !{divtree} -b !{truetree} -n empirical -o div_1_transfer_true.nw -@ !{task.cpus}
	gotree stats edges -i div_1_transfer_true.nw > transfer_to_true.txt
	'''
}

truetranstree.subscribe{
	file -> file.copyTo(resultDir.resolve(file.name))
}

/* Which bipartitions are really true? */
process classicalDistTrue {
	input:
	file divtree from div1tree4
	file truetree from truetree3.first()

	output:
	file "class_to_true.txt"   into trueclass
	file "div_1_class_true.nw" into trueclasstree

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

/* Compute transfer distance of branches of inferred tree using rand trees and one boot tree */
process transferSupport {
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
	booster -i !{divtree} -b boot.nw -@ !{task.cpus} -o div_1_transfer.nw
	gotree stats edges -i div_1_transfer.nw > transfer_values.txt
	'''
}
transtree.subscribe{
	file -> file.copyTo(resultDir.resolve(file.name))
}

/* Classical support to bootstrap trees */
process classicalSupport {
	input:
	file(divtree) from div1tree5
	file(boottrees) from boottrees2

	output:
	file("class_values.txt") into classvalues
	file("div_1_class.nw") into classtree

	shell:
	'''
	#!/bin/bash
	gotree compute support classical -i !{divtree} -b !{boottrees} -t !{task.cpus} > div_1_class.nw
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
	sort -k 1 !{quartetsanalyze} | cut -d ' ' -f 2,3,4,5,6 > quarts
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

#!/usr/bin/env nextflow

/* realignment */
params.alignment= "$baseDir/results/realignment.fasta.gz"
params.resultDir= 'results'

alignment = file(params.alignment)
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

/* 1, 32x32, 256*256 */
divisionChan = Channel.from([1, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32,32, 32, 32, 32, 32, 32, 32, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256])

/* Corresponding seeds */
seedChan = Channel.from([691607951, 8183053866, 4893742101, 4077898797, 343704747, 4873405240, 5273707124, 1339632896, 3130080935, 1374534207, 6246062756, 1496776910, 8878338027, 6323448214, 3237974410, 2987189437, 6222563888, 56088712, 612987361, 831390561, 17703152, 552425820, 200121068, 258751598, 560588769, 314545918, 723189768, 307371629, 224047533, 441116112, 711059685, 888806141, 505667153, 928157282, 459961156, 424359858, 626709707, 327352905, 540238733, 703575873, 18013801, 544901799, 48647327, 53413397, 170583724, 229153240, 577138662, 460553759, 273708681, 999136595, 298170956, 907835386, 867919107, 938410850, 133955531, 10212400, 911634335, 405190210, 512739722, 899914067, 212023562, 33359075, 50269235, 316149740, 450375210, 526021390, 360524958, 309915236, 714833447, 295450890, 596005100, 798408871, 435963321, 104610141, 874637946, 740731692, 88720025, 93091561, 926761478, 65555515, 406257263, 980661693, 928531359, 490519023, 705535084, 115891987, 478479455, 857079815, 905486958, 682627787, 958082332, 416827380, 828382589, 94572953, 211968712, 937210488, 908813635, 173582474, 19729740, 769444205, 181849579, 632313299, 314635474, 49173697, 438756197, 713287555, 373176894, 849259030, 757011844, 36726750, 701708388, 718558927, 775702230, 567807168, 89261920, 518494736, 225349766, 504449470, 508050957, 973313469, 46502735, 454725055, 313445844, 635341092, 851144904, 418922266, 84832117, 849207817, 721155316, 839146388, 692352770, 828668062, 476369830, 217204213, 131582369, 991035493, 427546531, 661224550, 137284653, 110420448, 416312353, 911929898, 782276324, 581411472, 621689283, 298959787, 595610657, 791193008, 980257120, 929983766, 880516178, 560474044, 645791528, 36809869, 732594423, 971820762, 104604049, 616415785, 108924739, 90704696, 723098283, 654550380, 119279360, 804807932, 802112239, 298379574, 263431666, 774964507, 550637942, 595317080, 714234017, 661345998, 673581034, 861167084, 785368360, 831649326, 173732130, 53329833, 979388868, 764940019, 49439511, 377808072, 503715200, 510603386, 474456613, 235323880, 346585948, 990384170, 696651861, 5995641, 650119708, 134059766, 423401210, 104029252, 627695789, 899971208, 403981963, 865497865, 346821802, 849040204, 573800984, 51739847, 811772856, 991018882, 154091715, 418071399, 413656407, 520011961, 157845656, 8081362, 158047669, 410155019, 667343138, 822150692, 705482051, 893843219, 38120548, 996268646, 130692697, 185347075, 299992613, 306259321, 968752663, 550076128, 832725075, 459353671, 959910967, 37430954, 269632355, 392823732, 166485873, 973140245, 862843339, 454858158, 649349021, 533096117, 187533061, 440966524, 526417351, 717674744, 178623730, 633458993, 177897361, 959593262, 669445421, 401437524, 994328730, 702782541, 879671908, 100997127, 816116046, 954951512, 569451763, 127820826, 998746022, 749155383, 124276265, 270958216, 120834664, 725969715, 253494423, 364335179, 629431349, 13438681, 717620698, 140197071, 361316781, 916196999, 228846804, 523546204, 993591906, 850175983, 409448801, 379405299, 495529579, 73331916, 866294544, 534004794, 439518353, 736995423, 437546102, 342737613, 246667894, 898636085, 904663166, 446344056, 896870295, 817490261, 967346316, 342505245, 698137517, 411763343, 72769971, 314859382, 342651881, 565800321, 993468606, 924524179, 721317759, 508548434, 65842731, 917537618, 17307824, 321455555, 749024539])

/* Rename the taxa of the alignment */
process renameAlignment {

	input : 
	file alignment
	
	output:
	set file("original.fa.gz"),file("name_map.txt.gz") into originalAlignment, originalAlignmentCopy

	shell:
	'''
	#!/usr/bin/env bash
	goalign trim name -i !{alignment} -o original.fa -m name_map.txt -n 20
	gzip original.fa
	gzip name_map.txt
	'''
}

originalAlignmentCopy.subscribe{
	file,file2 -> file.copyTo(alignmentDir.resolve(file.name)); file2.copyTo(alignmentDir.resolve(file2.name))
}

/**
  Divides the original alignment
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
dividedAlignment.into{dividedAlignmentToBoot; refAlignmentFastTree}

/**
	The Process that will reconstruct the Reference trees
*/
process runRefFastTree {
	tag "${refAlign} : div ${div} - seed ${seed}"

	input:
	set val(div), val(seed), file(refAlign) from refAlignmentFastTree

	output:
	set val(div), val(seed), file("refTree.nw.gz") into refTreeOutput
	
	shell:
	'''
	#!/usr/bin/env bash
	gunzip -c !{refAlign} > al.fa
	FastTree -nopr -nosupport -gtr -nt -gamma al.fa > refTree_tmp.nw
	raxmlHPC-PTHREADS -f e  -p ${RANDOM} -g refTree_tmp.nw -m GTRGAMMA -s al.fa -n TEST
	mv RAxML_result.TEST refTree.nw
	gzip refTree.nw
	rm -f al.fa al.phylip* RAxML*
	'''
}

refTreeOutput.subscribe{
	div, seed, fasttree -> fasttree.copyTo(treeDir.resolve("ref_"+div+"_"+seed+".nw.gz"))
}

/**
  We build the 1000 bootstraps alignment for a given reference alignment in a set of gz files
*/
process bootstrapAlignments {
	tag "${refAlign} : div ${div} - seed ${seed}"

	input :
	set val(div), val(seed), file(refAlign) from dividedAlignmentToBoot

	output:
	set val(div), val(seed), stdout, file("boot.tar.gz") into bootAlignment

	shell:
	'''
	#!/usr/bin/env bash
	goalign build seqboot -i !{refAlign} -n 1000 -o boot -S -s ${RANDOM} --tar --gz -t !{task.cpus}
	printf `goalign stats nseq -i !{refAlign}`
	'''
}

dividedBootAlignment = Channel.create()
groupedBootAlignment = Channel.create()

bootAlignment.choice( dividedBootAlignment, groupedBootAlignment ) { item -> item[2].toInteger() > 180 ? 0 : 1 }

process runGroupedBootFastTree {
	tag "${bootFile} : div ${div} - seed ${seed}"

	input:
	set val(div), val(seed), val(nbtaxa), file(bootFile) from groupedBootAlignment

	output:
	set val(div), val(seed), file("${bootFile.baseName}.nw.gz") into groupedBootTreeOutput
	
	shell:
	template 'bootstrapfastree.py'
}

/**
	We first divide the bootstrap trees then For each bootstrap alignment, we run FastTree
*/
process divideBootAlign{
	tag "${bootFile} : div ${div} - seed ${seed}"

	input:
	set val(div), val(seed), val(nbtaxa), file(bootFile) from dividedBootAlignment

	output:
	set val(div), val(seed), file("boot*.fa.gz") into dividedBootAlignToTree mode flatten
	
	shell:
	'''
	#!/usr/bin/env bash
	tar -xzf !{bootFile}
	gzip boot*.fa
	'''
}

process runBootFastTree {
	tag "${bootFile} : div ${div} - seed ${seed}"

	input:
	set val(div), val(seed), file(bootFile) from dividedBootAlignToTree

	output:
	set val(div), val(seed), file("${bootFile.baseName}.nw.gz") into bootTreeOutput
	
	shell:
	'''
	#!/usr/bin/env bash
	gunzip -c !{bootFile}                        \
         | FastTree -nopr -nosupport -gtr -nt -gamma \
         | gzip -c -  > !{bootFile.baseName}.nw.gz
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
and Concatenate all bootstrap trees into a single file
*/
process concatBootTreeFiles {
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

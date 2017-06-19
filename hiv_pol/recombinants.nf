#!/usr/bin/env nextflow

params.datadir  = "$baseDir/results/"
params.resultDir= 'results'
params.alignment= "$baseDir/results/realignment.fasta.gz"

/* Output directories */
resultDir    = file([params.resultDir, "recombinants"].join(File.separator))

resultDir.with {
    mkdirs()
}

sequences = Channel.fromPath(params.alignment).splitFasta( by: 5, file: true )

process printSequences{

	scratch true

	input:
	file(chunk) from sequences

	output:
	file("recomb.txt") into recombinants

	shell:
	'''
	!{task.jphmmpath}/src/jpHMM -s !{chunk} -v HIV -P !{task.jphmmpath}/priors -I !{task.jphmmpath}/input
	grep -v "^#" output/recombination.txt > recomb.txt
	rm -rf output
	'''
}

recombinants.collectFile(name: 'recombinants.txt').into{recombtocopy ; recombtoformat}

recombtocopy.subscribe{
	file -> file.copyTo(resultDir.resolve(file.name))
}


process formatrecombinants{

	input:
	file(recomb) from recombtoformat

	output:
	file("formated.txt") into recombformated
	
	shell:
	template 'formatrecomb.pl'
}


recombformated.subscribe{
	file -> file.copyTo(resultDir.resolve(file.name))
}

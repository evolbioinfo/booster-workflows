#!/usr/bin/env nextflow

params.alignment= "$baseDir/../data/vih/All.pol.aln.fasta.gz"  /* Full alignment */
params.hxb2     = "$baseDir/../data/vih/HXB2.fasta"            /* hxb2 sequence, useful to remove drug resistance positions */
params.taxafile = "$baseDir/../data/vih/pol_nonrecombinant.txt"/* File that contains taxa to keep from
                                                                 the alignment (non recombinant subtypes) */
params.resultDir= 'results'

alignment = file(params.alignment)
taxafile  = file(params.taxafile)
hxb2      = file(params.hxb2)
resultDir = file(params.resultDir)

resultDir.with {mkdirs()}

/* Selects good taxa and unalign the alignment */
/* And adds HXB2 sequence */
process selectTaxa {
	input :
	file alignment
	file taxafile
	file hxb2
	
	output :
	file("seq_selected.fa.gz") into selectedAlignment

	shell:
	'''
	#!/usr/bin/env bash
	sed 's/>.*/>HXB2/g' !{hxb2} > tmp.fa

        goalign subset -i !{alignment} -f !{taxafile} \
              | goalign unalign \
              > seqs.fa
	cat seqs.fa tmp.fa | gzip -c - > seq_selected.fa.gz
	rm tmp.fa seqs.fa
	'''
}

/* Realigns the sequences */
process realignSequences {
	input:
	file(sequences) from selectedAlignment
	
	output:
	file("alignment.fasta.gz") into realignment

	shell:
	'''
	gunzip -c !{sequences} > seqs.fa
	mafft --thread !{task.cpus} seqs.fa > alignment.fasta
	gzip alignment.fasta
	rm seqs.fa
	'''
}

/* Remove Drug Resistance Mutations */
process removeDRM {
	input:
	file(align) from realignment
	
	output:
	file("alignment_drm.fasta.gz") into drmrealignment

	scratch true

	shell:
	'''
	#!/usr/bin/env Rscript
	library(ape)
	library(big.phylo)
	system('gunzip -c !{align} > al.fa')
	ali=read.dna(file="al.fa", format="fa")
	out=seq.rm.drugresistance(ali, outfile="alignment_drm.RData")
	write.dna(file="alignment_drm.fasta",out$nodr.seq, format='fasta', colsep='', nbcol=-1)
	system('gzip alignment_drm.fasta')
	'''
}

process removeHXB2 {
	input:
	file(align) from drmrealignment

	output:
	file("realignment.fasta.gz") into outrealignment

	shell:
	'''
	goalign subset -r -i !{align} HXB2 | gzip -c - > realignment.fasta.gz
	'''
}

outrealignment.subscribe{
	file -> file.copyTo(resultDir.resolve(file.name))
}

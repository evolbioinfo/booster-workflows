# Realign without known recombinants and remove drug resistance mutations
nextflow run realign.nf --alignment "../data/vih/All.pol.aln.fasta.gz"  \
                        --hxb2 "../data/vih/HXB2.fasta"                 \
                        --taxafile "../data/vih/pol_nonrecombinant.txt" \
                        --resultDir "results"

# Infer trees (reference & bootstraps)
nextflow run trees.nf --alignment "results/realignment.fasta.gz" \
                      --resultDir "results"

# Compute support
nextflow run supports.nf --datadir "results"   \
                         --resultDir "results" 

# Detection of recombinants with jpHMM
nextflow run recombinants.nf --datadir "results"                        \
                             --resultDir "results"                      \
                             --alignment "results/realignment.fasta.gz"


# Plots
nextflow run plots.nf --tbetree "results/supports/1_691607951_tbe.nw"     \
	              --fbptree "results/supports/1_691607951_fbp.nw"     \
		      --namemap "results/aligns/name_map.txt.gz"          \
		      --resultdir "results/plots"                         \
		      --recomb "results/formated.txt"                     \
		      --boosterlog "results/supports/1_691607951_tbe.log" \
		      --allsupports "results/allsupports.txt"

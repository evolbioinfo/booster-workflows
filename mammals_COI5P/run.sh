# We infer trees (reference & bootstraps) using fasttree for
# ref and bootstrap alignments
nextflow run trees.nf -c nextflow_configs/trees.config --datadir ../data/mammals --resultDir results/fasttree

# We also infer ref and boot trees using RAxML (+rapid bootstrap)
nextflow run raxmltrees.nf -c nextflow_configs/raxmltrees.config --resultDir results/raxml

# We compute support for FastTree trees
nextflow run supports.nf -c nextflow_configs/supports.config --datadir results/fasttree --resultdir results/fasttree

# We compute support for RAxML rapid bootstrap
nextflow run supports.nf -c nextflow_configs/supports.config --datadir results/raxml --resultdir results/raxml

# We compare FastTree inferred tree with ncbi tree (quartets)
nextflow run compare_quartets.nf -c nextflow_configs/quartets.config                   \
                                 --div1tree results/fasttree/trees/ref_1_31144.nw.gz   \
                                 --truetree  ../data/ncbitax/ncbi.nw                   \
                                 --boottrees results/fasttree/trees/boot_1_31144.nw.gz \
                                 --mapfile   results/fasttree/aligns/name_map.txt.gz   \
                                 --resultdir results/fasttree/quartets

# We compare RAxML inferred (rapid bootstrap) tree with ncbi tree (quartets)
nextflow run compare_quartets.nf -c nextflow_configs/quartets.config                \
                                 --div1tree  results/raxml/trees/ref_1_31144.nw.gz  \
                                 --truetree  ../data/ncbitax/ncbi.nw                \
                                 --boottrees results/raxml/trees/boot_1_31144.nw.gz \
                                 --mapfile   results/raxml/aligns/name_map.txt.gz   \
                                 --resultdir results/raxml/quartets

# We generate figures
nextflow run plot.nf -c nextflow_configs/plots.config                                \
                     --tbetreefasttree "results/fasttree/support/1_31144_tbe.nw"     \
                     --fbptreefasttree "results/fasttree/support/1_31144_fbp.nw"     \
                     --namemapfasttree "results/fasttree/aligns/name_map.txt.gz"     \
                     --tbetreeraxml "results/raxml/support/1_31144_tbe.nw"           \
                     --fbptreeraxml "results/raxml/support/1_31144_fbp.nw"           \
                     --namemapraxml "results/raxml/aligns/name_map.txt.gz"           \
                     --resultdir "results/plots"                                     \
                     --ncbi "../data/ncbitax/ncbi_labels.nw"                         \
                     --allsupportfasttree "results/fasttree/allsupports.txt"         \
                     --allsupportraxml "results/raxml/allsupports.txt"               \
                     --conflictfasttree "results/fasttree/quartets/groupedinfos.txt" \
                     --conflictraxml "results/raxml/quartets/groupedinfos.txt"       \
		     --raxmltree "results/raxml/trees/ref_1_31144.nw.gz"             \
		     --raxmlboottrees "results/raxml/trees/boot_1_31144.nw.gz"       \
		     --fasttreetree "results/fasttree/trees/ref_1_31144.nw.gz"       \
		     --fasttreeboottrees "results/fasttree/trees/boot_1_31144.nw.gz"

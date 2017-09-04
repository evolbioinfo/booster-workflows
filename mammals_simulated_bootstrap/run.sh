########################################################
#               Simulations with INDELible             #
########################################################

nextflow run trees.nf -c nextflow_configs/trees.config                    \
	              --datadir '../mammals_COI5P/results/fasttree/trees' \
	              --resultdir 'results/indelible'                     \
	              --rateshuffle 0.0                                   \
                      --raterecombi 0.0                                   \
                      --lengthrecombi 0                                   \
                      --raterogue 0.00                                    \
                      --simulator "indelible"

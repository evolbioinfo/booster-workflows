########################################################
# Simulations with INDELible, without additional noise #
########################################################

nextflow run trees.nf -c nextflow_configs/trees.config \
	              --resultdir 'results/indelible'  \
	              --rateshuffle 0.0                \
                      --raterecombi 0.0                \
                      --lengthrecombi 0                \
                      --raterogue 0.00                 \
                      --simulator "indelible"

nextflow run supports.nf -c nextflow_configs/supports.config                                          \
                         --datadir 'results/indelible'                                                \
                         --resultdir 'results/indelible'                                              \
                         --truetree '../mammals_COI5P/results/fasttree/trees/ref_phyml_1_31144.nw.gz' \
                         --roguefile 'results/indelible/aligns/rogues.txt'

nextflow run compare_quartets.nf -c nextflow_configs/quartets.config                                         \
                                 --div1tree 'results/indelible/trees/ref_1_65512081.nw.gz'                   \
                                 --truetree '../mammals_COI5P/results/fasttree/trees/ref_phyml_1_31144.nw.gz'\
                                 --boottrees 'results/indelible/trees/boot_1_65512081.nw.gz'                 \
                                 --resultdir 'results/indelible/quartets'                                    \
                                 --roguefile 'results/indelible/aligns/rogues.txt'

###############################################################################################
# Simulations with INDELible with 50% additional noise (and 50% more noise on 5% of the taxa) #
###############################################################################################

nextflow run trees.nf  -c nextflow_configs/trees.config \
                      --resultdir 'results/indelible0.5'\
                      --rateshuffle 0.5                 \
                      --raterecombi 0.0                 \
                      --lengthrecombi 0                 \
                      --raterogue 0.05                  \
                      --simulator "indelible"

nextflow run supports.nf -c nextflow_configs/supports.config                                          \
                         --datadir 'results/indelible0.5'                                             \
                         --resultdir 'results/indelible0.5'                                           \
                         --truetree '../mammals_COI5P/results/fasttree/trees/ref_phyml_1_31144.nw.gz' \
                         --roguefile 'results/indelible0.5/aligns/rogues.txt'

nextflow run compare_quartets.nf -c nextflow_configs/quartets.config                                          \
                                 --div1tree 'results/indelible0.5/trees/ref_1_65512081.nw.gz'                 \
                                 --truetree '../mammals_COI5P/results/fasttree/trees/ref_phyml_1_31144.nw.gz' \
                                 --boottrees 'results/indelible0.5/trees/boot_1_65512081.nw.gz'               \
                                 --resultdir 'results/indelible0.5/quartets'                                  \
                                 --roguefile 'rogue'

###############################################################
#                      Plots Figures                          #
###############################################################

# Figure S12
nextflow run plot.nf -c nextfow_configs/plots.config                                   \
                      --tbetree "results/indelible/supports/1_65512081_tbe.nw"         \
                      --fbptree "results/indelible/supports/1_65512081_fbp.nw"         \
                      --noisytbetree "results/indelible0.5/supports/1_65512081_tbe.nw" \
                      --noisyfbptree "results/indelible0.5/supports/1_65512081_tbe.nw" \
                      --resultdir "results/plots"                                      \
                      --allsupports "results/indelible/allsupports.txt"                \
                      --conflicts "results/indelible/quartets/groupedinfos.txt"        \
                      --noisyallsupports "results/indelible/allsupports.txt"           \
                      --noisyconflicts "results/indelible/quartets/groupedinfos.txt"   \
                      --boosterlog "results/indelible0.5/supports/1_65512081_tbe.log"  \
                      --roguefile "results/indelible0.5/aligns/rogues.txt"

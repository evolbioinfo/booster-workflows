# BOOSTER - BOOtstrap Support by TransfER - Workflows

This git repository contains the workflows used in the study "Smooting Felsenstein's Phylogenetic Bootstrap in the Era of Big Data".

They are implemented in [Nextflow](https://www.nextflow.io/).

* Original data is located in the `data` folder.
* Following folders correspond to one analysis:
    * hiv_pol: analysis of 9147 sequences of hiv pol
    * mammals_COI5P: analysis of 1449 sequences of COI-5P protein in mammals;
    * mammals_simulated: analysis of simulated data
    * transfer_distance: analysis of transfer distance as a function of branch depth

Each folder contains a `run.sh` script to launch the analysis.

To run all the pipelines, the following dependencies are necessary:

* [goalign](https://github.com/fredericlemoine/goalign)
* [gotree](https://github.com/fredericlemoine/gotree)
* [booster](https://github.com/fredericlemoine/booster)

* [PhyML 20120412](http://www.atgc-montpellier.fr/phyml/download.php)
* [FastTree 2.1.8](http://www.microbesonline.org/fasttree/)
* [raxmlHPC-PTHREADS v8.2.8](http://sco.h-its.org/exelixis/software.html)
* python 2.x
* perl 5.x
* [mafft 7](http://mafft.cbrc.jp/alignment/software/)
* R 3.2.3 with the following packages:
    * [ape](https://cran.rstudio.com/web/packages/ape/index.html)
    * [big.phylo](https://github.com/olli0601/big.phylo)
    * [reshape2](https://cran.r-project.org/web/packages/reshape2/index.html)
    * [ggplot2](http://ggplot2.org/)
* [jpHMM](http://jphmm.gobics.de/)
* [tqdist](http://users-cs.au.dk/cstorm/software/tqdist/)
* [Seq-Gen](https://github.com/rambaut/Seq-Gen)
* [INDELible](http://abacus.gene.ucl.ac.uk/software/indelible/)
* [Nextflow](https://www.nextflow.io/)

#!/usr/bin/env bash

# Yultrees
echo "yule 16"
nextflow run transfer.nf --nbrand=1000 --nbref=100 --size=16   --divisions 1   --reftree="caterpillartree" --randtree="yuletree" -w work01 -resume
echo "yule 128"
nextflow run transfer.nf --nbrand=1000 --nbref=100 --size=128  --divisions 1   --reftree="caterpillartree" --randtree="yuletree" -w work02 -resume
echo "yule 1024"
nextflow run transfer.nf --nbrand=1000 --nbref=100 --size=1024 --divisions 10  --reftree="caterpillartree" --randtree="yuletree" -w work03 -resume
echo "yule 8192"
nextflow run transfer.nf --nbrand=1000 --nbref=100 --size=8192 --divisions 100 --reftree="caterpillartree" --randtree="yuletree" -w work04 -resume

# Caterpillars
echo "caterpillar 16"
nextflow run transfer.nf --nbrand=1000 --nbref=100 --size=16   --divisions 1   --reftree="caterpillartree" --randtree="caterpillartree" -w work05 -resume
echo "caterpillar 128"
nextflow run transfer.nf --nbrand=1000 --nbref=100 --size=128  --divisions 1   --reftree="caterpillartree" --randtree="caterpillartree" -w work06 -resume
echo "caterpillar 1024"
nextflow run transfer.nf --nbrand=1000 --nbref=100 --size=1024 --divisions 10  --reftree="caterpillartree" --randtree="caterpillartree" -w work07 -resume
echo "caterpillar 8192"
nextflow run transfer.nf --nbrand=1000 --nbref=100 --size=8192 --divisions 100 --reftree="caterpillartree" --randtree="caterpillartree" -w work08 -resume

# Balanced
echo "balanced 16"
nextflow run transfer.nf --nbrand=1000 --nbref=100 --size=16   --divisions 1   --reftree="caterpillartree" --randtree="balancedtree" -w work09 -resume
echo "balanced 128"
nextflow run transfer.nf --nbrand=1000 --nbref=100 --size=128  --divisions 1   --reftree="caterpillartree" --randtree="balancedtree" -w work10 -resume
echo "balanced 1024"
nextflow run transfer.nf --nbrand=1000 --nbref=100 --size=1024 --divisions 10  --reftree="caterpillartree" --randtree="balancedtree" -w work11 -resume
echo "balanced 8192"
nextflow run transfer.nf --nbrand=1000 --nbref=100 --size=8192 --divisions 100 --reftree="caterpillartree" --randtree="balancedtree" -w work16 -resume

# Uniform
echo "uniform 16"
nextflow run transfer.nf --nbrand=1000 --nbref=100 --size=16   --divisions 1   --reftree="caterpillartree" --randtree="uniformtree" -w work12 -resume
echo "uniform 128"
nextflow run transfer.nf --nbrand=1000 --nbref=100 --size=128  --divisions 1   --reftree="caterpillartree" --randtree="uniformtree" -w work13 -resume
echo "uniform 1024"
nextflow run transfer.nf --nbrand=1000 --nbref=100 --size=1024 --divisions 10  --reftree="caterpillartree" --randtree="uniformtree" -w work14 -resume
echo "uniform 8192"
nextflow run transfer.nf --nbrand=1000 --nbref=100 --size=8192 --divisions 100 --reftree="caterpillartree" --randtree="uniformtree" -w work15 -resume

# Plots
Rscript plot.Rscript 16
Rscript plot.Rscript 128
Rscript plot.Rscript 1024
Rscript plot.Rscript 8192

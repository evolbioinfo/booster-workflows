#!/usr/bin/env bash

gunzip -c !{tree} > tree.nh

for i in {1..1000}
do
    cat > control.txt <<EOF
[TYPE] AMINOACID 2
[SETTINGS]
EOF
    echo "  [randomseed] $((!{seed}+$i))" >> control.txt
    cat >> control.txt <<EOF
[MODEL]    modelname
  [submodel] !{genmodel}
  [indelmodel]  POW  1.5 5
  [rates] 0.0 0.441 4
  [indelrate]   0.02
  [statefreq] 0.046349 0.007288 0.019407 0.015423 0.000257 0.008818 0.006220 0.048443 0.013156 0.039273 0.070358 0.005750 0.030595 0.035916 0.031813 0.029183 0.034057 0.013934 0.013829 0.034391
EOF
    TREE=`cat tree.nh`
    echo "[TREE] treename $TREE" >> control.txt
    
    cat >> control.txt <<EOF
[PARTITIONS] partitionname
  [treename modelname 250]
[EVOLVE] partitionname 1 outputname
EOF
    
    indelible

    goalign reformat fasta -p -i outputname_TRUE.phy | gzip -c -  > bootsimu_${i}.fa.gz
    rm -f outputname* control.txt trees.txt
done
tar -cf bootsimu.tar bootsimu_* > /dev/null
rm -f bootsimu_*

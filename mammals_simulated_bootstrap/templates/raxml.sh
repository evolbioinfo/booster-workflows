#!/usr/bin/env bash
goalign reformat phylip -i !{inalign} -o al.phy
raxmlHPC-PTHREADS -f d -p ${RANDOM} -m PROTGAMMAWAG -c 6 -s al.phy -n TEST -T !{task.cpus}
mv RAxML_result.TEST !{outtree}
gzip !{outtree}
rm -f al.phy_* al.phy

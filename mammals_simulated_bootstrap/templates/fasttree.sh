#!/bin/env bash
goalign reformat fasta -i !{inalign} -o al.fa
FastTree -nopr -nosupport -wag -gamma al.fa > !{outtree}
gzip !{outtree}
rm -f al.fa

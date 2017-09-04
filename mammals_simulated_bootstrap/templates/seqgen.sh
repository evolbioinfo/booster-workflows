#!/usr/bin/env bash

gunzip -c !{tree} > tree.nh
seq-gen -of -a0.441 -g4  \
    -m!{genmodel}        \
    -l!{seqlen}          \
    -z ${RANDOM} tree.nh \
    | gzip -c -          \
    > original.fa.gz
rm tree.nh

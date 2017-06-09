#!/usr/bin/env bash

fbpurl=$(gotree rename -m !{namemap} -i !{fbptree} -r | gotree clear pvalues | gotree upload itol !{itolannot})
tbeurl=$(gotree rename -m !{namemap} -i !{tbetree} -r | gotree clear pvalues | gotree upload itol !{itolannot})

gotree dlimage itol -f svg -i $(basename $tbeurl) -c !{configtbe} -o !{outprefix}_tree_tbe.svg
gotree dlimage itol -f svg -i $(basename $fbpurl) -c !{configfbp} -o !{outprefix}_tree_fbp.svg


#!/usr/bin/env bash

# Just for tree pruning afterwards (avoid doing each time the rename...)
gotree rename -r -i !{tbetree} -m !{namemap} > topo

## We create a tree with one internal branch per subtype and keep the same
## Set of taxa than the trees we will look at
#allsubtypetree.pl !{recomb} | gotree prune -c topo > allsubtypetree.nw
#
## We compare branches of the tree with this subtype tree
#gotree rename -i !{tbetree} -m !{namemap} -r                    \
#    | gotree compare edges -m --moved-taxa -c allsubtypetree.nw \
#    | awk -F "\t" '{if($7>1){print $0}}'                        \
#    | grep "_subtype" > topo_compare_tbe
#
#gotree rename -i !{fbptree} -m !{namemap} -r \
#    | gotree compare edges -m --moved-taxa -c allsubtypetree.nw \
#    | awk -F "\t" '{if($7>1){print $0}}' \
#    | grep "_subtype" > topo_compare_fbp

# For each subtype, we take the branch of the inferred tree having the least transfer distance
# We then look at its support
#tail -n+2 topo_compare_tbe | sort -u -k12 -k10,10n -t$'\t' | awk -F"\t" '!_[$12]++' > !{outprefix}_support_subtypes_tbe.txt
#tail -n+2 topo_compare_fbp | sort -u -k12 -k10,10n -t$'\t' | awk -F"\t" '!_[$12]++' > !{outprefix}_support_subtypes_fbp.txt

# We look at only branches that are supported at more than 70% and take the min distance to each subtype
# We first generate one branch trees
subtypetree.pl !{recomb} A1 A2 | gotree prune -c topo > A.nh
subtypetree.pl !{recomb} B | gotree prune -c topo > B.nh
subtypetree.pl !{recomb} C | gotree prune -c topo > C.nh
subtypetree.pl !{recomb} D | gotree prune -c topo > D.nh
subtypetree.pl !{recomb} F1 F2 | gotree prune -c topo > F.nh
subtypetree.pl !{recomb} K | gotree prune -c topo > K.nh
subtypetree.pl !{recomb} G | gotree prune -c topo > G.nh
subtypetree.pl !{recomb} H | gotree prune -c topo > H.nh
subtypetree.pl !{recomb} J | gotree prune -c topo > J.nh

# And compare the two trees to each of the nine subtypes
# to search for the closest supported branches (>0.7) or any branches to each subtype
subtypes=( A B C D F K G H J )
for s in "${subtypes[@]}"
do

    gotree compare distances -i ${s}.nh                       \
	   -c <(gotree rename -i !{tbetree} -m !{namemap} -r) \
	   > t1

    gotree compare distances -i ${s}.nh                       \
	   -c <(gotree rename -i !{fbptree} -m !{namemap} -r) \
	   > t2
		 
    cat t1 | tail -n+2 | sort -k 4 -n \
	| head -n 1  \
	       > !{outprefix}_support_subtype_${s}_tbe.txt
    
    cat t1 | tail -n+2 | sort -k 4 -n \
	| awk -F "\t" '{if($6>0.7){print $0}}' \
	| head -n 1  \
	      > !{outprefix}_high_support_subtype_${s}_tbe.txt

    cat t2 | tail -n+2 | sort -k 4 -n              \
	| head -n 1                                \
	       > !{outprefix}_support_subtype_${s}_fbp.txt
    
    cat t2 | tail -n+2 | sort -k 4 -n                   \
	| awk -F "\t" '{if($6>0.7){print $0}}'          \
	| head -n 1                                     \
	> !{outprefix}_high_support_subtype_${s}_fbp.txt
done

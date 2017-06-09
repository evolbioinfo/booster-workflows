#!/usr/bin/env bash

# We create annotation file for simian clade using data from ncbi taxonomy
gotree subtree -i !{ncbi} -n "^Anthropoidea.Simiiformes" | gotree stats tips | tail -n+2 | cut -f 4 > simians.txt
echo "simians:$(cat simians.txt | tr '\n' ',')" > map.txt

# We take the approximate simian clade from inferred tree (1 or 2 errors)
gotree rename -i !{tree} -m !{namemap} -r \
    | gotree reroot outgroup Ornithorhynchus_anatinus Zaglossus_bruijni Tachyglossus_aculeatus \
    | gotree annotate -m map.txt \
    | gotree subtree -n simians \
    | gotree stats tips \
    | tail -n+2 | cut -f 4 \
		      > simians_inferred.txt

// Real simian bipartition
gotree compute bipartitiontree -i !{ncbi} -f simians.txt > simians.nw
// Inferred simian bipartition
gotree rename -m !{namemap} -r -i !{tree} | gotree compute bipartitiontree -f simians_inferred.txt > simians_inferred.nw

// Function to compare the edges
function extracttaxa(){
    echo $1 | gotree rename -m !{namemap} -r \
	| gotree prune -c $2 \
	| gotree compare edges -m --moved-taxa -c $2 \
	| awk '{if($5=="false" && $11==""){print $0}}' \
	| sort -u -k10,10n -t$'\t' | head -n 1
}
export -f extracttaxa

// Exact simian clade
zcat !{boottrees} | parallel extracttaxa {} simians.nw > simiiformes_boot.txt
CANBMOVE=$(grep 'Canis_adustus' simiiformes_boot.txt | wc -l)
MRNBMOVE=$(grep 'Maxomys_rajah' simiiformes_boot.txt | wc -l)
# Number of taxa found in Simiens in at least 1 bootstrap tree
NBTAXMOVEONE=$(awk -F "\t" '{print $11}' simiiformes_boot.txt | sed 's/,/\n/g' | sed 's/[+-]//g' | sort -u  | wc)
# Number of taxa found in Simiens in at least 1% of bootstrap trees
NBTAXMOVEONEPERCENT=$(awk -F "\t" '{print $11}' simiiformes_boot.txt | sed 's/,/\n/g' | sed 's/[+-]//g' | sort | uniq -c | sort -n | awk '{if($1 > 10){print $0}}' | wc)

echo "Exact Simian clade:"
echo "    - Canis adustus moves=$CANBMOVE"
echo "    - Maxomys Rajah moves=$MRNBMOVE"
echo "    - Nb taxa moving in > 1  bootstrap tree=$NBTAXMOVEONE"
echo "    - Nb taxa moving in > 1% bootstrap tree=$NBTAXMOVEONEPERCENT"

// Inferred simian clade
zcat !{boottrees} | parallel extracttaxa {} simians_inferred.nw > simiiformes_boot_inferred.txt
CANBMOVE=$(grep 'Canis_adustus' simiiformes_boot_inferred.txt | wc -l)
MRNBMOVE=$(grep 'Maxomys_rajah' simiiformes_boot_inferred.txt | wc -l)
# Number of taxa found in Simiens in at least 1 bootstrap tree
NBTAXMOVEONE=$(awk -F "\t" '{print $11}' simiiformes_boot_inferred.txt | sed 's/,/\n/g' | sed 's/[+-]//g' | sort -u  | wc)
# Number of taxa found in Simiens in at least 1% of bootstrap trees
NBTAXMOVEONEPERCENT=$(awk -F "\t" '{print $11}' simiiformes_boot_inferred.txt | sed 's/,/\n/g' | sed 's/[+-]//g' | sort | uniq -c | sort -n | awk '{if($1 > 10){print $0}}' | wc)

echo "Inferred Simian clade:"
echo "    - Canis adustus moves=$CANBMOVE"
echo "    - Maxomys Rajah moves=$MRNBMOVE"
echo "    - Nb taxa moving in > 1  bootstrap tree=$NBTAXMOVEONE"
echo "    - Nb taxa moving in > 1% bootstrap tree=$NBTAXMOVEONEPERCENT"

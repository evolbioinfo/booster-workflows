#!/usr/bin/env bash

# Just for tree pruning afterwards (avoid doing each time the rename...)
gotree rename -r -i !{tbetree} -m !{namemap} > topo

# We look at only branches that are supported at more than 70% and take the min distance to each subtype
# We first generate one branch trees per interesting subtype / geographical area
geographicaltree.pl !{recomb} BI,DJ,ET,KE,SD,SO,TZ,UG C | gotree prune -c topo > C_EastAfrica
geographicaltree.pl !{recomb} MM,IN,NP C | gotree prune -c topo > C_India
geographicaltree.pl !{recomb} BR,UY,AR C | gotree prune -c topo > C_SouthAmerica
#geographicaltree.pl !{recomb} BY,GE,KZ,RU,UA,UZ A1 > A1_Russia

# Remove SouthAmerican sequences from east africa tree
grep -E "(UY|AR|BR)" formated.txt  | grep " C$" | cut -f 1 -d ' ' | sed 's/_/./g' > SA_sequences
gotree prune -i C_EastAfrica -f SA_sequences > C_EastAfrica_clean
gotree rename -i !{fbptree} -m !{namemap} -r | gotree prune -f SA_sequences > fbp_clean
gotree rename -i !{tbetree} -m !{namemap} -r | gotree prune -f SA_sequences > tbe_clean

# And compare the two trees to each of these clades
# to search for the closest supported branches (>0.7) or any branches to each subtype
subtypes=( C_EastAfrica C_India C_SouthAmerica )
for s in "${subtypes[@]}"
do
    gotree compare distances -i ${s}                       \
	   -c <(gotree rename -i !{tbetree} -m !{namemap} -r) \
	   > t1

    gotree compare distances -i ${s}                       \
	   -c <(gotree rename -i !{fbptree} -m !{namemap} -r) \
	   > t2
		 
    cat t1 | tail -n+2 | sort -k 4 -n \
	| head -n 1  \
	       > !{outprefix}_support_subtype_geo_${s}_tbe.txt
    
    cat t1 | tail -n+2 | sort -k 4 -n \
	| awk -F "\t" '{if($6>0.7){print $0}}' \
	| head -n 1  \
	      > !{outprefix}_high_support_subtype_geo_${s}_tbe.txt

    cat t2 | tail -n+2 | sort -k 4 -n              \
	| head -n 1                                \
	       > !{outprefix}_support_subtype_geo_${s}_fbp.txt
    
    cat t2 | tail -n+2 | sort -k 4 -n                   \
	| awk -F "\t" '{if($6>0.7){print $0}}'          \
	| head -n 1                                     \
	> !{outprefix}_high_support_subtype_geo_${s}_fbp.txt
done

# Same with clean trees
s=C_EastAfrica_clean
gotree compare distances -i ${s} -c tbe_clean > t1
gotree compare distances -i ${s} -c fbp_clean > t2
		 
cat t1 | tail -n+2 | sort -k 4 -n \
       | head -n 1                \
       > !{outprefix}_support_subtype_geo_${s}_tbe.txt
    
cat t1 | tail -n+2 | sort -k 4 -n \
       | awk -F "\t" '{if($6>0.7){print $0}}' \
       | head -n 1  \
       > !{outprefix}_high_support_subtype_geo_${s}_tbe.txt

cat t2 | tail -n+2 | sort -k 4 -n \
       | head -n 1                \
       > !{outprefix}_support_subtype_geo_${s}_fbp.txt

cat t2 | tail -n+2 | sort -k 4 -n             \
       | awk -F "\t" '{if($6>0.7){print $0}}' \
       | head -n 1                            \
       > !{outprefix}_high_support_subtype_geo_${s}_fbp.txt

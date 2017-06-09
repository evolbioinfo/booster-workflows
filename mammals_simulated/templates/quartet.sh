#!/usr/bin/env bash
CE=`quartet_dist -v !{truetree} !{truetree} | cut -f 7`
Q=`quartet_dist -v !{staredge} !{truetree}`
# Number of quartet resolved in both trees and that agree 
A=`echo $Q | cut -f 5 -d ' '`
# Quartet Distance  B+C+D
BCD=`echo $Q | cut -f 3 -d ' '`
# Number of quartets unresolved in both trees
E=`echo $Q | cut -f 7 -d ' '`
# Number of tips
TIPS=`echo $Q | cut -f 1 -d ' '`
# Topo depth
DEPTH=`gotree stats edges -i !{staredge} | awk '{if($5=="false"){print $7}}'`
# Number of resolved quartets in div1
ABC=$((ABC=DEPTH*(DEPTH-1)/2*(TIPS-DEPTH)*(TIPS-DEPTH-1)/2))
C=$((C=CE-E))
B=$((B=ABC-A-C))
ARAND=0
BRAND=0
for j in {1..10}
do
    gotree shuffletips -i !{staredge} -o tmp -s $RANDOM
    QTEMP=`quartet_dist -v tmp !{truetree}`
    ETEMP=`echo $QTEMP | cut -d ' ' -f 7`
    ATEMP=`echo $QTEMP | cut -d ' ' -f 5`
    ARAND=$((ARAND=ARAND+ATEMP))
    CTEMP=$((CTEMP=CE-ETEMP))
    BTEMP=$((BTEMP=ABC-ATEMP-CTEMP))
    BRAND=$((BRAND=BRAND+BTEMP))
done
echo !{staredge} $ABC $A $ARAND $B $BRAND > quartet

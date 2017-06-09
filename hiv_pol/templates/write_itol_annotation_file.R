#!/usr/bin/env Rscript

library(dplyr)
## annots.txt
a=read.table("!{recomb}")
## Out file
out="annot_itol.txt"

types=c("A1","A2","B","01_AE","C","D","F1","F2","G","H","J","K")
c=c("#ff0000","#950000","#0000ff","#00ff00","#996633","#00ffff","#800080","#cf00cf","#008000","#ff8000","#ff00ff","#ffff00")
a$col="#000000"

i=1
for(t in types){
    a$col[a$V2==t]=c[i]
    i=i+1
}

a$V1=gsub("_",".",a$V1)
a$t="branch"
a$style="normal"
a$width="5"
write("TREE_COLORS\nSEPARATOR SPACE\nDATASET_LABEL Subtypes\nCOLOR #ff0000\nDATA",file=out)
write.table(file=out,a[,c("V1","t","col","style","width")],col.names=F,row.names=F,quote = F,sep=" ",append=T)

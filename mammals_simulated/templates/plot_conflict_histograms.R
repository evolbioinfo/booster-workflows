#!/usr/bin/env Rscript
library(ggplot2)
library(reshape2)

cbPalette <- c("#e79f00", "#9ad0f3", "#0072B2", "#D55E00", "#CC79A7", "#F0E442", "#009E73")

conflict=read.table("!{conflict}")
colnames(conflict)=c("id","length","Depth","ncbi","totalquartets","agree","agreerand","disagree","disagreerand","true","TBE","FBP")

conflict$normquartet=((conflict$disagree)/(conflict$totalquartets))
conflict$c=cut(conflict$normquartet,breaks=c(0,0.2,1),include.lowest = T)
dat=melt(conflict, id.vars = c("c","normquartet","Depth"),measure.vars = c("FBP","TBE"))
dat$supcutoff=dat$value>!{cutoff}

dat1 = dcast(dat, c ~ variable*supcutoff, fun.aggregate = length,fill = 0, drop = FALSE)
dat1$total=dat1$TBE_TRUE+dat1$TBE_FALSE
dat2=dat1
dat1[,2:5]=dat1[,2:5]*100/dat1$total
dat.melt = melt(dat1, id.vars = "c", measure.vars = c("FBP_FALSE","TBE_FALSE","FBP_TRUE","TBE_TRUE"))
dat.melt$supcutoff=(dat.melt$variable=="FBP_TRUE" | dat.melt$variable=="TBE_TRUE")
dat.melt$type=substr(dat.melt$variable,1,3)

dat2.melt = melt(dat2, id.vars = "c", measure.vars = c("FBP_FALSE","TBE_FALSE","FBP_TRUE","TBE_TRUE"))

datmerged=merge(dat.melt,dat2.melt,by.x=c("c","variable"),by.y=c("c","variable"))


svg("!{outprefix}_!{cutoff}_conflict.svg")
ggplot(datmerged[datmerged$supcutoff,], aes(x = type, y = value.x, fill = type, label=value.y)) +
    geom_bar(position = "stack", stat = "identity") +
    geom_text()+
    facet_wrap( ~ c)+
    theme_bw()+
    scale_fill_manual(values=cbPalette)+
    scale_colour_manual(values=cbPalette)
dev.off()

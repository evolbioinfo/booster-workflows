#!/usr/bin/env Rscript
library(ggplot2)
library(reshape2)

cbPalette <- c("#e79f00", "#9ad0f3", "#0072B2", "#D55E00", "#CC79A7", "#F0E442", "#009E73")

supports=read.table("!{allsupport}")
colnames(supports)=c("div", "seed","collapse", "type","support","depth")

supports=supports[supports$div == 1 & supports$collapse=="not_collapsed",]
supports$supcutoff=supports$support>!{cutoff}

supports$d=cut(supports$depth,breaks=c(0, 4,16,max(supports$depth)))
dat = dcast(supports, d ~ type*supcutoff, fun.aggregate = length)

dat$total=dat$TBE_TRUE+dat$TBE_FALSE
dat2=dat
dat[,2:5]=dat[,2:5]*100/dat$total

dat.melt = melt(dat, id.vars = "d", measure.vars = c("FBP_FALSE","TBE_FALSE","FBP_TRUE","TBE_TRUE"))
dat.melt$supcutoff=(dat.melt$variable=="FBP_TRUE" | dat.melt$variable=="TBE_TRUE")
dat.melt$type=substr(dat.melt$variable,1,3)

dat2.melt = melt(dat2, id.vars = "d", measure.vars = c("FBP_FALSE","TBE_FALSE","FBP_TRUE","TBE_TRUE"))

datmerged=merge(dat.melt,dat2.melt,by.x=c("d","variable"),by.y=c("d","variable"))

svg("!{outprefix}_!{cutoff}_depth.svg")
ggplot(datmerged[datmerged$supcutoff,], aes(x = type, y = value.x, fill = type,label=round(value.y))) +
    geom_bar(position = "stack", stat = "identity") +
    geom_text()+
    facet_wrap( ~ d)+
    theme_bw()+
    scale_fill_manual(values=cbPalette)+
    scale_colour_manual(values=cbPalette)+
    xlab("Support type") +
    ylab("% branches")
dev.off()

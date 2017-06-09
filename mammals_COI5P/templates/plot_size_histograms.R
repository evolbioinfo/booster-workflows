#!/usr/bin/env Rscript
library(ggplot2)
library(reshape2)

cbPalette <- c("#e79f00", "#9ad0f3", "#0072B2", "#D55E00", "#CC79A7", "#F0E442", "#009E73")

supports=read.table("!{allsupport}")
colnames(supports)=c("div", "seed","collapse", "type","support","depth")

supports=supports[supports$collapse=="not_collapsed",]
supports$supcutoff=supports$support>!{cutoff}

dat = dcast(supports[supports$div%in%c(1,8,64),], div*seed ~ type*supcutoff, fun.aggregate = length)
dat$total=dat$TBE_TRUE+dat$TBE_FALSE
dat2=dat
dat2=aggregate(. ~ div, dat2[,c(1,3,4,5,6)], mean)

dat[,3:6]=dat[,3:6]*100/dat$total
dat$TBE_TRUE/dat$FBP_TRUE
dat=aggregate(. ~ div, dat[,c(1,3,4,5,6)], mean)

dat.melt = melt(dat, id.vars = "div", measure.vars = c("FBP_FALSE","TBE_FALSE","FBP_TRUE","TBE_TRUE"))
dat.melt$supcutoff=(dat.melt$variable=="FBP_TRUE" | dat.melt$variable=="TBE_TRUE")
dat.melt$type=substr(dat.melt$variable,1,3)

dat2.melt = melt(dat2, id.vars = "div", measure.vars = c("FBP_FALSE","TBE_FALSE","FBP_TRUE","TBE_TRUE"))

datmerged=merge(dat.melt,dat2.melt,by.x=c("div","variable"),by.y=c("div","variable"))

svg("!{outprefix}_!{cutoff}_size.svg")
ggplot(datmerged[datmerged$supcutoff,], aes(x = type, y = value.x, fill = type, label=floor(value.y))) +
    geom_bar(position = "stack", stat = "identity") +
    geom_text()+
    facet_wrap( ~ floor(1449/div))+
    theme_bw()+
    scale_fill_manual(values=cbPalette)+
    scale_colour_manual(values=cbPalette)+
    xlab("Support type") +
    ylab("% branches")
dev.off()

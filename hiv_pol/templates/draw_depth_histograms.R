#!/usr/bin/env Rscript

supports=read.table("!{allsupports}",header=F)
colnames(supports)=c("div", "seed","collapse", "type","support","depth")

tbe=supports[supports$type=="TBE" & supports$support>0.7 & supports$depth>1 & supports$div==1,]
fbp=supports[supports$type=="FBP" & supports$support>0.7 & supports$depth>1 & supports$div==1,]

svg("Figure1_histogram_TBE.svg")
h=hist(tbe$depth,breaks=c(0,2,4,8,32,max(tbe$depth)),plot=F)
bp=barplot(h$counts,names.arg=c("[2]","]2,4]","]4,8]","]8,32]",paste0("]32,",max(tbe$depth),"]")))
text(x=bp, y=h$counts, labels=round(h$counts,0), pos=3, xpd=NA)
dev.off()

svg("Figure1_histogram_FBP.svg")
h=hist(fbp$depth,breaks=c(0,2,4,8,32,max(fbp$depth)),plot=F)
bp=barplot(h$counts,names.arg=c("[2]","]2,4]","]4,8]","]8,32]",paste0("]32,",max(fbp$depth),"]")))
text(x=bp, y=h$counts, labels=round(h$counts,0), pos=3, xpd=NA)
dev.off()

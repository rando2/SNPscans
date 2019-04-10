library("readr")
library("ggplot2")
stringsAsFactors = FALSE

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

wildfile="karlsonfst.quantfilt.window500000.TA.tsv"
deepfile="karlsonfst.quantfilt.window500000.tsv"
deepTA <- read_delim(deepfile, delim="\t")
wildTA <- read_delim(paste("~/wild/2018/vv2align/fst/",wildfile, sep=""), delim="\t")

ovlap <- merge(deepTA, wildTA, by=c("CHROM","STARTPOS","ENDPOS"))

jpeg("hist-deepFST.jpg", width = 4, height = 4, units = 'in', res = 300)
p <- ggplot(deepTA, aes(x=KARLSSON_FST)) + geom_histogram(position="identity", alpha=0.5) + xlab("FST") + ggtitle("FST Distribution, Deep Tame & Aggr")+ theme(plot.title = element_text(hjust = 0.5)) + geom_vline(xintercept=quantile(x=deepTA$KARLSSON_FST, probs=0.995))
print(p)
dev.off()

jpeg("hist-poolFST.jpg", width = 4, height = 4, units = 'in', res = 300)
p <- ggplot(wildTA, aes(x=KARLSSON_FST)) + geom_histogram(position="identity", alpha=0.5) + xlab("FST") + ggtitle("FST Distribution, Pooled Tame & Aggr")+ theme(plot.title = element_text(hjust = 0.5)) + geom_vline(xintercept=quantile(x=wildTA$KARLSSON_FST, probs=0.995))
print(p)
dev.off()

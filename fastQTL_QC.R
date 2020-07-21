d = read.table("Naive/permutations.all.chunks.txt", hea=F, stringsAsFactors=F)
colnames(d) = c("pid", "nvar", "shape1", "shape2", "dummy", "sid", "dist", "npval", "slope", "ppval", "bpval")
plot(d$ppval, d$bpval, xlab="Direct method", ylab="Beta approximation", main="Check plot")
abline(0, 1, col="red")
d$bonferroni = p.adjust(d$bpval, method="bonferroni")
sign <- d[which(d$bonferroni <= 0.05), c(1,6,12)]
#1027 associations p-value < 0.05

write.table(x = sign, file = "Bonferroni_sign_naive.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

library(illuminaHumanv4.db)

x <- illuminaHumanv4ENTREZREANNOTATED
# Get the manufacturer identifiers that are mapped to a GO ID
mapped_genes <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_genes])

unlist(xx[sign$pid])

library(clusterProfiler)
#GO
c5 <- read.gmt("c5.bp.v7.1.entrez.gmt")
enrich_c5 <- enricher(unlist(xx[sign$pid]),
                          TERM2GENE=c5, pvalueCutoff = 0.05, pAdjustMethod = "BH", universe = unlist(xx[d$pid]))
enrich_c5_df<- as.data.frame(enrich_c5)
#enrich_c5_subset<- subset(enrich_c5, Count>5)

library(grDevices)
win.graph(width=40, height=50)

dotplot(enrich_c5, font.size=6)

##########################################################################################
d = read.table("LPS/permutations.all.chunksLPS2.txt", hea=F, stringsAsFactors=F)
colnames(d) = c("pid", "nvar", "shape1", "shape2", "dummy", "sid", "dist", "npval", "slope", "ppval", "bpval")
plot(d$ppval, d$bpval, xlab="Direct method", ylab="Beta approximation", main="Check plot")
abline(0, 1, col="red")
d$bonferroni = p.adjust(d$bpval, method="bonferroni")
sign <- d[which(d$bonferroni <= 0.05), c(1,6,12)]
#560 associations p-value < 0.05

write.table(x = sign, file = "Bonferroni_sign_LPS2.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)


library(illuminaHumanv4.db)

x <- illuminaHumanv4ENTREZREANNOTATED
# Get the manufacturer identifiers that are mapped to a GO ID
mapped_genes <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_genes])

library(clusterProfiler)
c5 <- read.gmt("c5.bp.v7.1.entrez.gmt")
enrich_c5 <- enricher(unlist(xx[sign$pid]),
                      TERM2GENE=c5, pvalueCutoff = 0.05, pAdjustMethod = "BH", universe = unlist(xx[d$pid]))
enrich_c5_df<- as.data.frame(enrich_c5)
enrich_c5_subset<- subset(enrich_c5, Count>5)

library(grDevices)
win.graph(width=40, height=50)

dotplot(enrich_c5, font.size=6)
#########################################################################################
d = read.table("LPS/permutations.all.chunksLPS24.txt", hea=F, stringsAsFactors=F)
colnames(d) = c("pid", "nvar", "shape1", "shape2", "dummy", "sid", "dist", "npval", "slope", "ppval", "bpval")
plot(d$ppval, d$bpval, xlab="Direct method", ylab="Beta approximation", main="Check plot")
abline(0, 1, col="red")
d$bonferroni = p.adjust(d$bpval, method="bonferroni")
sign <- d[which(d$bonferroni <= 0.05), c(1,6,12)]
#720 associations p-value < 0.05

write.table(x = sign, file = "Bonferroni_sign_LPS24.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

library(illuminaHumanv4.db)

x <- illuminaHumanv4ENTREZREANNOTATED
# Get the manufacturer identifiers that are mapped to a GO ID
mapped_genes <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_genes])


library(clusterProfiler)
c5 <- read.gmt("c5.bp.v7.1.entrez.gmt")
enrich_c5 <- enricher(unlist(xx[sign$pid]),
                      TERM2GENE=c5, pvalueCutoff = 0.05, pAdjustMethod = "BH", universe = unlist(xx[d$pid]))
enrich_c5_df<- as.data.frame(enrich_c5)
#enrich_c5_subset<- subset(enrich_c5, Count>5)

library(grDevices)
win.graph(width=40, height=50)

dotplot(enrich_c5, font.size=6)

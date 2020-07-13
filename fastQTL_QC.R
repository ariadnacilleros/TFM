d = read.table("Naive/permutations.all.chunks.txt", hea=F, stringsAsFactors=F)
colnames(d) = c("pid", "nvar", "shape1", "shape2", "dummy", "sid", "dist", "npval", "slope", "ppval", "bpval")
plot(d$ppval, d$bpval, xlab="Direct method", ylab="Beta approximation", main="Check plot")
abline(0, 1, col="red")
d$bonferroni = p.adjust(d$bpval, method="bonferroni")
sign <- d[which(d$bonferroni <= 0.05), c(1,6,12)]
#1027 associations p-value < 0.05

write.table(x = sign, file = "Bonferroni_naive_raw.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


d = read.table("LPS/permutations.all.chunksLPS2.txt", hea=F, stringsAsFactors=F)
colnames(d) = c("pid", "nvar", "shape1", "shape2", "dummy", "sid", "dist", "npval", "slope", "ppval", "bpval")
plot(d$ppval, d$bpval, xlab="Direct method", ylab="Beta approximation", main="Check plot")
abline(0, 1, col="red")
d$bonferroni = p.adjust(d$bpval, method="bonferroni")
sign <- d[which(d$bonferroni <= 0.05), c(1,6,12)]
#560 associations p-value < 0.05


d = read.table("LPS/permutations.all.chunksLPS24.txt", hea=F, stringsAsFactors=F)
colnames(d) = c("pid", "nvar", "shape1", "shape2", "dummy", "sid", "dist", "npval", "slope", "ppval", "bpval")
plot(d$ppval, d$bpval, xlab="Direct method", ylab="Beta approximation", main="Check plot")
abline(0, 1, col="red")
d$bonferroni = p.adjust(d$bpval, method="bonferroni")
sign <- d[which(d$bonferroni <= 0.05), c(1,6,12)]
#720 associations p-value < 0.05
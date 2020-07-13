#########################
#
#STEP:PREPARE .EPI FOR SMR
#
##########################

library(illuminaHumanv4.db)
library(stringr)
NA_LPS2 <- read.table("C:/Users/ariad/OneDrive/Desktop/Illumina_annotation/LPS/NA_LPS2.txt", quote="\"", comment.char="")
phenotype <- read.delim("C:/Users/ariad/OneDrive/Desktop/Illumina_annotation/LPS/phenotype_LPS2.txt")

probes <- setdiff(phenotype$PROBE_ID, NA_LPS2$V1)

#get chr, start curated
df <- phenotype[phenotype$PROBE_ID %in% probes, c(1,2,4)]

#get strand
x <- illuminaHumanv4GENOMICLOCATION
# Get the probe identifiers that are mapped to any cytoband
mapped_probes <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_probes])

strand_vector <- c()
for (probe in df$PROBE_ID){
  strand <- str_split(xx[probe], ":", n = , simplify = FALSE)[[1]][4]
  strand_vector <- append(strand_vector, strand)
}

df <- cbind(df, strand_vector)

write.table(x = df, file = "probes_info_LPS2.txt", sep = "\t", quote = FALSE, row.names = FALSE)


#ILMN_3288752
x <- illuminaHumanv4SYMBOL
# Get the probe identifiers that are mapped to a gene symbol
mapped_probes <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_probes])

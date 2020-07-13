naive <- read.table("SMR_Immuno_Imp/smr_naive_immuno.smr", header = TRUE, sep ="\t")

naive[order(naive$p_SMR), ] #order by p-value 

library(illuminaHumanv4.db)

x <- illuminaHumanv4SYMBOLREANNOTATED
# Get the probe identifiers that are mapped to a gene symbol
mapped_probes <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_probes])

gene_symbol <- xx[naive$probeID]

vector <- c()
for (probe in naive$probeID){
  if (is.na(names(gene_symbol[probe]))){
    vector <- append(vector, NA)
  }
  else{
    vector <- append(vector, gene_symbol[probe])
  }
}

df_naive<- cbind(naive, unlist(vector))
df_naive <- df_naive[df_naive$p_SMR < 10e-5, ]
df_naive <- df_naive[order(df_naive$p_SMR), ]

write.table(df_naive, file = "annot_naive_immuno_imp.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep ="\t")


########################################################################################

LPS2 <- read.table("SMR_Immuno_Imp/smr_LPS2_immuno.smr", header = TRUE, sep ="\t")

LPS2[order(LPS2$p_SMR), ] #order by p-value 

library(illuminaHumanv4.db)

x <- illuminaHumanv4SYMBOLREANNOTATED
# Get the probe identifiers that are mapped to a gene symbol
mapped_probes <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_probes])

gene_symbol <- xx[LPS2$probeID]

vector <- c()
for (probe in LPS2$probeID){
  if (is.na(names(gene_symbol[probe]))){
    vector <- append(vector, NA)
  }
  else{
    vector <- append(vector, gene_symbol[probe])
  }
}

df_LPS2<- cbind(LPS2, unlist(vector))
df_LPS2 <- df_LPS2[df_LPS2$p_SMR < 10e-5, ]
df_LPS2 <- df_LPS2[order(df_LPS2$p_SMR), ]


write.table(df_LPS2, file = "annot_LPS2_immuno_imp.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep ="\t")


########################################################################################
LPS24 <- read.table("SMR_Immuno_Imp/smr_LPS24_immuno.smr", header = TRUE, sep ="\t")

LPS24[order(LPS24$p_SMR), ] #order by p-value 

library(illuminaHumanv4.db)

x <- illuminaHumanv4SYMBOLREANNOTATED
# Get the probe identifiers that are mapped to a gene symbol
mapped_probes <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_probes])

gene_symbol <- xx[LPS24$probeID]

vector <- c()
for (probe in LPS24$probeID){
  if (is.na(names(gene_symbol[probe]))){
    vector <- append(vector, NA)
  }
  else{
    vector <- append(vector, gene_symbol[probe])
  }
}

df_LPS24<- cbind(LPS24, unlist(vector))
df_LPS24 <- df_LPS24[df_LPS24$p_SMR < 10e-5, ]
df_LPS24 <- df_LPS24[order(df_LPS24$p_SMR), ]

write.table(df_LPS24, file = "annot_LPS24_immuno_imp.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep ="\t")


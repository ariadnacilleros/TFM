#LPS 2H
##########################
#
#STEP1: REMOVE SAMPLES
#
##########################
library(readr)
LPS2 <- read_delim("LPS/LPS2.47231.261.b.txt", 
                    "\t", escape_double = FALSE, trim_ws = TRUE)

miss_pheno <- setdiff(seq(1, 432), as.numeric(colnames(LPS2)[2:262])) #samples missing in phenotype
remov_geno <- setdiff(miss_pheno, c(14, 82, 111, 165, 215, 291, 320, 375, 386, 389)) #samples to be removed from genotype

#write .txt for PLINK
write.table(x = cbind(remov_geno, remov_geno), file = "remove_pheno.txt", quote = FALSE,sep = "\t", row.names = FALSE, col.names = FALSE )

#From phenotype we only need to remove the sample 111 (column 34)
LPS2_clean <- subset(LPS2, select = -c(34)) 
dim(LPS2_clean) #260 samples + column of probes


##########################
#
#STEP2: REMOVE PROBES
#
##########################

library(illuminaHumanv4.db)
library(stringr)

################################# By genomic location ####################################
x <- illuminaHumanv4GENOMICLOCATION 
mapped_probes <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_probes])

llista<- list() #store +1 locus
llista2<- list() #store sex chr
llista3 <- list() #store non chr

for (probe in LPS2_clean$PROBE_ID){
  if (identical(grep(",",xx[probe]), integer(0)) == FALSE){ #if they have more than one locus
    llista <- append(x = llista, values = probe)
  }
  else{
    chr <- (str_match(xx[probe], "chr(.*?):"))[,2]
    if (is.na(chr) == TRUE ){ #if there is no information about chr
      llista3 <- append(x=llista3, values = probe)
    }
    else if (chr %in% c("X" , "Y")){ #if the probe is in sexual chr
      llista2 <- append(x = llista2, values = probe)
    }
  }
}

length(llista) #3223 probes mapped in more than one position
length(llista2) #2169 probes in sex chr
length(llista3) #350 probes without chr

list_positions <- unlist(llista)
list_sex <- unlist(llista2)
list_non_chr <- unlist(llista3)


################################# By SNP presence ####################################
x <- illuminaHumanv4OVERLAPPINGSNP
mapped_probes <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_probes])

SNP <- colnames(data.frame(xx))

llista_SNP <- intersect(LPS2_clean$PROBE_ID, SNP) #14.899 probes with SNPs


################################# By quality ####################################
x <- illuminaHumanv4PROBEQUALITY
mapped_probes <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_probes])

No_match <- colnames(data.frame(xx[grep("No match", xx)]))

llista_no_match <- intersect(No_match, LPS2_clean$PROBE_ID) #376 probes with No match

Bad <- colnames(data.frame(xx[grep("Bad", xx)])) 

llista_bad <- intersect(Bad, LPS2_clean$PROBE_ID) #12.379 probes with Bad


################################# Disccard probes ####################################
LPS2_clean <- subset(LPS2_clean, !(PROBE_ID %in% list_positions)) #remove probes +1 locus
dim(LPS2_clean) #44008 probes

LPS2_clean <- subset(LPS2_clean, !(PROBE_ID %in% list_non_chr))
dim(LPS2_clean) #43658 probes

LPS2_clean <- subset(LPS2_clean, !(PROBE_ID %in% list_sex))
dim(LPS2_clean) #41489 probes

LPS2_clean <- subset(LPS2_clean, !(PROBE_ID %in% llista_SNP))
dim(LPS2_clean) #28029 probes

LPS2_clean <- subset(LPS2_clean, !(PROBE_ID %in% llista_no_match))
dim(LPS2_clean) #28008 probes

LPS2_clean <- subset(LPS2_clean, !(PROBE_ID %in% llista_bad))
dim(LPS2_clean) #20808 probes 


##########################
#
#STEP3: COMPARE PROBE IDs OF LAST BLAST 
#
##########################

blast_seq <- read_delim("blast_seq.txt", 
                        "\t", escape_double = FALSE, col_names = FALSE, 
                        trim_ws = TRUE, )
setdiff(LPS2_clean$PROBE_ID, blast_seq$X1)
#we can use the same results as the naive blast search

##########################
#
#STEP4: REMOVE CROSS-HYBRIDIZATION PROBES
#
##########################
cross_hyb <- read_csv("cross_hyb.txt", col_names = FALSE)

LPS2_clean <- subset(LPS2_clean, !(PROBE_ID %in% cross_hyb$X1))
dim(LPS2_clean) #14512 final probes 

##########################
#
#STEP5: ADD CHR, START & END
#
##########################

x <- illuminaHumanv4GENOMICLOCATION
# Get the probe identifiers that are mapped to any cytoband
mapped_probes <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_probes])

df_positions <- data.frame()
for (probe in LPS2_clean$PROBE_ID){
  chr <- str_split(xx[probe], ":", n = , simplify = FALSE)[[1]][1]
  start_p <- str_split(xx[probe], ":", n = , simplify = FALSE)[[1]][2]
  end_p <- start <- str_split(xx[probe], ":", n = , simplify = FALSE)[[1]][3]
  df_positions <- rbind(df_positions, cbind(chr, start_p, end_p))
}

rownames(df_positions) <- LPS2_clean$PROBE_ID

head(df_positions)


final_df <- cbind(df_positions, LPS2_clean)

colnames_df <- list()

for (i in colnames(final_df)){
  if (i != "chr" && i != "start_p" && i != "end_p" && i != "PROBE_ID"){
    colnames_df <- append(colnames_df, values = paste0(i,"_",i))
  }
  else {
    colnames_df <- append(colnames_df, values = i)
  }
}

colnames(final_df) <- unlist(colnames_df)


##########################
#
#STEP6: NON-CANONICAL CHR
#
##########################
table(final_df$chr)

library(stringr)
final_df$chr <- str_replace_all(final_df$chr, pattern =  "chr17_ctg5_hap1", replacement = "chr17")
final_df$chr <- str_replace_all(final_df$chr, pattern =  "chr19_gl000209_random", replacement = "chr19")
final_df$chr <- str_replace_all(final_df$chr, pattern =  "chr4_ctg9_hap1", replacement = "chr4")
final_df$chr <- str_replace_all(final_df$chr, pattern =  "chr4_gl000194_random", replacement = "chr4")
final_df$chr <- str_replace_all(final_df$chr, pattern =  "chr6_apd_hap1", replacement = "chr6")
final_df$chr <- str_replace_all(final_df$chr, pattern =  "chr6_cox_hap2", replacement = "chr6")
final_df$chr <- str_replace_all(final_df$chr, pattern =  "chr6_dbb_hap3", replacement = "chr6")
final_df$chr <- str_replace_all(final_df$chr, pattern =  "chr6_mann_hap4", replacement = "chr6")
final_df$chr <- str_replace_all(final_df$chr, pattern =  "chr6_mcf_hap5", replacement = "chr6")
final_df$chr <- str_replace_all(final_df$chr, pattern =  "chr6_qbl_hap6", replacement = "chr6")
final_df$chr <- str_replace_all(final_df$chr, pattern =  "chr6_ssto_hap7", replacement = "chr6")
final_df$chr <- str_replace_all(final_df$chr, pattern =  "chr7_gl000195_random", replacement = "chr7")

disscard <- c(final_df[final_df$chr == "chrUn_gl000211", ]$PROBE_ID, final_df[final_df$chr == "chrUn_gl000212", ]$PROBE_ID, final_df[final_df$chr == "chrUn_gl000223", ]$PROBE_ID) 

final_df <- subset(final_df, !(PROBE_ID %in% disscard))

final_df$chr <- substring(final_df$chr, 4)


write.table(x = final_df, file = "phenotype_LPS2.txt", sep = "\t", quote = FALSE, row.names = FALSE)


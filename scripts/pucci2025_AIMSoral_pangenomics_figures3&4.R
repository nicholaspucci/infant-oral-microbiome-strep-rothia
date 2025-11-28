#Scripts provide the input matrices for heatmaps in Figure 3 and Figure 4.
#These were further customized in iTOL: Interactive Tree of Life (Letunic and Bork (2024) Nucleic Acids Res doi: 10.1093/nar/gkae268)) 
#and Adobe Illustrator (v28.7.5, ©1984-2024 Adobe)
setwd("/Users/path/to/directory")

install.packages("tidyr")
install.packages("dplyr")
install.packages("FSA")
install.packages("purrr")
library(tidyr)
library(dplyr)
library(FSA)
library(purrr)
####################################Metabolism heatmap input matrix Figure 3 and Figure 4
modr <- read.csv("ROTHIAgenomic75-metabolism2_modules.txt", sep = '\t', header=TRUE, stringsAsFactors = FALSE)
mods <- read.csv("STREPTOgenomic75-metabolism2_modules.txt", sep = '\t', header=TRUE, stringsAsFactors = FALSE)
matchr <- read.csv("metadata_rothia_metabolism.tsv", sep = '\t', header=TRUE, stringsAsFactors = FALSE)
matchs <- read.csv("metadata_strepto_metabolism.tsv", sep = '\t', header=TRUE, stringsAsFactors = FALSE)

matchr$db_name <- matchr$name

MODS <- merge(mods,matchs,by='db_name',all.x=T)
MODR <- merge(modr,matchr,by='db_name',all.x=T)


#Module category STREPTO
MODS2 <- aggregate(MODS$pathwise_copy_number, by=list(MODS$db_name,MODS$group,MODS$completeness,MODS$module_category), FUN=sum)
colnames(MODS2) <- c('db_name','group','completeness','module_category','pathwise_copy_number')
#Module category ROTHIA
MODR2 <- aggregate(MODR$pathwise_copy_number, by=list(MODR$db_name,MODR$group,MODR$completeness,MODR$module_category), FUN=sum)
colnames(MODR2) <- c('db_name','group','completeness','module_category','pathwise_copy_number')

MODS2 <- subset(MODS2, (!grepl('Biosynthesis of other secondary metabolites|Gene set|Module set|Xenobiotics biodegradation', module_category)))
MODR2 <- subset(MODR2, (!grepl('Biosynthesis of other secondary metabolites|Gene set|Glycan metabolism|Module set', module_category)))

#normalize metabolic module pathwise copy-number by genome/MAG completeness
MODS2$normalized_path <- MODS2$pathwise_copy_number / (MODS2$completeness / 100)
MODR2$normalized_path <- MODR2$pathwise_copy_number / (MODR2$completeness / 100)

#Filter genomes/MAGs by completeness (i.e. ≥90%)
MODS3 <- MODS2[MODS2$completeness >= 90, ]
MODR3 <- MODR2[MODR2$completeness >= 90, ]


df_complete_strepto <- MODS3 %>%
  complete(db_name, module_category, 
           fill = list(pathwise_copy_number = 0, normalized_path = 0)) %>%
  group_by(db_name) %>%
  fill(group, completeness, .direction = "downup") %>%
  ungroup()

df_complete_rothia <- MODR3 %>%
  complete(db_name, module_category, 
           fill = list(pathwise_copy_number = 0, normalized_path = 0)) %>%
  group_by(db_name) %>%
  fill(group, completeness, .direction = "downup") %>%
  ungroup()

# Dunn's test for pairwise comparisons with FDR correction
df_complete_strepto$group <- as.factor(df_complete_strepto$group)
df_complete_rothia$group <- as.factor(df_complete_rothia$group)

dunn_results_strepto <- df_complete_strepto %>%
  group_by(module_category) %>%
  summarise(
    dunn_test = list(dunnTest(normalized_path ~ group, 
                              data = cur_data(), 
                              method = "bh"))
  ) %>%
  mutate(comparisons = map(dunn_test, ~.x$res)) %>%
  unnest(comparisons) %>%
  select(module_category, Comparison, Z, P.unadj, P.adj)

dunn_results_rothia <- df_complete_rothia %>%
  group_by(module_category) %>%
  summarise(
    dunn_test = list(dunnTest(normalized_path ~ group, 
                              data = cur_data(), 
                              method = "bh"))
  ) %>%
  mutate(comparisons = map(dunn_test, ~.x$res)) %>%
  unnest(comparisons) %>%
  select(module_category, Comparison, Z, P.unadj, P.adj)


##Convert long- to wide-format
MODS2reformat <- MODS2[,c(-2,-3,-5)]
MODR2reformat <- MODR2[,c(-2,-3,-5)]
matr <- spread(MODR2reformat, key = module_category, value = normalized_path)
mats <- spread(MODS2reformat, key = module_category, value = normalized_path)
rownames(matr) <- matr$db_name
rownames(mats) <- mats$db_name
matr <- matr[,c(-1)]
mats <- mats[,c(-1)]
matr[is.na(matr)] <- 0
mats[is.na(mats)] <- 0

matr <- cbind(matr, "Glycan metabolism" = 0)

#write.csv(matr, "matrix-metabolism-fig4-rothia.csv")
#write.csv(mats, "matrix-metabolism-fig3-strepto.csv")


###################################Metabolic modules presence/absence heatmap input matrix Figure 3/4
#Selection based on output metabolic enrichment analysis (focus on S. sp000187445 and R. sp902373285)
strepmod <- read.csv("STREPTOgenomic-metabolism2_matrix75-module_pathwise_presence-MATRIX.txt", sep = '\t', header=TRUE, stringsAsFactors = FALSE)
rothmod <- read.csv("ROTHIAgenomic-metabolism2_matrix75-module_pathwise_presence-MATRIX.txt", sep = '\t', header=TRUE, stringsAsFactors = FALSE)

module_strep <- c("M00010","M00020","M00023","M00026","M00028","M00052","M00053","M00140", "M00844","M00916")
module_roth <- c("M00006","M00016","M00082","M00116","M00126","M00432","M00526","M00527", "M00916","M00926","M00364")

filtstrepto <- strepmod %>%
  filter(module %in% module_strep)

filtrothia <- rothmod %>%
  filter(module %in% module_roth)

#write.csv(filtstrepto, "matrix-metabolic-modules-fig3.csv")
#write.csv(filtrothia, "matrix-metabolic-modules-fig4.csv")

###################################Streptococcus CAZymes: generation of heatmap input matrix Figure 3
strepto <- read.csv("CAZYmes-streptoccus.tsv", sep = '\t', header=TRUE, stringsAsFactors = FALSE)
misc <- read.csv("metadata_strepto_cazymes.tsv", sep = '\t', header=TRUE, stringsAsFactors = FALSE)

#normalize CAZyme-encoding gene copy-number by genome/MAG completeness
strepto$normalized_cazyme <- strepto$copy_number / (strepto$genome_completeness / 100)

#Filter genomes/MAGs by completeness (i.e. ≥90%)
strepto2 <- strepto[strepto$genome_completeness >= 90, ]

df_complete <- strepto2 %>%
  complete(bin, cazyme, 
           fill = list(copy_number = 0, normalized_cazyme = 0)) %>%
  group_by(bin) %>%
  fill(name, group, genome_completeness, .direction = "downup") %>%
  ungroup()

# Dunn's test for pairwise comparisons with FDR correction
df_complete$group <- as.factor(df_complete$group)
df_complete <- subset(df_complete, (!grepl('GH163', cazyme)))

dunn_cazymes_strepto <- df_complete %>%
  group_by(cazyme) %>%
  summarise(
    dunn_test = list(dunnTest(normalized_cazyme ~ group, 
                              data = cur_data(), 
                              method = "bh"))
  ) %>%
  mutate(comparisons = map(dunn_test, ~.x$res)) %>%
  unnest(comparisons) %>%
  select(cazyme, Comparison, Z, P.unadj, P.adj)


#Filter CAZymes for figure 3 heatmap
#A small selection CAZyme-encoding genes found to be significantly more abundant in 
#S. sp000187445 vs. S. salivarius and S. sp000187445 and S. salivarius vs. other 
#infant-related Streptococcus species are selected
list_mat <- c("GH36","GH66","GH68","GH8","GH87","GT2_Glyco_trans_2_3","GT45","GT81","GT4","GT35")
longfig <- strepto %>%
  filter(cazyme %in% list_mat)

mat_input<- longfig[,c(-3,-6)]
CAZYm<- spread(mat_input, key = cazyme, value = normalized_cazyme)

CAZYm[is.na(CAZYm)] <- 0
CAZYm <- as.matrix(CAZYm)
rownames(CAZYm) <- CAZYm[,1]
CAZYm <- CAZYm[,-1:-3]
CAZYm <- as.data.frame(CAZYm)

#write.csv(CAZYm, "matrix-cazymes-fig3-heatmap.csv")


#Rothia CAZymes: generation of heatmap input matrix Figure 4

rothia <- read.csv("CAZYmes-rothia.tsv", sep = '\t', header=TRUE, stringsAsFactors = FALSE)
misc <- read.csv("metadata_rothia_cazymes.tsv", sep = '\t', header=TRUE, stringsAsFactors = FALSE)

#normalize CAZyme-encoding gene copy-number by genome/MAG completeness
rothia$normalized_cazyme <- rothia$copy_number / (rothia$genome_completeness / 100)

#Filter genomes/MAGs by completeness (i.e. ≥90%)
rothia2 <- rothia[rothia$genome_completeness >= 90, ]

df_complete <- rothia2 %>%
  complete(bin, cazyme, 
           fill = list(copy_number = 0, normalized_cazyme = 0)) %>%
  group_by(bin) %>%
  fill(name, group, genome_completeness, .direction = "downup") %>%
  ungroup()

# Dunn's test for pairwise comparisons with FDR correction
df_complete$group <- as.factor(df_complete$group)
dunn_cazymes_rothia <- df_complete %>%
  group_by(cazyme) %>%
  summarise(
    dunn_test = list(dunnTest(normalized_cazyme ~ group, 
                              data = cur_data(), 
                              method = "bh"))
  ) %>%
  mutate(comparisons = map(dunn_test, ~.x$res)) %>%
  unnest(comparisons) %>%
  select(cazyme, Comparison, Z, P.unadj, P.adj)

#write.csv(dunn_results, "rothia_CAZYME_updated.csv")

#Filter CAZymes for figure 3 heatmap
#A small selection CAZyme-encoding genes found to be significantly more abundant in 
#R. sp000187445 vs. . salivarius and S. sp000187445 and S. salivarius vs. other 
#infant-related Streptococcus species are selected
list_mat <- c("AA1","AA1_1","GH13_8","GT11","GT2_Glyco_tranf_2_2","GT2_Glyco_tranf_2_4","GT25","GT4", "GT76")

longfig <- rothia %>%
  filter(cazyme %in% list_mat)

mat_input<- longfig[,c(-3,-6)]
CAZYm<- spread(mat_input, key = cazyme, value = normalized_cazyme)

CAZYm[is.na(CAZYm)] <- 0
CAZYm <- as.matrix(CAZYm)
rownames(CAZYm) <- CAZYm[,1]
CAZYm <- CAZYm[,-1:-3]
CAZYm <- as.data.frame(CAZYm)

#write.csv(CAZYm, "matrix-cazymes-fig4-heatmap.csv")




####################################
#Enterosalivary pathway: normalized gene copy-number heatmap Figure 4
mat <- read.csv("ROTHIAgenomic-metabolism2_matrix75-enzyme_hits-MATRIX.txt", sep = '\t', header=TRUE, stringsAsFactors = FALSE)
ko <- read.csv("KOlist-NOSreduction.txt", sep = '\t', header=TRUE, stringsAsFactors = FALSE)
misc <- read.csv("CAZYmes-rothia.tsv", sep = '\t', header=TRUE, stringsAsFactors = FALSE)

misc <- misc[,c(1,3,4,6)]
misc <- unique(misc)
mat2 <- mat %>%
  filter(enzyme %in% ko$KO)

mat2$KO <- mat2$enzyme
mat2$enzyme <- NULL

mat3 <- merge(ko,mat2,by='KO',all.x=T)
mat3 <- na.omit(mat3)
mat3 <- mat3[,c(-1,-3,-4)]

mat3 <- as.data.frame(t(mat3))
colnames(mat3) <- mat3[1,]
mat3 <- mat3[-1,]
mat3$name <- rownames(mat3)
long <- gather(mat3, enzymes, copy_number, "nirB":"narK", factor_key = T)
long <- unique(long)

LONG <- merge(long, misc, by='name', all.x=T)
LONG$copy_number <- as.numeric(LONG$copy_number)
LONG$normalized_enzyme <- LONG$copy_number / (LONG$genome_completeness / 100)
#Filter genomes/MAGs by completeness (i.e. ≥90%)
LONG2 <- LONG[LONG$genome_completeness >= 90, ]

df_complete_entsalv <- LONG2 %>%
  complete(bin, enzymes, 
           fill = list(copy_number = 0, normalized_enzyme = 0)) %>%
  group_by(bin) %>%
  fill(name, group, genome_completeness, .direction = "downup") %>%
  ungroup()

# Dunn's test for pairwise comparisons with FDR correction
df_complete_entsalv$group <- as.factor(df_complete_entsalv$group)
df_complete_entsalv$copy_number <- as.numeric(df_complete_entsalv$copy_number)

dunn_entsalv_rothia <- df_complete_entsalv %>%
  dplyr::group_by(enzymes) %>%
  dplyr::summarise(
    dunn_test = list(dunnTest(normalized_enzyme ~ group, 
                              data = cur_data(), 
                              method = "bh"))
  ) %>%
  mutate(comparisons = map(dunn_test, ~.x$res)) %>%
  unnest(comparisons) %>%
  select(enzymes, Comparison, Z, P.unadj, P.adj)


LONG <- LONG[,c(-3,-4,-6)]
NOmat<- spread(LONG, key = enzymes, value = normalized_enzyme)
NOmat[is.na(NOmat)] <- 0
NOmat <- as.matrix(NOmat)
rownames(NOmat) <- NOmat[,1]
NOmat <- as.data.frame(NOmat)

#write.csv(NOmat,"matrix-NOpathway-fig4-heatmap_check.csv")





################################
#Lectins/adhesins: normalized gene copy-number heatmap Figure 3
lectin <- read.csv("lectins_strepto.tsv", sep = '\t', header=TRUE, stringsAsFactors = FALSE)

lectin$normalized_count <- lectin$copy_number / (lectin$genome_completeness / 100)
lectin <- lectin[lectin$genome_completeness >= 90, ]

LEC<- spread(lectin, key = description, value = normalized_count)
LEC[is.na(LEC)] <- 0
LEC <- as.matrix(LEC)
rownames(LEC) <- LEC[,1]
LEC <- as.data.frame(LEC)

lectin_complete <- gather(LEC, description,normalized_count, "accessory Sec-dependent serine-rich glycoprotein adhesin, partial":"YSIRK-type signal peptide-containing protein", factor_key = T)

# Dunn's test for pairwise comparisons with FDR correction
lectin_complete$group <- as.factor(lectin_complete$group)
lectin_complete$normalized_count <- as.numeric(lectin_complete$normalized_count)
lectin_complete$copy_number <- as.numeric(lectin_complete$copy_number)

dunn_results <- lectin_complete %>%
  group_by(description) %>%
  summarise(
    dunn_test = list(dunnTest(normalized_count ~ group, 
                              data = cur_data(), 
                              method = "bh"))
  ) %>%
  mutate(comparisons = map(dunn_test, ~.x$res)) %>%
  unnest(comparisons) %>%
  select(description, Comparison, Z, P.unadj, P.adj)


LEC2 <- LEC[,c(-1:-5)]

#write.csv(LECstats2,"adhesins-fig3-heatmap.csv")


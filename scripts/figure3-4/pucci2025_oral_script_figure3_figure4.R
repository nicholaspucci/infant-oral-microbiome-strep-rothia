#Scripts provide the input matrices for heatmaps in Figure 3 and Figure 4.
#These were further customized in iTOL: Interactive Tree of Life (Letunic and Bork (2024) Nucleic Acids Res doi: 10.1093/nar/gkae268)) 
#and Adobe Illustrator (v28.7.5, ©1984-2024 Adobe)
setwd("/Users/path/to/directory")

library(tidyr)
library(dplyr)
####################################Metabolism heatmap input matrix Figure 3 and Figure 4
modr <- read.csv("ROTHIAgenomic75-metabolism2_modules.txt", sep = '\t', header=TRUE, stringsAsFactors = FALSE)
mods <- read.csv("STREPTOgenomic75-metabolism2_modules.txt", sep = '\t', header=TRUE, stringsAsFactors = FALSE)
matchr <- read.csv("rothia-isolates-match-cazymes.tsv", sep = '\t', header=TRUE, stringsAsFactors = FALSE)
matchs <- read.csv("strepto-isolates-match-cazymes.tsv", sep = '\t', header=TRUE, stringsAsFactors = FALSE)

matchs$db_name <- matchs$name
matchr$db_name <- matchr$name

MODS <- merge(mods,matchs,by='db_name',all.x=T)
MODR <- merge(modr,matchr,by='db_name',all.x=T)


#Module category STREPTO
MODS2 <- aggregate(MODS$pathwise_copy_number, by=list(MODS$db_name,MODS$group,MODS$module_category), FUN=sum)
colnames(MODS2) <- c('db_name','group','module_category','pathwise_copy_number')
#Module category ROTHIA
MODR2 <- aggregate(MODR$pathwise_copy_number, by=list(MODR$db_name,MODR$group,MODR$module_category), FUN=sum)
colnames(MODR2) <- c('db_name','group','module_category','pathwise_copy_number')

matr <- spread(MODR2, key = module_category, value = pathwise_copy_number)
mats <- spread(MODS2, key = module_category, value = pathwise_copy_number)
rownames(matr) <- matr$db_name
rownames(mats) <- mats$db_name
matr <- matr[,c(-1:-2)]
mats <- mats[,c(-1:-2)]
matr[is.na(matr)] <- 0
mats[is.na(mats)] <- 0

matr <- matr[,c(-2,-6,-10)]

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
misc <- read.csv("strepto-isolates-match-cazymes.tsv", sep = '\t', header=TRUE, stringsAsFactors = FALSE)

#normalize CAZyme-encoding gene copy-number by genome/MAG completeness
strepto$normalized_cazyme <- strepto$copy_number / (strepto$genome_completeness / 100)

#Filter genomes/MAGs by completeness (i.e. ≥90%)
strepto2 <- strepto[strepto$genome_completeness >= 90, ]

#Carry out pairwise wilcoxon tests to test differences in CAZyme-encoding gene copy-number among infant-related species
phoc <- pairwise.wilcox.test(strepto2$normalized_cazyme, interaction(strepto2$group, strepto2$cazyme), p.adjust.method = "fdr")
phoc_p <- as.data.frame(phoc$p.value)
phoc_p$column1 <- rownames(phoc_p)

#re-order dataframe
longphoc <- gather(phoc_p, column2,p_val, "LACTARIUS.CE1":"SALIVARIUS.PL8_3", factor_key = T)

longphoc2 <- longphoc %>%
  mutate(
    group1 = sub("\\..*", "", `column1`),
    group2 = sub("\\..*", "", `column2`)
  )

longphoc2 <- longphoc2 %>%
  mutate(
    cazyme1 = sub(".*\\.", "", column1),
    cazyme2 = sub(".*\\.", "", column2)
  )

longphoc3 <- longphoc2[longphoc2$group1 == "Sp000187445" & longphoc2$group2 == "SALIVARIUS", ]

longphoc3 <- longphoc3 %>%
  drop_na()

longphoc3 <- longphoc3 %>%
  filter(cazyme1 == cazyme2)

longphoc4 <- longphoc3 %>%
  filter(p_val < 0.05)

longphoc5 <- longphoc2 %>%
  filter(p_val < 0.05)

longphoc5 <- longphoc5 %>%
  filter(cazyme1 == cazyme2)



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
misc <- read.csv("rothia-isolates-match-cazymes.tsv", sep = '\t', header=TRUE, stringsAsFactors = FALSE)

#normalize CAZyme-encoding gene copy-number by genome/MAG completeness
rothia$normalized_cazyme <- rothia$copy_number / (rothia$genome_completeness / 100)

#Filter genomes/MAGs by completeness (i.e. ≥90%)
rothia2 <- rothia[rothia$genome_completeness >= 90, ]

#Carry out pairwise wilcoxon tests to test differences in CAZyme-encoding gene copy-number among infant-related species
phoc <- pairwise.wilcox.test(rothia2$normalized_cazyme, interaction(rothia2$group, rothia2$cazyme), p.adjust.method = "fdr")
phoc_p <- as.data.frame(phoc$p.value)
phoc_p$column1 <- rownames(phoc_p)

longphoc <- gather(phoc_p, column2,p_val, "R. mucilaginosa.AA1":"R. mucilaginosa_A.SLH", factor_key = T)

longphoc_rmuci <- subset(longphoc, (!grepl('sp9023', column1)))
longphoc_rmuci <- subset(longphoc_rmuci, (!grepl('sp9023', column2)))
longphoc_sp <- subset(longphoc, (grepl('sp9023', column1)))
longphoc_sp <- subset(longphoc_sp, (grepl('sp9023', column2)))

longphoc_rmuci$cazyme1 <- sub("^[^.]*\\.[^.]*\\.", "", longphoc_rmuci$column1)
longphoc_rmuci$cazyme2 <- sub("^[^.]*\\.[^.]*\\.", "", longphoc_rmuci$column2)
longphoc_sp$cazyme1 <- sub("^[^.]*\\.[^.]*\\.", "", longphoc_sp$column1)
longphoc_sp$cazyme2 <- sub("^[^.]*\\.[^.]*\\.", "", longphoc_sp$column2)

# Extract everything before the second '.' for group1 and group2
longphoc_rmuci$group1 <- sub("\\.[^.]*$", "", longphoc_rmuci$column1)  # Remove everything after the last '.'
longphoc_rmuci$group2 <- sub("\\.[^.]*$", "", longphoc_rmuci$column2)  # Remove everything after the last '.'

longphoc_sp <- longphoc_sp %>%
  mutate(
    group1 = sub("\\..*", "", `column1`),
    group2 = sub("\\..*", "", `column2`)
  )


longphoc2 <- rbind(longphoc_rmuci, longphoc_sp)
## across all species
longphoc3 <- longphoc2 %>%
  drop_na()

longphoc3 <- longphoc3 %>%
  filter(cazyme1 == cazyme2)

longphoc4 <- longphoc3 %>%
  filter(p_val < 0.05)


#Filter CAZymes for figure 3 heatmap
#A small selection CAZyme-encoding genes found to be significantly more abundant in 
#S. sp000187445 vs. S. salivarius and S. sp000187445 and S. salivarius vs. other 
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
LONG <- LONG[,c(-3:-6)]

NOmat<- spread(LONG, key = enzymes, value = normalized_enzyme)
NOmat[is.na(NOmat)] <- 0
NOmat <- as.matrix(NOmat)
rownames(NOmat) <- NOmat[,1]
NOmat <- as.data.frame(NOmat)

#write.csv(NOmat,"matrix-NOpathway-fig4-heatmap.csv")

phoc <- pairwise.wilcox.test(LONG$normalized_enzyme, interaction(LONG$group, LONG$enzymes), p.adjust.method = "fdr")
phoc_p <- as.data.frame(phoc$p.value)
phoc_p$column1 <- rownames(phoc_p)

longphoc <- gather(phoc_p, column2,p_val, "R. mucilaginosa.nirB":"R. mucilaginosa_A.narK", factor_key = T)


longphoc <- longphoc %>%
  mutate(
    no1 = sub(".*\\.", "", column1),
    no2 = sub(".*\\.", "", column2)
  )

longphoc2 <- longphoc %>%
  drop_na()
longphoc2 <- longphoc2 %>%
  filter(no1 == no2)

result <- LONG %>%
  group_by(group, enzymes) %>%
  summarize(meanARG = mean(normalized_enzyme, na.rm = TRUE),
            standard_error = sd(normalized_enzyme, na.rm = TRUE) / sqrt(n()))



ggplot(result, aes(x = reorder(interaction(enzymes,group), -meanARG), y=meanARG, fill=group)) + 
  geom_bar(position='dodge',stat = 'identity', col='black', size=0.1)+
  geom_errorbar(aes(ymin=meanARG-standard_error, ymax=meanARG+standard_error), width=0.9, position=position_dodge2(), size=0.1)+
  theme_classic()+
  theme(legend.title=element_blank(),axis.title.x=element_blank(),
        panel.spacing.x = unit(0,"lines"))+
  theme(legend.position="bottom")+
  facet_wrap(~enzymes, scales = 'free_x')+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=4))






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
LEC2 <- LEC[,c(-1:-5)]

#write.csv(LECstats2,"adhesins-fig3-heatmap.csv")



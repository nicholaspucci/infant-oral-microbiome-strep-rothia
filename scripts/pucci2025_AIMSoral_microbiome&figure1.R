setwd("/Users/nicholaspucci/Documents/AUMC")
library(tidyr)
library(magrittr)
library(dplyr)
library(corrplot)
library(vegan)
library(ggplot2)
library(readxl)

gtdbclade <- read.csv("gtdb_genome_clade_id.csv", sep = ',', header=TRUE, stringsAsFactors = FALSE)
oral <- read.csv("aims-oral-sylph.tsv", sep = '\t', header=TRUE, stringsAsFactors = FALSE)
taxa <- read.csv("aims-oral-r220.tsv", sep = '\t', header=TRUE, stringsAsFactors = FALSE)
meta <- read.csv("metadata_Pucci2025.tsv", sep = '\t', header=TRUE, stringsAsFactors = FALSE)

gtdbclade <- gtdbclade[,-1]

oral2 <- merge(oral,gtdbclade,by='Genome_file',all.x=T)

ORAL <- merge(oral2, meta, by='Sample_file', all.x = T)


#Make FIGURE 1A (further customed in Adobe Illustrator)
#Create column containing genus annotations for each species
ORAL$genus <- sapply(ORAL$clade_name, function(x) {
  if (grepl("g__", x)) {
    sub(".*(g__[^|]*).*", "\\1", x)
  } else if (grepl("undescribed", x, ignore.case = TRUE)) {
    x
  } else {
    NA
  }
})


#For each timepoint and sample type, re-order the dataframe to obtain the relative abundances 
#of all bacterial species per subject (i.e. 'Pair' column) 
mom <- ORAL %>%
  filter(grepl("pregnancy", Measurement, ignore.case = TRUE) & 
           grepl("tongue swab", sample_type, ignore.case = TRUE))

momt <- ORAL %>%
  filter(grepl("pregnancy", Measurement, ignore.case = TRUE) & 
           grepl("tooth plaque", sample_type, ignore.case = TRUE))

tongue1 <- ORAL %>%
  filter(grepl("1 month", Measurement, ignore.case = TRUE) & 
           grepl("tongue swab", sample_type, ignore.case = TRUE))

tongue6 <- ORAL %>%
  filter(grepl("6 months", Measurement, ignore.case = TRUE) & 
           grepl("tongue swab", sample_type, ignore.case = TRUE))

tooth6 <- ORAL %>%
  filter(grepl("6 months", Measurement, ignore.case = TRUE) & 
           grepl("tooth plaque", sample_type, ignore.case = TRUE))

tonguemagg <- aggregate(mom$Taxonomic_abundance, by=list(mom$genus,mom$Pair,mom$Sample_file), FUN=sum)
colnames(tonguemagg) <- c('genus','Pair','Sample_file','Taxonomic_abundance')
tonguemtagg <- aggregate(momt$Taxonomic_abundance, by=list(momt$genus,momt$Pair,momt$Sample_file), FUN=sum)
colnames(tonguemtagg) <- c('genus','Pair','Sample_file','Taxonomic_abundance')
tongue1agg <- aggregate(tongue1$Taxonomic_abundance, by=list(tongue1$genus,tongue1$Pair,tongue1$Sample_file), FUN=sum)
colnames(tongue1agg) <- c('genus','Pair','Sample_file','Taxonomic_abundance')
tongue6agg <- aggregate(tongue6$Taxonomic_abundance, by=list(tongue6$genus,tongue6$Pair,tongue6$Sample_file), FUN=sum)
colnames(tongue6agg) <- c('genus','Pair','Sample_file','Taxonomic_abundance')
tooth6agg <- aggregate(tooth6$Taxonomic_abundance, by=list(tooth6$genus,tooth6$Pair,tooth6$Sample_file), FUN=sum)
colnames(tooth6agg) <- c('genus','Pair','Sample_file','Taxonomic_abundance')

tonguemwide <- pivot_wider(data = tonguemagg, 
                           id_cols = Pair,  # Specify the columns to keep as identifiers
                           names_from = genus,  # Specify the column containing the new column names
                           values_from = Taxonomic_abundance) 
tonguemwide <- as.matrix(tonguemwide)
tonguemwide[is.na(tonguemwide)] <- 0
rownames(tonguemwide) <- tonguemwide[,1]
tonguemwide <- as.data.frame(tonguemwide)
tonguemwide$Pair <- NULL

tonguemwide[,1:104] %<>% lapply(function(x) as.numeric(as.character(x)))
tonguemwide <- tonguemwide[,1:104] %>% dplyr::select(which(!colSums(tonguemwide[,1:104], na.rm=TRUE) %in% 0))
tonguemwide <- tonguemwide/100


tonguemtwide <- pivot_wider(data = tonguemtagg, 
                            id_cols = Pair,  # Specify the columns to keep as identifiers
                            names_from = genus,  # Specify the column containing the new column names
                            values_from = Taxonomic_abundance) 
tonguemtwide <- as.matrix(tonguemtwide)
tonguemtwide[is.na(tonguemtwide)] <- 0
rownames(tonguemtwide) <- tonguemtwide[,1]
tonguemtwide <- as.data.frame(tonguemtwide)
tonguemtwide$Pair <- NULL

tonguemtwide[,1:109] %<>% lapply(function(x) as.numeric(as.character(x)))
tonguemtwide <- tonguemtwide[,1:109] %>% dplyr::select(which(!colSums(tonguemtwide[,1:109], na.rm=TRUE) %in% 0))
tonguemtwide <- tonguemtwide/100

tongue1wide <- pivot_wider(data = tongue1agg, 
                           id_cols = Pair,  # Specify the columns to keep as identifiers
                           names_from = genus,  # Specify the column containing the new column names
                           values_from = Taxonomic_abundance) 
tongue1wide <- as.matrix(tongue1wide)
tongue1wide[is.na(tongue1wide)] <- 0
rownames(tongue1wide) <- tongue1wide[,1]
tongue1wide <- as.data.frame(tongue1wide)
tongue1wide$Pair <- NULL

library(dplyr)
tongue1wide[,1:26] %<>% lapply(function(x) as.numeric(as.character(x)))
tongue1wide <- tongue1wide[,1:26] %>% dplyr::select(which(!colSums(tongue1wide[,1:26], na.rm=TRUE) %in% 0))
tongue1wide <- tongue1wide/100

tongue6wide <- pivot_wider(data = tongue6agg, 
                           id_cols = Pair,  # Specify the columns to keep as identifiers
                           names_from = genus,  # Specify the column containing the new column names
                           values_from = Taxonomic_abundance) 
tongue6wide <- as.matrix(tongue6wide)
tongue6wide[is.na(tongue6wide)] <- 0
rownames(tongue6wide) <- tongue6wide[,1]
tongue6wide <- as.data.frame(tongue6wide)
tongue6wide$Pair <- NULL

library(dplyr)
tongue6wide[,1:53] %<>% lapply(function(x) as.numeric(as.character(x)))
tongue6wide <- tongue6wide[,1:53] %>% dplyr::select(which(!colSums(tongue6wide[,1:53], na.rm=TRUE) %in% 0))
tongue6wide <- tongue6wide/100

tooth6wide <- pivot_wider(data = tooth6agg, 
                          id_cols = Pair,  # Specify the columns to keep as identifiers
                          names_from = genus,  # Specify the column containing the new column names
                          values_from = Taxonomic_abundance) 
tooth6wide <- as.matrix(tooth6wide)
tooth6wide[is.na(tooth6wide)] <- 0
rownames(tooth6wide) <- tooth6wide[,1]
tooth6wide <- as.data.frame(tooth6wide)
tooth6wide$Pair <- NULL

library(dplyr)
tooth6wide[,1:71] %<>% lapply(function(x) as.numeric(as.character(x)))
tooth6wide <- tooth6wide[,1:71] %>% dplyr::select(which(!colSums(tooth6wide[,1:71], na.rm=TRUE) %in% 0))
tooth6wide <- tooth6wide/100


tonguemwide$Pair <- rownames(tonguemwide)
tonguemtwide$Pair <- rownames(tonguemtwide)
tongue1wide$Pair <- rownames(tongue1wide)
tongue6wide$Pair <- rownames(tongue6wide)
tooth6wide$Pair <- rownames(tooth6wide)

tonguemlong <- tonguemwide%>%
  pivot_longer(
    cols = c(,1:104), # columns to pivot, excluding the identifier column(s)
    names_to = "genus", # name of the new column that will contain the names of the original columns
    values_to = "Taxonomic_abundance" # name of the new column that will contain the values of the original columns
  )

tonguemtlong <- tonguemtwide%>%
  pivot_longer(
    cols = c(,1:109), # columns to pivot, excluding the identifier column(s)
    names_to = "genus", # name of the new column that will contain the names of the original columns
    values_to = "Taxonomic_abundance" # name of the new column that will contain the values of the original columns
  )

tongue1long <- tongue1wide%>%
  pivot_longer(
    cols = c(,1:26), # columns to pivot, excluding the identifier column(s)
    names_to = "genus", # name of the new column that will contain the names of the original columns
    values_to = "Taxonomic_abundance" # name of the new column that will contain the values of the original columns
  )

tongue6long <- tongue6wide%>%
  pivot_longer(
    cols = c(,1:53), # columns to pivot, excluding the identifier column(s)
    names_to = "genus", # name of the new column that will contain the names of the original columns
    values_to = "Taxonomic_abundance" # name of the new column that will contain the values of the original columns
  )

tooth6long <- tooth6wide%>%
  pivot_longer(
    cols = c(,1:71), # columns to pivot, excluding the identifier column(s)
    names_to = "genus", # name of the new column that will contain the names of the original columns
    values_to = "Taxonomic_abundance" # name of the new column that will contain the values of the original columns
  )



#For each timepoint and sample type, calculate the average abundance of each genus
averageTM <- tonguemlong %>%
  group_by(genus) %>%
  dplyr::summarize(mean_abundance = mean(Taxonomic_abundance),
                   standard_error = sd(Taxonomic_abundance, na.rm = TRUE) / sqrt(n()))

averageTOM <- tonguemtlong %>%
  group_by(genus) %>%
  dplyr::summarize(mean_abundance = mean(Taxonomic_abundance),
                   standard_error = sd(Taxonomic_abundance, na.rm = TRUE) / sqrt(n()))

averageT1 <- tongue1long %>%
  group_by(genus) %>%
  dplyr::summarize(mean_abundance = mean(Taxonomic_abundance),
                   standard_error = sd(Taxonomic_abundance, na.rm = TRUE) / sqrt(n()))
averageT6 <- tongue6long %>%
  group_by(genus) %>%
  dplyr::summarize(mean_abundance = mean(Taxonomic_abundance),
                   standard_error = sd(Taxonomic_abundance, na.rm = TRUE) / sqrt(n()))
averageTO6 <- tooth6long %>%
  group_by(genus) %>%
  dplyr::summarize(mean_abundance = mean(Taxonomic_abundance),
                   standard_error = sd(Taxonomic_abundance, na.rm = TRUE) / sqrt(n()))


#For each timepoint and sample type, identify the 5 most abundant genera
df_top5tonguem <- averageTM %>%
  slice_max(mean_abundance, n = 5) %>%       # Retain the top 5 highest values in column D
  ungroup()

df_top5tonguemt <- averageTOM %>%
  slice_max(mean_abundance, n = 5) %>%       # Retain the top 5 highest values in column D
  ungroup()

df_top5tongue1 <- averageT1 %>%
  slice_max(mean_abundance, n = 5) %>%       # Retain the top 5 highest values in column D
  ungroup()

df_top5tongue6 <- averageT6 %>%
  slice_max(mean_abundance, n = 5) %>%       # Retain the top 5 highest values in column D
  ungroup()

df_top5tooth6 <- averageTO6 %>%
  slice_max(mean_abundance, n = 5) %>%       # Retain the top 5 highest values in column D
  ungroup()


#For each timepoint and sample type, calculate the relative abundance of the remaining genera ('Other') 
#as 100 minus the total
df_totalm <- sum(df_top5tonguem$mean_abundance)
other_abundancem <- 1 - df_totalm
new_rowm <- c("Other bacteria genus",0.310928363636364,NA) 
testm <- rbind(df_top5tonguem,new_rowm)

df_totalmt <- sum(df_top5tonguemt$mean_abundance)
other_abundanceTOm <- 1 - df_totalmt
new_rowmt <- c("Other bacteria genus",0.234536333333333,NA) 
testmt <- rbind(df_top5tonguemt,new_rowmt)

df_total1 <- sum(df_top5tongue1$mean_abundance)
other_abundance1 <- 1 - df_total1
new_row1 <- c("Other bacteria genus",0.182370444444444,NA) 
test1 <- rbind(df_top5tongue1,new_row1)

df_total6 <- sum(df_top5tongue6$mean_abundance)
other_abundance6 <- 1 - df_total6
new_row6 <- c("Other bacteria genus",0.1732249444444,NA) 
test6 <- rbind(df_top5tongue6,new_row6)

df_totalTO6 <- sum(df_top5tooth6$mean_abundance)
other_abundanceTO6 <- 1 - df_totalTO6
new_rowTO6 <- c("Other bacteria genus",0.1873585,NA) 
testTO6 <- rbind(df_top5tooth6,new_rowTO6)

#For each timepoint and sampletype, add a column with sample type identity
testm$Oral_niche <- as.factor("Tongue dorsum")
testmt$Oral_niche <- as.factor("Dental plaque")
test1$Oral_niche <- as.factor("Tongue dorsum")
test6$Oral_niche <- as.factor("Tongue dorsum")
testTO6$Oral_niche <- as.factor("Dental plaque")

#For each timepoint and sampletype, add a column with sample type timepoint identity
testm$Measurement_point <- as.factor("Mother (34wks pregnancy)")
testmt$Measurement_point <- as.factor("Mother (34wks pregnancy)")
test1$Measurement_point <- as.factor("Infant (1 month)")
test6$Measurement_point <- as.factor("Infant (6 months)")
testTO6$Measurement_point <- as.factor("Infant (6 months)")

#merge all dataframes together
GENUSM <- rbind(testm,testmt)
GENUS1 <- rbind(GENUSM,test1)
GENUS2 <- rbind(GENUS1,test6)
GENUS3 <- rbind(GENUS2,testTO6)

GENUS3$mean_abundance <- as.numeric(GENUS3$mean_abundance)
GENUS3$mean_abundance2 <- GENUS3$mean_abundance*100

#manually assign colors to genus
genus_colors <- c("g__Streptococcus" = "#FFC107", 
                  "g__Rothia" ='#D81B60',
                  "g__Veillonella"='#1E88E5',
                  "g__Prevotella"='#004D40',
                  "g__Pauljensenia"='#FEC181',
                  "g__Neisseria" = "#2C8D40", 
                  "g__Lactobacillus" = "#7A55D7", 
                  "g__Actinomyces" = "#5E9FD0",
                  "Other bacteria genus"='light gray'
)


ggplot(GENUS3, aes(x=Measurement_point, y=mean_abundance2, fill=genus)) +
  geom_bar(position="stack", stat="identity")+theme_classic()+
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5, size = 7))+
  ylab("Relative abundance (%)")+
  facet_wrap(~Oral_niche)+
  scale_fill_manual(values = genus_colors)




###Make FIGURE 1B (further customed in Adobe Illustrator)
ORAL2 <- ORAL %>%
  filter(Taxonomic_abundance > 0.1)

ORAL2$genus <- sapply(ORAL2$clade_name, function(x) {
  if (grepl("g__", x)) {
    sub(".*(g__[^|]*).*", "\\1", x)
  } else if (grepl("undescribed", x, ignore.case = TRUE)) {
    x
  } else {
    NA
  }
})

ORAL2$species <- sapply(ORAL2$clade_name, function(x) {
  # Check if the clade_name contains 's__'
  if (grepl("s__", x)) {
    # Extract substring starting from 's__'
    sub(".*(s__[^|]*).*", "\\1", x)
  } else if (grepl("undescribed", x, ignore.case = TRUE)) {
    # If 'undescribed' is found, retain the entire string
    x
  } else {
    # Otherwise, return NA or empty string
    NA
  }
})

#Remove unnecessary columns
ORALspecies2 <- ORAL2[,c(-2,-4:-17)]

mom_figure1B <- ORALspecies2 %>%
  filter(grepl("pregnancy", Measurement, ignore.case = TRUE) & 
           grepl("tongue swab", sample_type, ignore.case = TRUE))

momt_figure1B <- ORALspecies2 %>%
  filter(grepl("pregnancy", Measurement, ignore.case = TRUE) & 
           grepl("tooth plaque", sample_type, ignore.case = TRUE))

tongue1_figure1B <- ORALspecies2 %>%
  filter(grepl("1 month", Measurement, ignore.case = TRUE) & 
           grepl("tongue swab", sample_type, ignore.case = TRUE))

tongue6_figure1B <- ORALspecies2 %>%
  filter(grepl("6 months", Measurement, ignore.case = TRUE) & 
           grepl("tongue swab", sample_type, ignore.case = TRUE))

tooth6_figure1B <- ORALspecies2 %>%
  filter(grepl("6 months", Measurement, ignore.case = TRUE) & 
           grepl("tooth plaque", sample_type, ignore.case = TRUE))


tonguemwide_figure1B <- pivot_wider(data = mom_figure1B, 
                                    id_cols = Pair,  # Specify the columns to keep as identifiers
                                    names_from = species,  # Specify the column containing the new column names
                                    values_from = Taxonomic_abundance) 
tonguemwide_figure1B <- as.matrix(tonguemwide_figure1B)
tonguemwide_figure1B[is.na(tonguemwide_figure1B)] <- 0
rownames(tonguemwide_figure1B) <- tonguemwide_figure1B[,1]
tonguemwide_figure1B <- as.data.frame(tonguemwide_figure1B)
tonguemwide_figure1B$Pair <- NULL

tonguemwide_figure1B[,1:361] %<>% lapply(function(x) as.numeric(as.character(x)))
tonguemwide_figure1B <- tonguemwide_figure1B[,1:361] %>% dplyr::select(which(!colSums(tonguemwide_figure1B[,1:361], na.rm=TRUE) %in% 0))
tonguemwide_figure1B <- tonguemwide_figure1B/100


tonguemtwide_figure1B <- pivot_wider(data = momt_figure1B, 
                                     id_cols = Pair,  # Specify the columns to keep as identifiers
                                     names_from = species,  # Specify the column containing the new column names
                                     values_from = Taxonomic_abundance) 
tonguemtwide_figure1B <- as.matrix(tonguemtwide_figure1B)
tonguemtwide_figure1B[is.na(tonguemtwide_figure1B)] <- 0
rownames(tonguemtwide_figure1B) <- tonguemtwide_figure1B[,1]
tonguemtwide_figure1B <- as.data.frame(tonguemtwide_figure1B)
tonguemtwide_figure1B$Pair <- NULL

tonguemtwide_figure1B[,1:431] %<>% lapply(function(x) as.numeric(as.character(x)))
tonguemtwide_figure1B <- tonguemtwide_figure1B[,1:431] %>% dplyr::select(which(!colSums(tonguemtwide_figure1B[,1:431], na.rm=TRUE) %in% 0))
tonguemtwide_figure1B <- tonguemtwide_figure1B/100

tongue1wide_figure1B <- pivot_wider(data = tongue1_figure1B, 
                                    id_cols = Pair,  # Specify the columns to keep as identifiers
                                    names_from = species,  # Specify the column containing the new column names
                                    values_from = Taxonomic_abundance) 
tongue1wide_figure1B <- as.matrix(tongue1wide_figure1B)
tongue1wide_figure1B[is.na(tongue1wide_figure1B)] <- 0
rownames(tongue1wide_figure1B) <- tongue1wide_figure1B[,1]
tongue1wide_figure1B <- as.data.frame(tongue1wide_figure1B)
tongue1wide_figure1B$Pair <- NULL

library(dplyr)
#tongue1wide[,1:56] %<>% lapply(function(x) as.numeric(as.character(x)))
#tongue1wide <- tongue1wide[,1:56] %>% dplyr::select(which(!colSums(tongue1wide[,1:56], na.rm=TRUE) %in% 0))
tongue1wide_figure1B[,1:90] %<>% lapply(function(x) as.numeric(as.character(x)))
tongue1wide_figure1B <- tongue1wide_figure1B[,1:90] %>% dplyr::select(which(!colSums(tongue1wide_figure1B[,1:90], na.rm=TRUE) %in% 0))
tongue1wide_figure1B <- tongue1wide_figure1B/100

tongue6wide_figure1B <- pivot_wider(data = tongue6_figure1B, 
                                    id_cols = Pair,  # Specify the columns to keep as identifiers
                                    names_from = species,  # Specify the column containing the new column names
                                    values_from = Taxonomic_abundance) 
tongue6wide_figure1B <- as.matrix(tongue6wide_figure1B)
tongue6wide_figure1B[is.na(tongue6wide_figure1B)] <- 0
rownames(tongue6wide_figure1B) <- tongue6wide_figure1B[,1]
tongue6wide_figure1B <- as.data.frame(tongue6wide_figure1B)
tongue6wide_figure1B$Pair <- NULL

library(dplyr)
tongue6wide_figure1B[,1:216] %<>% lapply(function(x) as.numeric(as.character(x)))
tongue6wide_figure1B <- tongue6wide_figure1B[,1:216] %>% dplyr::select(which(!colSums(tongue6wide_figure1B[,1:216], na.rm=TRUE) %in% 0))
tongue6wide_figure1B <- tongue6wide_figure1B/100

tooth6wide_figure1B <- pivot_wider(data = tooth6_figure1B, 
                                   id_cols = Pair,  # Specify the columns to keep as identifiers
                                   names_from = species,  # Specify the column containing the new column names
                                   values_from = Taxonomic_abundance) 
tooth6wide_figure1B <- as.matrix(tooth6wide_figure1B)
tooth6wide_figure1B[is.na(tooth6wide_figure1B)] <- 0
rownames(tooth6wide_figure1B) <- tooth6wide_figure1B[,1]
tooth6wide_figure1B <- as.data.frame(tooth6wide_figure1B)
tooth6wide_figure1B$Pair <- NULL

library(dplyr)
tooth6wide_figure1B[,1:210] %<>% lapply(function(x) as.numeric(as.character(x)))
tooth6wide_figure1B <- tooth6wide_figure1B[,1:210] %>% dplyr::select(which(!colSums(tooth6wide_figure1B[,1:210], na.rm=TRUE) %in% 0))
tooth6wide_figure1B <- tooth6wide_figure1B/100


tonguemwide_figure1B$Pair <- rownames(tonguemwide_figure1B)
tonguemtwide_figure1B$Pair <- rownames(tonguemtwide_figure1B)
tongue1wide_figure1B$Pair <- rownames(tongue1wide_figure1B)
tongue6wide_figure1B$Pair <- rownames(tongue6wide_figure1B)
tooth6wide_figure1B$Pair <- rownames(tooth6wide_figure1B)

tonguemlong_figure1B <- tonguemwide_figure1B%>%
  pivot_longer(
    cols = c(,1:361), # columns to pivot, excluding the identifier column(s)
    names_to = "species", # name of the new column that will contain the names of the original columns
    values_to = "Taxonomic_abundance" # name of the new column that will contain the values of the original columns
  )

tonguemtlong_figure1B <- tonguemtwide_figure1B%>%
  pivot_longer(
    cols = c(,1:431), # columns to pivot, excluding the identifier column(s)
    names_to = "species", # name of the new column that will contain the names of the original columns
    values_to = "Taxonomic_abundance" # name of the new column that will contain the values of the original columns
  )

tongue1long_figure1B <- tongue1wide_figure1B%>%
  pivot_longer(
    cols = c(,1:90), # columns to pivot, excluding the identifier column(s)
    names_to = "species", # name of the new column that will contain the names of the original columns
    values_to = "Taxonomic_abundance" # name of the new column that will contain the values of the original columns
  )

tongue6long_figure1B <- tongue6wide_figure1B%>%
  pivot_longer(
    cols = c(,1:216), # columns to pivot, excluding the identifier column(s)
    names_to = "species", # name of the new column that will contain the names of the original columns
    values_to = "Taxonomic_abundance" # name of the new column that will contain the values of the original columns
  )

tooth6long_figure1B <- tooth6wide_figure1B%>%
  pivot_longer(
    cols = c(,1:210), # columns to pivot, excluding the identifier column(s)
    names_to = "species", # name of the new column that will contain the names of the original columns
    values_to = "Taxonomic_abundance" # name of the new column that will contain the values of the original columns
  )



averageM_figure1B <- tonguemlong_figure1B %>%
  group_by(species) %>%
  dplyr::summarize(mean_abundance = mean(Taxonomic_abundance),
                   standard_error = sd(Taxonomic_abundance, na.rm = TRUE) / sqrt(n()))
averageMT_figure1B <- tonguemtlong_figure1B %>%
  group_by(species) %>%
  dplyr::summarize(mean_abundance = mean(Taxonomic_abundance),
                   standard_error = sd(Taxonomic_abundance, na.rm = TRUE) / sqrt(n()))
averageT1_figure1B <- tongue1long_figure1B %>%
  group_by(species) %>%
  dplyr::summarize(mean_abundance = mean(Taxonomic_abundance),
                   standard_error = sd(Taxonomic_abundance, na.rm = TRUE) / sqrt(n()))
averageT6_figure1B <- tongue6long_figure1B %>%
  group_by(species) %>%
  dplyr::summarize(mean_abundance = mean(Taxonomic_abundance),
                   standard_error = sd(Taxonomic_abundance, na.rm = TRUE) / sqrt(n()))
averageTO6_figure1B <- tooth6long_figure1B %>%
  group_by(species) %>%
  dplyr::summarize(mean_abundance = mean(Taxonomic_abundance),
                   standard_error = sd(Taxonomic_abundance, na.rm = TRUE) / sqrt(n()))

averageM_figure1B$Measurement_point <- as.factor("Mother")
averageMT_figure1B$Measurement_point <- as.factor("Mother")
averageT1_figure1B$Measurement_point <- as.factor("1 month")
averageT6_figure1B$Measurement_point <- as.factor("6 months")
averageTO6_figure1B$Measurement_point <- as.factor("6 months")
averageM_figure1B$sample_type <- as.factor("tongue dorsum")
averageMT_figure1B$sample_type <- as.factor("dental plaque")
averageT1_figure1B$sample_type <- as.factor("tongue dorsum")
averageT6_figure1B$sample_type <- as.factor("tongue dorsum")
averageTO6_figure1B$sample_type <- as.factor("dental plaque")

######### Figure 5A (bargraph, further customized in Adobe Illustrator)
filteraverageT6 <- subset(averageT6_figure1B, (grepl('Pauljensenia sp900541895|Rothia sp902373285|Streptococcus sp000187445', species)))

ggplot(filteraverageT6,aes(reorder(species,-mean_abundance),mean_abundance))+
  geom_bar(position="dodge",stat = "identity",col='black',size=0.5) +
  geom_errorbar(position=position_dodge2(),aes(ymin=mean_abundance-standard_error, ymax=mean_abundance+standard_error), size=0.4,width=0.6)+
  theme_classic()+
  ylab("Relative abundance (%)")

######################################################################
df_top10tongue1_figure1B <- averageT1_figure1B %>%
  slice_max(mean_abundance, n = 10) %>%       # Retain the top 5 highest values in column D
  ungroup()

df_top10tongue6_figure1B <- averageT6_figure1B %>%
  slice_max(mean_abundance, n = 10) %>%       # Retain the top 5 highest values in column D
  ungroup()

df_top10mom_figure1B <- averageM_figure1B %>%
  slice_max(mean_abundance, n = 10) %>%       # Retain the top 5 highest values in column D
  ungroup()

df_top10momTO_figure1B <- averageMT_figure1B %>%
  slice_max(mean_abundance, n = 10) %>%       # Retain the top 5 highest values in column D
  ungroup()

df_top10tongueTO6_figure1B <- averageTO6_figure1B %>%
  slice_max(mean_abundance, n = 10) %>%       # Retain the top 5 highest values in column D
  ungroup()



species_abundance <- rbind(df_top10mom_figure1B, df_top10tongue1_figure1B)
species_abundance2 <- rbind(species_abundance, df_top10tongue6_figure1B)
species_abundance3 <- rbind(species_abundance2, df_top10tongueTO6_figure1B)
species_abundance4 <- rbind(species_abundance3, df_top10momTO_figure1B)

add_genus <- ORALspecies2[,c(6,7)]
add_genus <- distinct(add_genus)
species_abundance5 <- merge(species_abundance4,add_genus, by='species',all.x=T)

genus_colors <- c("g__Streptococcus" = "#FFC107", 
                  "g__Rothia" ='#D81B60',
                  "g__Veillonella"='#1E88E5',
                  "g__Prevotella"='#004D40',
                  "g__Pauljensenia"='#FEC181',
                  "g__Neisseria" = "#2C8D40", 
                  "g__Lactobacillus" = "#7A55D7", 
                  "g__Actinomyces" = "#5E9FD0",
                  "g__Scardovia"='coral',
                  "g__Haemophilus_D"="khaki",
                  "g__Alloprevotella"="burlywood"
                  
)

#Figure 1B
ggplot(species_abundance5,aes(reorder(species,mean_abundance),mean_abundance, fill=genus))+
  geom_bar(position="dodge",stat = "identity",col='black',size=0.5) +
  geom_errorbar(position=position_dodge2(),aes(ymin=mean_abundance-standard_error, ymax=mean_abundance+standard_error), size=0.5,width=0.9)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 0 , hjust=0.5, size=6),
        axis.text.y = element_text(size=8))+
  facet_wrap(~Measurement_point*sample_type,scales='fixed', ncol=5)+
  coord_flip()+
  theme(legend.position="right")+
  guides(fill=guide_legend(nrow=20,byrow=TRUE))+
  scale_fill_manual(values = genus_colors)

#Make FIGURE 1C
ORAL$species <- sapply(ORAL$clade_name, function(x) {
  # Check if the clade_name contains 's__'
  if (grepl("s__", x)) {
    # Extract substring starting from 's__'
    sub(".*(s__[^|]*).*", "\\1", x)
  } else if (grepl("undescribed", x, ignore.case = TRUE)) {
    # If 'undescribed' is found, retain the entire string
    x
  } else {
    # Otherwise, return NA or empty string
    NA
  }
})

#Remove unnecessary columns
ORALspecies <- ORAL[,c(-2,-4:-17,-21)]


#Re-order dataframe to matrix to be used as input for PCoA
df_wide_figure1b <- ORALspecies %>%
  pivot_wider(names_from = species, values_from = Taxonomic_abundance, values_fill = 0)

df_wide_figure1b <- as.matrix(df_wide_figure1b)
rownames(df_wide_figure1b) <- df_wide_figure1b[,c(1)]
df_wide_figure1b <- as.data.frame(df_wide_figure1b)

df_wide_figure1b[,5:903] %<>% lapply(function(x) as.numeric(as.character(x)))
df_wide_figure1b_2 <- df_wide_figure1b[,5:903] %>% dplyr::select(which(!colSums(df_wide_figure1b[,5:903], na.rm=TRUE) %in% 0))


#Exclude non-numeric columns like Sample_file, timepoint etc. from df_wide_figure1b_2
abundance_matrix <- df_wide_figure1b_2[ , !(names(df_wide_figure1b_2) %in% c("Sample_file", "Measurement", "Pair", "sample_type"))]
abundance_matrix2 <- as.data.frame(t(abundance_matrix))

#Compute the Bray-Curtis dissimilarity matrix and PCoA
distance_matrix <- vegdist(abundance_matrix, method = "bray")
pcoa_output <- cmdscale(distance_matrix, eig = TRUE, k = 2)
pcoa_scores <- as.data.frame(scores(pcoa_output))

#Gather PCoA output
pcoa_data <- data.frame(Sample_file = df_wide_figure1b$Sample_file,
                      sample_type = df_wide_figure1b$sample_type,
                      Measurement_point = df_wide_figure1b$Measurement,
                      check = rownames(pcoa_scores),
                      PCoA1 = pcoa_scores$Dim1, 
                      PCoA2 = pcoa_scores$Dim2)

#Assign custom colors to timepoint
custom_colors <- c("pregnancy" = "orchid2", 
                   "1 month" ='#B91C1C',
                   "6 months"='#0D9488'
                   
)

#Calculate percentage of variance explained by each PC
pcoa_output$eig[1]/sum(pcoa_output$eig)
pcoa_output$eig[2]/sum(pcoa_output$eig)

ggplot(pcoa_data, aes(x = PCoA1, y = PCoA2, shape=sample_type)) +
  geom_point(size=2.2)+
  geom_point(aes(x = PCoA1, y = PCoA2, color = Measurement_point, shape = sample_type), size=2) +
  labs(x = "PCoA1 (20%)", y = "PCoA2 (16%)")+
  theme_classic()+
  scale_color_manual(values = custom_colors)





#Make Figure 1C - Species Richness bargraphs
##RICHNESS
mom$frequency <- as.numeric(c(1))
momt$frequency <- as.numeric(c(1))
tongue1$frequency <- as.numeric(c(1))
tongue6$frequency <- as.numeric(c(1))
tooth6$frequency <- as.numeric(c(1))

sum_speciesM <- mom %>%
  group_by(Pair) %>%
  dplyr::summarize(sum_frequency = sum(frequency))
sum_speciesMTO <- momt %>%
  group_by(Pair) %>%
  dplyr::summarize(sum_frequency = sum(frequency))
sum_species1 <- tongue1 %>%
  group_by(Pair) %>%
  dplyr::summarize(sum_frequency = sum(frequency))
sum_species6 <- tongue6 %>%
  group_by(Pair) %>%
  dplyr::summarize(sum_frequency = sum(frequency))
sum_speciesTO6 <- tooth6 %>%
  group_by(Pair) %>%
  dplyr::summarize(sum_frequency = sum(frequency))


average_speciesM <- sum_speciesM %>%
  dplyr::summarize(mean_frequency = mean(sum_frequency),
                   standard_error = sd(sum_frequency, na.rm = TRUE) / sqrt(n()))
average_speciesMTO <- sum_speciesMTO %>%
  dplyr::summarize(mean_frequency = mean(sum_frequency),
                   standard_error = sd(sum_frequency, na.rm = TRUE) / sqrt(n()))
average_species1 <- sum_species1 %>%
  dplyr::summarize(mean_frequency = mean(sum_frequency),
                   standard_error = sd(sum_frequency, na.rm = TRUE) / sqrt(n()))
average_species6 <- sum_species6 %>%
  dplyr::summarize(mean_frequency = mean(sum_frequency),
                   standard_error = sd(sum_frequency, na.rm = TRUE) / sqrt(n()))
average_speciesTO6 <- sum_speciesTO6 %>%
  dplyr::summarize(mean_frequency = mean(sum_frequency),
                   standard_error = sd(sum_frequency, na.rm = TRUE) / sqrt(n()))


average_speciesM$Measurement <- as.character(c("Mother"))
average_speciesMTO$Measurement <- as.character(c("Mother"))
average_species1$Measurement <- as.character(c("Infant (1 month)"))
average_species6$Measurement <- as.character(c("Infant (6 months)"))
average_speciesTO6$Measurement <- as.character(c("Infant (6 months)"))


# Make species richness (Figure 1B) bargraphs (further customized in Adobe Illustrator)
ggplot(average_speciesM,aes(Measurement,mean_frequency))+
  geom_bar(position="dodge",stat = "identity",col='black',size=0.5) +
  geom_errorbar(position=position_dodge2(),aes(ymin=mean_frequency-standard_error, ymax=mean_frequency+standard_error), size=0.4,width=0.6)+
  theme_classic()+
  ylab("Alpha diversity (species richness)")+
  theme(axis.text.y=element_blank())+
  coord_flip()+
  ylim(0,180)

ggplot(average_speciesMTO,aes(Measurement,mean_frequency))+
  geom_bar(position="dodge",stat = "identity",col='black',size=0.5) +
  geom_errorbar(position=position_dodge2(),aes(ymin=mean_frequency-standard_error, ymax=mean_frequency+standard_error), size=0.4,width=0.6)+
  theme_classic()+
  ylab("Alpha diversity (species richness)")+
  theme(axis.text.y=element_blank())+
  coord_flip()+
  ylim(0,180)


ggplot(average_species1,aes(Measurement,mean_frequency))+
  geom_bar(position="dodge",stat = "identity",col='black',size=0.5) +
  geom_errorbar(position=position_dodge2(),aes(ymin=mean_frequency-standard_error, ymax=mean_frequency+standard_error), size=0.4,width=0.6)+
  theme_classic()+
  ylab("Alpha diversity (species richness)")+
  theme(axis.text.y=element_blank())+
  coord_flip()+
  ylim(0,180)

ggplot(average_species6,aes(Measurement,mean_frequency))+
  geom_bar(position="dodge",stat = "identity",col='black',size=0.5) +
  geom_errorbar(position=position_dodge2(),aes(ymin=mean_frequency-standard_error, ymax=mean_frequency+standard_error), size=0.4,width=0.6)+
  theme_classic()+
  ylab("Alpha diversity (species richness)")+
  theme(axis.text.y=element_blank())+
  coord_flip()+
  ylim(0,180)

ggplot(average_speciesTO6,aes(Measurement,mean_frequency))+
  geom_bar(position="dodge",stat = "identity",col='black',size=0.5) +
  geom_errorbar(position=position_dodge2(),aes(ymin=mean_frequency-standard_error, ymax=mean_frequency+standard_error), size=0.4,width=0.6)+
  theme_classic()+
  ylab("Alpha diversity (species richness)")+
  theme(axis.text.y=element_blank())+
  coord_flip()+
  ylim(0,180)








######## supplementary PCoA's
tongue <- ORALspecies %>%
  filter(grepl("tongue swab", sample_type, ignore.case = TRUE))

#Re-order dataframe to matrix to be used as input for PCoA
df_wide_figure1b <- tongue %>%
  pivot_wider(names_from = species, values_from = Taxonomic_abundance, values_fill = 0)

df_wide_figure1b <- as.matrix(df_wide_figure1b)
rownames(df_wide_figure1b) <- df_wide_figure1b[,c(1)]
df_wide_figure1b <- as.data.frame(df_wide_figure1b)

df_wide_figure1b[,5:695] %<>% lapply(function(x) as.numeric(as.character(x)))
df_wide_figure1b_2 <- df_wide_figure1b[,5:695] %>% dplyr::select(which(!colSums(df_wide_figure1b[,5:695], na.rm=TRUE) %in% 0))


#Exclude non-numeric columns like Sample_file, timepoint etc. from df_wide_figure1b_2
abundance_matrix <- df_wide_figure1b_2[ , !(names(df_wide_figure1b_2) %in% c("Sample_file", "Measurement", "Pair", "sample_type"))]
abundance_matrix2 <- as.data.frame(t(abundance_matrix))

#Compute the Bray-Curtis dissimilarity matrix and PCoA
distance_matrix <- vegdist(abundance_matrix, method = "bray")
pcoa_output <- cmdscale(distance_matrix, eig = TRUE, k = 2)
pcoa_scores <- as.data.frame(scores(pcoa_output))

#Gather PCoA output
pcoa_data <- data.frame(Sample_file = df_wide_figure1b$Sample_file,
                        sample_type = df_wide_figure1b$sample_type,
                        Measurement_point = df_wide_figure1b$Measurement,
                        check = rownames(pcoa_scores),
                        PCoA1 = pcoa_scores$Dim1, 
                        PCoA2 = pcoa_scores$Dim2)

#Assign custom colors to timepoint
custom_colors <- c("pregnancy" = "#B17A6B", 
                   "1 month" ='#7A9B3D',
                   "6 months"='#B92428'
                   
)

#Calculate percentage of variance explained by each PC
pcoa_output$eig[1]/sum(pcoa_output$eig)
pcoa_output$eig[2]/sum(pcoa_output$eig)

ggplot(pcoa_data, aes(x = PCoA1, y = PCoA2, shape=sample_type)) +
  geom_point(size=2.2)+
  geom_point(aes(x = PCoA1, y = PCoA2, color = Measurement_point, shape = sample_type), size=2) +
  labs(x = "PCoA1 (26%)", y = "PCoA2 (10%)")+
  theme_classic()+
  scale_color_manual(values = custom_colors)



## INFANT NUTRITION
######## supplementary PCoA's
infant <- ORALspecies %>%
  filter(!grepl("pregnancy|1 month", Measurement, ignore.case = TRUE))

nutri <- read.csv("PROJECTS/CLEAN_oral_microbiome/ISME_JOURNAL/ISME_REVIEWS/PLOSbiology/REVIEW_PLOTCOMPBIO/GITHUBoral/Figure 1/metadata_nutrition.tsv", sep = '\t', header=TRUE, stringsAsFactors = FALSE)
nutri <- nutri %>%
  filter(!grepl("pregnancy|1 month", Measurement, ignore.case = TRUE))

nutri <- nutri[,c(1,5,7)]

infant2 <- merge(infant,nutri,by='Sample_file',all.x=T)


#Re-order dataframe to matrix to be used as input for PCoA
df_wide_figure1b <- infant2 %>%
  pivot_wider(names_from = species, values_from = Taxonomic_abundance, values_fill = 0)

df_wide_figure1b <- as.matrix(df_wide_figure1b)
rownames(df_wide_figure1b) <- df_wide_figure1b[,c(1)]
df_wide_figure1b <- as.data.frame(df_wide_figure1b)

df_wide_figure1b[,7:359] %<>% lapply(function(x) as.numeric(as.character(x)))
df_wide_figure1b_2 <- df_wide_figure1b[,7:359] %>% dplyr::select(which(!colSums(df_wide_figure1b[,7:359], na.rm=TRUE) %in% 0))


#Exclude non-numeric columns like Sample_file, timepoint etc. from df_wide_figure1b_2
abundance_matrix <- df_wide_figure1b_2[ , !(names(df_wide_figure1b_2) %in% c("Sample_file", "Measurement", "Pair", "sample_type"))]
abundance_matrix2 <- as.data.frame(t(abundance_matrix))

#Compute the Bray-Curtis dissimilarity matrix and PCoA
distance_matrix <- vegdist(abundance_matrix, method = "bray")
pcoa_output <- cmdscale(distance_matrix, eig = TRUE, k = 2)
pcoa_scores <- as.data.frame(scores(pcoa_output))

#Gather PCoA output
pcoa_data <- data.frame(Sample_file = df_wide_figure1b$Sample_file,
                        sample_type = df_wide_figure1b$sample_type,
                        Measurement_point = df_wide_figure1b$Measurement.x,
                        Milk_feeding = df_wide_figure1b$Milk_feeding_6mnt,
                        check = rownames(pcoa_scores),
                        PCoA1 = pcoa_scores$Dim1, 
                        PCoA2 = pcoa_scores$Dim2)

#Assign custom colors to timepoint
custom_colors <- c("Breastmilk" = "orchid2", 
                   "Mixed" ='#B91C1C',
                   "Formula"='#0D9488'
                   
)

#Calculate percentage of variance explained by each PC
pcoa_output$eig[1]/sum(pcoa_output$eig)
pcoa_output$eig[2]/sum(pcoa_output$eig)

ggplot(pcoa_data, aes(x = PCoA1, y = PCoA2, shape=sample_type)) +
  geom_point(size=2.2)+
  geom_point(aes(x = PCoA1, y = PCoA2, color = Milk_feeding, shape = sample_type)) +
  labs(x = "PCoA1 (23%)", y = "PCoA2 (13%)")+
  theme_classic()+
  scale_color_manual(values = custom_colors)


adonis_result <- adonis2(distance_matrix ~ Milk_feeding_6mnt, data = nutri)



#CHECK SHARED SPECIES AMONG TIMEPOINTS
# Get unique species per time point
species_by_time <- ORALspecies2 %>%
  dplyr::select(Measurement, species) %>%
  distinct() %>%
  dplyr::group_by(Measurement) %>%
  dplyr::summarise(species_list = list(unique(species)))

# Extract species sets
pregnancy_species <- species_by_time$species_list[[which(species_by_time$Measurement == "pregnancy")]]
month1_species     <- species_by_time$species_list[[which(species_by_time$Measurement == "1 month")]]
month6_species     <- species_by_time$species_list[[which(species_by_time$Measurement == "6 months")]]

# Intersections
shared_preg_1mo   <- intersect(pregnancy_species, month1_species)
shared_preg_6mo   <- intersect(pregnancy_species, month6_species)
shared_1mo_6mo    <- intersect(month1_species, month6_species)

# Counts
length(shared_preg_1mo)   # species shared between pregnancy and 1 month
length(shared_preg_6mo)   # species shared between pregnancy and 6 months
length(shared_1mo_6mo)    # species shared between 1 month and 6 months

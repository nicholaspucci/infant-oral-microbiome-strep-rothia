# **Infant Oral Microbiome Development** 
This repository contains the code and data used to generate the main figures in 
**"Functional adaptations and metabolic interactions of undescribed *Streptococcus* and *Rothia* species during infant oral microbiome development"**.

## **Overview** 
This study analyzes oral microbiomes from a longitudinal cohort of 24 mother-infant dyads at 1 and 6 months postpartum using shotgun metagenomics. We identified two previously undescribed *Streptococcus* and *Rothia* species that are prevalent, abundant, and strongly co-occurring members of the oral microbiome of six-month-old infants. The repository provides code to reproduce the main figures presented in the paper.

## **Repository Structure**

/

├── README.md               \# This file

├── scripts/                \# Analysis and visualization scripts  
│   ├── figure1/            \# Code for Figure 1  
│   │   └── pucci2025\_oral\_script\_figure1.R  
│   ├── figure2/            \# Code for Figure 2  
│   │   └── pucci2025\_oral\_script\_figure2.sh  
│   ├── figure3-4/            \# Code for Figure 3  
│   │   └── pucci2025\_oral\_script\_figure3\_figure4.R  
│   └── figure5/            \# Code for Figure 5  
│   │   └── pucci2025\_oral\_script\_figure5.R

├── data/                   \# Data files used in the analysis  
│   ├── metadata\_Pucci2025.tsv  \# Study metadata   
│   ├── aims-oral-sylph.tsv     \# Taxonomic abundance data  
│   ├── aims-oral-r220.tsv      \# Taxa information  
│   └── gtdb\_genome\_clade\_id.csv \# GTDB taxonomy information  
│   ├── CAZymes-rothia.tsv      \# CAZyme copy-number Rothia  
│   └── rothia-isolates-match-cazymes.tsv      \#MAG information Rothia  
│   ├── CAZymes-streptococcus.tsv      \# CAZyme copy-number Streptococcus  
│   └── strepto-isolates-match-cazymes.tsv      \#MAG information Streptococcus  
│   ├── concatenated\_AAseqs-rothia.fa      \# Protein sequences Rothia tree  
│   ├── concatenated\_AAseqs-streptococcus.fa      \# Protein sequences Streptococcus tree  
│   └── genus\_phylo\_analysis.txt      \# genus information for pucci2025\_oral\_script\_figure2.sh  
│   ├── KOlist-NOSreduction.txt      \# List KOs in NO3-NO2-NO pathway  
│   ├── lectins\_strepto.tsv      \#Adhesin/lectin copy-numbers Streptococcus  
│   ├── network4\_PTM\_attributes.tsv      \# Metabolite information  
│   ├── network4\_PTM\_network.tsv      \# Genome-metabolite Node-edge list  
│   ├── PhyloMInt\_result.tsv      \# PhyloMInt metabolic complementarity indices  
│   ├── ROTHIAgenomic-metabolism2\_matrix75-module\_pathwise\_presence-MATRIX.txt      \# Module abundance matrix Rothia  
│   ├── ROTHIAgenomic75-metabolism2\_modules.txt      \# Module information Rothia  
│   ├── STREPTOgenomic-metabolism2\_matrix75-module\_pathwise\_presence-MATRIX.txt      \# Module abundance matrix Streptococcus  
│   ├── STREPTOgenomic75-metabolism2\_modules.txt      \# Module information Streptococcus

└── LICENSE                 \# License information

## **Data Description**

### **Main Data Files**

* **metadata\_Pucci2025.tsv**: Contains sample metadata including subject identifiers, timepoints (pregnancy, 1 month, 6 months), and sample types (tongue swab, tooth plaque)  
* **aims-oral-sylph.tsv**: Contains taxonomic abundance data for each sample  
* **aims-oral-r220.tsv**: Contains taxonomic information  
* **gtdb\_genome\_clade\_id.csv**: Contains GTDB taxonomy classifications for genomes

### **Phylogenomic Analysis Files**

* **genus\_phylo\_analysis.txt**: List of genera (Streptococcus and Rothia) for phylogenomic analysis  
* **concatenated\_AAseqs-Streptococcus.fa**: Concatenated amino acid sequences of single-copy core genes from Streptococcus genomes  
* **concatenated\_AAseqs-Rothia.fa**: Concatenated amino acid sequences of single-copy core genes from Rothia genomes

**Functional and Genomic Data for Figures 3–5**

* **CAZymes-streptococcus.tsv:** Copy-number matrix of carbohydrate-active enzymes (CAZymes) in *Streptococcus* genomes.  
* **Strepto-isolates-match-cazymes.tsv:** Mapping between *Streptococcus* genome names and isolate/species identifiers.  
* **CAZymes-rothia.tsv:** Copy-number matrix of CAZymes in *Rothia* genomes.  
* **Rothia-isolates-match-cazymes.tsv:** Mapping between *Rothia* genome names and isolate/species identifiers.  
* **STREPTOgenomic75-metabolism2\_modules.txt:** KEGG module abundance profiles for *Streptococcus* genomes.  
* **ROTHIAgenomic75-metabolism2\_modules.txt:** KEGG module abundance profiles for *Rothia* genomes.  
* **STREPTOgenomic-metabolism2\_matrix75-module\_pathwise\_presence-MATRIX.txt:** Binary matrix (presence/absence) of KEGG modules in *Streptococcus* genomes.  
* **ROTHIAgenomic-metabolism2\_matrix75-module\_pathwise\_presence-MATRIX.txt:** Binary matrix (presence/absence) of KEGG modules in *Rothia* genomes.  
* **Lectins\_strepto.tsv:** Copy-numbers of predicted adhesin and lectin genes in *Streptococcus* genomes.  
* **KOlist-NOSreduction.txt:** List of KEGG Orthologs involved in nitrate/nitrite reduction pathways (used to assess enterosalivary nitrate metabolism).  
* **ROTHIAgenomic-metabolism2\_matrix75-enzyme\_hits-MATRIX.txt:** Copy-number matrix of nitrate/nitrite pathway enzyme genes in *Rothia* genomes.

### **Data for PhyloMInt-Based Interaction Modeling (Figure 5\)**

* **PhyloMInt\_result.tsv:** Output from PhyloMInt software: pairwise metabolic complementarity scores between genome pairs.  
* **network4\_PTM\_network.tsv:** Node-edge list for species–metabolite interaction predictions.  
* **network4\_PTM\_attributes.tsv:** Metadata including superclasses for each metabolite in the interaction network.

## **Figures Overview**

### **Figure 1** Shows the taxonomic composition of the oral microbiome at different timepoints:

* **1A**: Genus-level relative abundance in maternal and infant samples  
* **1B**: Species-level composition of the most abundant genera  
* **1C**: Beta diversity analysis using PCoA of Bray-Curtis dissimilarities  
* **1D**: Alpha diversity (species richness) across sample types and timepoints

### **Figure 2** Phylogenomic analysis of *Streptococcus* and *Rothia* species:

* Generates phylogenomic trees based on concatenated amino acid sequences of single-copy core genes  
* Includes 33 *Streptococcus* and 63 *Rothia* genomes/MAGs  
* Trees were visualized and customized using iTOL (Interactive Tree Of Life)

### **Figure 3–4** Genomic and functional characterization of novel *Streptococcus* and *Rothia* species:

* **Metabolic modules (Fig. 3 and 4\)** : Presence/absence and copy number heatmaps derived from KEGG module annotations and enrichment analyses.  
* **CAZymes**: Significant CAZyme families (e.g., GH, GT classes) visualized in heatmaps.  
* **Lectins/Adhesins** (***Streptococcus***, **Fig. 3**): Gene copy-numbers of predicted lectins/adhesins.  
* **Enterosalivary nitrate reduction (*Rothia*, Fig. 4\)**: Gene content of the NO₃⁻–NO₂⁻ reduction pathway assessed in *Rothia* species.

Outputs were further customized using iTOL and Adobe Illustrator.

### **Figure 5** Metabolic interactions between co-occurring *Streptococcus* and *Rothia* species:

* **5B**: Metabolic complementarity scores computed with PhyloMInt.  
* **5C**: Species–metabolite heatmap based on predicted secretion/uptake traits.

Final visualizations were refined in Adobe Illustrator.

## **Requirements**

### **R Packages** The R analysis scripts require the following packages:

* tidyr  
* magrittr  
* dplyr  
* corrplot  
* vegan  
* ggplot2  
* readxl  
* pheatmap  
* stringr  

### **External Software**

* anvi'o (v8.0) \- for phylogenomic analysis  
* iTOL (Interactive Tree Of Life) \- for tree visualization  
* Adobe Illustrator (v28.7.5) \- for figure customization

## **Usage**

### **To reproduce Figure 1 from the paper:**
#Update the working directory path in the script
#setwd("/your/path/here") 
#source("scripts/figure1/pucci2025\_oral\_script\_figure1.R")

### **To reproduce Figure 2 (phylogenomic trees):**
#Make sure you have anvi'o v8.0 installed and activated
#Prepare a text file 'genus\_phylo\_analysis.txt' with genera names (Streptococcus and Rothia)
#Ensure concatenated amino acid sequence files exist in the working directory

./scripts/figure2/pucci2025\_oral\_script\_figure2.sh

### **To reproduce Figure 3 and 4 from the paper:**
#Update the working directory path in the script
#setwd("/your/path/here") 
#source("scripts/figure3-4/pucci2025\_oral\_script\_figure3\_figure4.R")

### **To reproduce Figure 5 from the paper:**
#Update the working directory path in the script
#setwd("/your/path/here") 
#source("scripts/figure5/pucci2025\_oral\_script\_figure5.R")

Note: The scripts produce raw plots/trees that were further customized in Adobe Illustrator and iTOL for the final manuscript figures.

## **Data Access** 
The raw metagenomic sequencing data has been deposited in the European Nucleotide Archive (ENA) under accession number PRJEB88622). The dataset will be publicly available on 2025-09-01

## **Citation** 
If you use this code or data in your research, please cite:

Pucci, N., et al. (2025). Functional adaptations and metabolic interactions of undescribed Streptococcus and Rothia species during infant oral microbiome development. \[Journal Name\]. doi: \[doi\]

## **License** 
This project is licensed under the \[LICENSE\] \- see the LICENSE file for details.

## **Contact** 
For questions about the code or data, please contact: \[[n.pucci@amsterdamumc.nl](mailto:n.pucci@amsterdamumc.nl)\] 



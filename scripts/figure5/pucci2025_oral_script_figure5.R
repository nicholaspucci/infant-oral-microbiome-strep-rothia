setwd("/Users/path/to/directory")

# Figure 5
library(pheatmap)
library(stringr)

################# Figure 1 (bargraph, see script:'pucci2025_oral_script_figure1.R)

################# Figure 5B
mint <- read.csv("PhyloMInt_result.tsv", sep = '\t', header=TRUE, stringsAsFactors = FALSE)
mint2 <- mint[,-3]

mint2$comparison <- paste(mint2$A, "-", mint2$B)
mint2$reference_genome <- sub(" -.*", "", mint2$comparison)
mint2 <- mint2 %>%
  mutate(
    reference_genome = str_extract(comparison, "^[^-]+"),
    tested_genome = str_extract(comparison, "(?<=-).*")
  )


ggplot(mint2, aes(x = reorder(tested_genome,-Complementarity), y=Complementarity)) + 
  geom_bar(position='dodge',stat = 'identity', col='black', size=0.1)+
  theme_classic()+
  theme(legend.title=element_blank(),axis.title.x=element_blank(),
        panel.spacing.x = unit(0,"lines"))+
  theme(legend.position="bottom")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  facet_wrap(~reference_genome, scales='free_x')


#####################Figure 5C
net <- read.csv("network4_PTM_network.tsv", sep = '\t', header=TRUE, stringsAsFactors = FALSE)
attr <- read.csv("network4_PTM_attributes.tsv", sep = '\t', header=TRUE, stringsAsFactors = FALSE)
attr <- attr[,c(5,6)]

# Identify species and metabolites based on attributes
species <- unique(c(net$node1[net$attributes == "secrete"], net$node2[net$attributes == "uptake"]))
metabolites <- unique(c(net$node2[net$attributes == "secrete"], net$node1[net$attributes == "uptake"]))

# Create an empty interaction matrix
interaction_matrix <- matrix(0, nrow = length(species), ncol = length(metabolites),
                             dimnames = list(species, metabolites))

# Populate the matrix based on secretion or uptake
for (i in 1:nrow(net)) {
  species_name <- ifelse(net$attributes[i] == "secrete", net$node1[i], net$node2[i])
  metabolite_name <- ifelse(net$attributes[i] == "secrete", net$node2[i], net$node1[i])
  
  if (species_name %in% species && metabolite_name %in% metabolites) {
    if (net$attributes[i] == "secrete") {
      interaction_matrix[species_name, metabolite_name] <- 
        ifelse(interaction_matrix[species_name, metabolite_name] == -1, 2, 1)
    } else if (net$attributes[i] == "uptake") {
      interaction_matrix[species_name, metabolite_name] <- 
        ifelse(interaction_matrix[species_name, metabolite_name] == 1, 2, -1)
    }
  }
}

# Merge metabolite superclasses with column names
metabolite_annotations <- data.frame(metabolite_name = colnames(interaction_matrix))
metabolite_annotations <- merge(metabolite_annotations, attr, by = "metabolite_name", all.x = TRUE)

# Order metabolites by superclass
metabolite_annotations <- metabolite_annotations[order(metabolite_annotations$Super.class), ]

# Reorder the interaction matrix columns based on ordered metabolites
interaction_matrix <- interaction_matrix[, metabolite_annotations$metabolite_name]

# Set superclasses as row names for annotation
row.names(metabolite_annotations) <- metabolite_annotations$metabolite_name
metabolite_annotations$metabolite_name <- NULL  # Remove unnecessary column

# Define a color palette for superclasses
superclass_colors <- setNames(rainbow(length(unique(metabolite_annotations$Super.class))), 
                              unique(metabolite_annotations$Super.class))

# Create annotation object
annotation_col <- list(Super.class = superclass_colors)

# Define a custom color palette for the heatmap
heatmap_colors <- c("-1" = "blue", "0" = "white", "1" = "red", "2" = "purple")

#Figure 5C (further customized in Adobe Illustrator)
pheatmap(interaction_matrix,
         color = unname(heatmap_colors), 
         breaks = c(-1.5, -0.5, 0.5, 1.5, 2.5),  # Define color breakpoints
         cluster_rows = F, cluster_cols = T,  # Keep columns ordered by superclass
         annotation_col = metabolite_annotations,
         annotation_colors = annotation_col,
         treeheight_col = 0,
         main = "Species-Metabolite Interaction Heatmap (Ordered by Superclass)",
         legend_labels = c("Uptake (-1)", "No Interaction (0)", "Secretion (+1)", "Both (+1 & -1)"),
         fontsize = 7)

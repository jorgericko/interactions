#Picrust


# Install the devtools package if not already installed

#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")

#pkgs <- c("ALDEx2","lefser")

#for (pkg in pkgs) {
# if (!requireNamespace(pkg, quietly = TRUE))
#  BiocManager::install(pkg)}

#install.packages("ggpicrust2")
# If you want to analyze the abundance of KEGG pathways instead of KO within the pathway, please set `ko_to_kegg` to TRUE.
# KEGG pathways typically have more descriptive explanations.


library(readr)
library(ggpicrust2)
library(tibble)
library(tidyverse)
library(ggprism)
library(patchwork)


library(tidyr)
library(dplyr)

MetaCyc_pathway_map

load("MetaCyc_pathway_map.RData") #Metacyc Pathways ontology

setwd("~/Projects/BETO/Microbial Analysis/2024_Interactions/R_Picrust")

# Load necessary data: abundance data and metadata
abundance_file <- "../Picrust/pathabun_exported/feature-table.biom.tsv" #path


# Run ggpicrust2 with imported data.frame

abundance_data <- read_delim(abundance_file, delim = "\t", col_names = TRUE, trim_ws = TRUE)

names(abundance_data)[names(abundance_data) == "OUT_ID"] <- "pathway"

#data("metacyc_abundance")
#data("metadata")
colnames(abundance_data) #Filter blank, food waste and manure samples
abundance_data <- abundance_data[,-c(339,340,341)]
colnames(abundance_data)





#### older code for reference


fermentation_pathways <- read.csv("fermentation_pathways.csv")

# RWAD_Major <- read.csv("RWAD_Major.csv")
RWAD_Major <- read.csv("RWAD_path.csv")

RWAD_Major_expanded <- RWAD_Major %>%
  separate_rows(pathway, sep = " // ") %>%
  distinct()  # Ensure unique pathway-process pairs





metadata <- read_delim(
  "DeLong_metadata_final.txt",
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE)

metadata <- metadata[-1, ]
names(metadata)[names(metadata) == "sample-id"] <- "sample_id"


metadata$T_pH_I <- paste(metadata$temperature, metadata$pH,metadata$inoculum, sep = "_")
metadata$T_pH_I_D <- paste(metadata$temperature, metadata$pH, metadata$inoculum, metadata$day, sep = "_")
metadata$T_pH_D <- paste(metadata$temperature, metadata$pH,metadata$day, sep = "_")
metadata$D_pH_T <- paste(metadata$day, metadata$pH,metadata$temperature, sep = "_")
metadata$pH_T_I <- paste(metadata$pH, metadata$temperature,metadata$inoculum, sep = "_")
metadata$D_pH <- paste(metadata$day,metadata$pH, sep = "_")
metadata$D_pH_T <- paste(metadata$day,metadata$pH, metadata$temperature, sep = "_")
metadata$D_pH_T_I <- paste(metadata$day,metadata$pH, metadata$temperature, metadata$inoculum, sep = "_")
metadata$D_pH_I <- paste(metadata$day,metadata$pH,metadata$inoculum, sep = "_")

metadata_interactions <- metadata[metadata$sample_id %in% colnames(abundance_data), ]

metadata_interactions_clean <- metadata_interactions %>% filter(temperature != "55",day != "D0")


metadata_F <- metadata_interactions_clean %>% filter(feedstock == "F")
metadata_M <- metadata_interactions_clean %>% filter(feedstock == "M")

abundance_F <- abundance_data[, c("pathway", names(abundance_data)[names(abundance_data) %in% metadata_F$sample_id & names(abundance_data) != "pathway"])]
abundance_M <- abundance_data[, c("pathway", names(abundance_data)[names(abundance_data) %in% metadata_M$sample_id & names(abundance_data) != "pathway"])]


library(dplyr)
library(tibble)

# First, ensure 'pathway' is not removed or altered
# Remove non-numeric columns (except 'pathway'), check variance, then convert 'pathway' to row names
abundance_F_non0 <- abundance_F %>%
  select(pathway, where(is.numeric)) %>%
  mutate(across(where(is.numeric), ~ ifelse(is.na(.), 0, .))) %>%  # Replace NA with 0 in numeric columns
  filter_if(where(is.numeric), any_vars(. != 0 & !is.na(.))) 

abundance_M_non0 <- abundance_M %>%
  select(pathway, where(is.numeric)) %>%
  mutate(across(where(is.numeric), ~ ifelse(is.na(.), 0, .))) %>%  # Replace NA with 0 in numeric columns
  filter_if(where(is.numeric), any_vars(. != 0 & !is.na(.)))



abundance_non0 <- abundance_data %>%
  select(pathway, where(is.numeric)) %>%
  mutate(across(where(is.numeric), ~ ifelse(is.na(.), 0, .))) %>%  # Replace NA with 0 in numeric columns
  filter_if(where(is.numeric), any_vars(. != 0 & !is.na(.)))



RWAD_Major_abundance_F <- abundance_F_non0 %>%
  inner_join(RWAD_Major_expanded, by = "pathway") %>%
  group_by(process) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE))

RWAD_Major_abundance_M <- abundance_M_non0 %>%
  inner_join(RWAD_Major_expanded, by = "pathway") %>%
  group_by(process) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE))


RWAD_Major_abundance <- abundance_non0 %>%
  inner_join(RWAD_Major_expanded, by = "pathway") %>%
  group_by(process) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE))


library(ggh4x)

help("ggpicrust2")



fermentation_MetaCyc <- MetaCyc_pathway_map %>%
  # Filter for pathways that are in fermentation_pathways OR have Superclass2 as Fermentation
  filter(pathway %in% fermentation_pathways$Fermentation | Superclass2 == "Fermentation") %>%
  # Remove duplicate rows, if any, after filtering
  distinct()




#Explore the data with different pcoa
#pathway_pca(abundance = abundance_M_non0  %>% column_to_rownames("pathway") , metadata = metadata_M, group = "inoculum_pH")

M_metacyc_daa_results_df <- pathway_daa(abundance = abundance_M_non0%>% column_to_rownames("pathway"), metadata = metadata_M, group = "pH_T_I", daa_method = "LinDA", reference = "7_35_A")

M_metacyc_daa_annotated_results_dfs <- pathway_annotation(pathway = "MetaCyc", daa_results_df = M_metacyc_daa_results_df, ko_to_kegg = FALSE)




F_metacyc_daa_results_df <- pathway_daa(abundance = abundance_F_non0%>% column_to_rownames("pathway"), metadata = metadata_F, group = "pH_T_I", daa_method = "LinDA", reference = "7_35_A")
F_metacyc_daa_annotated_results_dfs <- pathway_annotation(pathway = "MetaCyc", daa_results_df = F_metacyc_daa_results_df, ko_to_kegg = FALSE)



metacyc_daa_results_df <- pathway_daa(abundance = abundance_non0%>% column_to_rownames("pathway"), metadata = metadata, group = "feedstock_pH", daa_method = "LinDA", reference = "F_7")

metacyc_daa_annotated_results_dfs <- pathway_annotation(pathway = "MetaCyc", daa_results_df = metacyc_daa_results_df, ko_to_kegg = FALSE)


#F_metacyc_with_p_0.05 <- F_metacyc_daa_annotated_results_dfs %>% filter(p_adjust < 0.05)

#filtered_abundance_data_0 <- filtered_abundance_data[filtered_abundance_data$OUT_ID %in% filtered_data$feature, ]
F_metacyc_with_p_0.05 <- F_metacyc_daa_annotated_results_dfs %>% filter(p_adjust < 0.05)
M_metacyc_with_p_0.05 <- M_metacyc_daa_annotated_results_dfs %>% filter(p_adjust < 0.05)



# pathway_errorbar(abundance = abundance_F_non0,
#                  daa_results_df = F_metacyc_daa_annotated_results_dfs,
#                  Group = metadata$temperature,
#                  ko_to_kegg = TRUE,
#                  p_values_threshold = 0.05,
#                  order = "pathway_class",
#                  select = NULL,
#                  p_value_bar = TRUE,
#                  colors = NULL,
#                  x_lab = NULL)



#Keep entries for feature column in df1 that match rownames in df2 assigned to entries in superclass2 in df2 

library(dplyr)
library(tibble)


metadata_F$pH_T_I <- paste(metadata_F$pH, metadata_F$temperature,metadata_F$inoculum, sep = "_")
metadata_M$pH_T_I <- paste(metadata_M$pH, metadata_M$temperature,metadata_M$inoculum, sep = "_")


#Food waste
Carbohydrate_Degradation_F <- F_metacyc_daa_annotated_results_dfs %>%
  inner_join(tibble::rownames_to_column(MetaCyc_pathway_map %>% filter(Superclass2 == "Carbohydrate Degradation"), var = "feature"), by = "feature")

Amino_Acid_Degradation_F <- F_metacyc_daa_annotated_results_dfs %>%
  inner_join(tibble::rownames_to_column(MetaCyc_pathway_map %>% filter(Superclass2 == "Amino Acid Degradation"), var = "feature"), by = "feature")

Fermentation_F <- F_metacyc_daa_annotated_results_dfs %>%
  inner_join(tibble::rownames_to_column(fermentation_MetaCyc, var = "feature"), by = "feature")

FA_Lipid_bios_F <- F_metacyc_daa_annotated_results_dfs %>%
  inner_join(tibble::rownames_to_column(MetaCyc_pathway_map %>% filter(Superclass2 == "Fatty Acid and Lipid Biosynthesis"), var = "feature"), by = "feature")

Aromatic_deg_F <- F_metacyc_daa_annotated_results_dfs %>%
  inner_join(tibble::rownames_to_column(MetaCyc_pathway_map %>% filter(Superclass2 == "Aromatic Compound Degradation"), var = "feature"), by = "feature")



#Manure
Carbohydrate_Degradation_M <- M_metacyc_daa_annotated_results_dfs %>%
  inner_join(tibble::rownames_to_column(MetaCyc_pathway_map %>% filter(Superclass2 == "Carbohydrate Degradation"), var = "feature"), by = "feature")

Amino_Acid_Degradation_M <- M_metacyc_daa_annotated_results_dfs %>%
  inner_join(tibble::rownames_to_column(MetaCyc_pathway_map %>% filter(Superclass2 == "Amino Acid Degradation"), var = "feature"), by = "feature")

Fermentation_M <- M_metacyc_daa_annotated_results_dfs %>%
  inner_join(tibble::rownames_to_column(fermentation_MetaCyc, var = "feature"), by = "feature")

FA_Lipid_bios_M <- M_metacyc_daa_annotated_results_dfs %>%
  inner_join(tibble::rownames_to_column(MetaCyc_pathway_map %>% filter(Superclass2 == "Fatty Acid and Lipid Biosynthesis"), var = "feature"), by = "feature")


Aromatic_deg_M <- M_metacyc_daa_annotated_results_dfs %>%
  inner_join(tibble::rownames_to_column(MetaCyc_pathway_map %>% filter(Superclass2 == "Aromatic Compound Degradation"), var = "feature"), by = "feature")



#all
Carbohydrate_Degradation <- metacyc_daa_annotated_results_dfs %>%
  inner_join(tibble::rownames_to_column(MetaCyc_pathway_map %>% filter(Superclass2 == "Carbohydrate Degradation"), var = "feature"), by = "feature")

Amino_Acid_Degradation <- metacyc_daa_annotated_results_dfs %>%
  inner_join(tibble::rownames_to_column(MetaCyc_pathway_map %>% filter(Superclass2 == "Amino Acid Degradation"), var = "feature"), by = "feature")

Fermentation <- metacyc_daa_annotated_results_dfs %>%
  inner_join(tibble::rownames_to_column(fermentation_MetaCyc, var = "feature"), by = "feature")

FA_Lipid_bios <- metacyc_daa_annotated_results_dfs %>%
  inner_join(tibble::rownames_to_column(MetaCyc_pathway_map %>% filter(Superclass2 == "Fatty Acid and Lipid Biosynthesis"), var = "feature"), by = "feature")

FA_Lipid_deg <- metacyc_daa_annotated_results_dfs %>%
   inner_join(tibble::rownames_to_column(MetaCyc_pathway_map %>% filter(Superclass2 == "Fatty Acid and Lipid Degradation"), var = "feature"), by = "feature")

Aromatic_deg <- metacyc_daa_annotated_results_dfs %>%
  inner_join(tibble::rownames_to_column(MetaCyc_pathway_map %>% filter(Superclass2 == "Aromatic Compound Degradation"), var = "feature"), by = "feature")




pathway_heatmap(abundance = abundance_non0 %>% right_join(Aromatic_deg %>% select(all_of(c("feature", "description"))) %>% distinct(), by = c("pathway" = "feature")) %>% select(-"pathway") %>% column_to_rownames("description"), metadata = metadata, group = "feedstock_pH")
pathway_heatmap(abundance = abundance_M_non0 %>% right_join(Aromatic_deg %>% select(all_of(c("feature", "description"))) %>% distinct(), by = c("pathway" = "feature")) %>% select(-"pathway") %>% column_to_rownames("description"), metadata = metadata_M, group = "pH_T_I")



# pathway_heatmap(abundance = abundance_F_non0 %>% right_join(Carbohydrate_Degradation_F %>% select(all_of(c("feature", "description"))) %>% distinct(), by = c("pathway" = "feature")) %>% select(-"pathway") %>% column_to_rownames("description"), metadata = metadata_F, group = "pH_T_I")

pdf("Aromatic_deg.pdf", width = 16, height =6) # Set width and height
pathway_heatmap(abundance = abundance_non0 %>% right_join(Aromatic_deg %>% select(all_of(c("feature", "description"))) %>% distinct(), by = c("pathway" = "feature")) %>% select(-"pathway") %>% column_to_rownames("description"), metadata = metadata, group = "feedstock_pH")
dev.off()


pdf("M_Aromatic_deg.pdf", width = 16, height =6) # Set width and height
pathway_heatmap(abundance = abundance_M_non0 %>% right_join(Aromatic_deg %>% select(all_of(c("feature", "description"))) %>% distinct(), by = c("pathway" = "feature")) %>% select(-"pathway") %>% column_to_rownames("description"), metadata = metadata_M, group = "pH_T_I")
dev.off()


pdf("F_Aromatic_deg.pdf", width = 16, height =6) # Set width and height
pathway_heatmap(abundance = abundance_F_non0 %>% right_join(Aromatic_deg %>% select(all_of(c("feature", "description"))) %>% distinct(), by = c("pathway" = "feature")) %>% select(-"pathway") %>% column_to_rownames("description"), metadata = metadata_F, group = "pH_T_I")
dev.off()



pdf("Carbohydrate_Degradation.pdf", width = 16, height =4.5) # Set width and height
pathway_heatmap(abundance = abundance_non0 %>% right_join(Carbohydrate_Degradation %>% select(all_of(c("feature", "description"))) %>% distinct(), by = c("pathway" = "feature")) %>% select(-"pathway") %>% column_to_rownames("description"), metadata = metadata, group = "feedstock_pH")
dev.off()

pdf("Amino_Acid_Degradation.pdf", width = 16, height =4.5) # Set width and height
pathway_heatmap(abundance = abundance_non0 %>% right_join(Amino_Acid_Degradation %>% select(all_of(c("feature", "description"))) %>% distinct(), by = c("pathway" = "feature")) %>% select(-"pathway") %>% column_to_rownames("description"), metadata = metadata, group = "feedstock_pH")
dev.off()

pdf("Fermentation.pdf", width = 16, height =4.5) # Set width and height
pathway_heatmap(abundance = abundance_non0 %>% right_join(Fermentation %>% select(all_of(c("feature", "description"))) %>% distinct(), by = c("pathway" = "feature")) %>% select(-"pathway") %>% column_to_rownames("description"), metadata = metadata, group = "feedstock_pH")
dev.off()

pdf("FA_Lipid_bios.pdf", width = 16, height =4.5) # Set width and height
pathway_heatmap(abundance = abundance_non0 %>% right_join(FA_Lipid_bios %>% select(all_of(c("feature", "description"))) %>% distinct(), by = c("pathway" = "feature")) %>% select(-"pathway") %>% column_to_rownames("description"), metadata = metadata, group = "feedstock_pH")
dev.off()



#print heatmap plots

pdf("Carbohydrate_Degradation_F.pdf", width = 16, height =4.5) # Set width and height
pathway_heatmap(abundance = abundance_F_non0 %>% right_join(Carbohydrate_Degradation_F %>% select(all_of(c("feature", "description"))) %>% distinct(), by = c("pathway" = "feature")) %>% select(-"pathway") %>% column_to_rownames("description"), metadata = metadata_F, group = "pH_T_I")
dev.off()

pdf("Amino_Acid_Degradation_F.pdf", width = 16, height =4.5) # Set width and height
pathway_heatmap(abundance = abundance_F_non0 %>% right_join(Amino_Acid_Degradation_F %>% select(all_of(c("feature", "description"))) %>% distinct(), by = c("pathway" = "feature")) %>% select(-"pathway") %>% column_to_rownames("description"), metadata = metadata_F, group = "pH_T_I")
dev.off()

pdf("Fermentation_F.pdf", width = 16, height =4.5) # Set width and height
pathway_heatmap(abundance = abundance_F_non0 %>% right_join(Fermentation_F %>% select(all_of(c("feature", "description"))) %>% distinct(), by = c("pathway" = "feature")) %>% select(-"pathway") %>% column_to_rownames("description"), metadata = metadata_F, group = "pH_T_I")
dev.off()

pdf("FA_Lipid_bios_F.pdf", width = 16, height =4.5) # Set width and height
pathway_heatmap(abundance = abundance_F_non0 %>% right_join(FA_Lipid_bios_F %>% select(all_of(c("feature", "description"))) %>% distinct(), by = c("pathway" = "feature")) %>% select(-"pathway") %>% column_to_rownames("description"), metadata = metadata_F, group = "pH_T_I")
dev.off()


pdf("Carbohydrate_Degradation_M.pdf", width = 16, height =4.5) # Set width and height
pathway_heatmap(abundance = abundance_M_non0 %>% right_join(Carbohydrate_Degradation_M %>% select(all_of(c("feature", "description"))) %>% distinct(), by = c("pathway" = "feature")) %>% select(-"pathway") %>% column_to_rownames("description"), metadata = metadata_M, group = "pH_T_I")
dev.off()

pdf("Amino_Acid_Degradation_M.pdf", width = 16, height =4.5) # Set width and height
pathway_heatmap(abundance = abundance_M_non0 %>% right_join(Amino_Acid_Degradation_M %>% select(all_of(c("feature", "description"))) %>% distinct(), by = c("pathway" = "feature")) %>% select(-"pathway") %>% column_to_rownames("description"), metadata = metadata_M, group = "pH_T_I")
dev.off()

pdf("Fermentation_M.pdf", width = 16, height =4.5) # Set width and height
pathway_heatmap(abundance = abundance_M_non0 %>% right_join(Fermentation_M %>% select(all_of(c("feature", "description"))) %>% distinct(), by = c("pathway" = "feature")) %>% select(-"pathway") %>% column_to_rownames("description"), metadata = metadata_M, group = "pH_T_I")
dev.off()

pdf("FA_Lipid_bios_M.pdf", width = 16, height =4.5) # Set width and height
pathway_heatmap(abundance = abundance_M_non0 %>% right_join(FA_Lipid_bios_M %>% select(all_of(c("feature", "description"))) %>% distinct(), by = c("pathway" = "feature")) %>% select(-"pathway") %>% column_to_rownames("description"), metadata = metadata_M, group = "pH_T_I")
dev.off()



#Print major
pdf("RWAD_Major_abundance_F.pdf", width = 16, height =6) # Set width and height
pathway_heatmap(abundance = RWAD_Major_abundance_F %>% column_to_rownames("process"), metadata = metadata_F, group = "pH_T_I")
dev.off()

pdf("RWAD_Major_abundance_M.pdf", width = 16, height =6) # Set width and height
pathway_heatmap(abundance = RWAD_Major_abundance_M %>% column_to_rownames("process"), metadata = metadata_M, group = "pH_T_I")
dev.off()


pdf("RWAD_Major_abundance.pdf", width = 16, height =6) # Set width and height
pathway_heatmap(abundance = RWAD_Major_abundance %>% column_to_rownames("process"), metadata = metadata, group = "feedstock_pH")
dev.off()


# 
# 
# Carbohydrate_Degradation_F <- F_metacyc_with_p_0.05 %>%
#   inner_join(tibble::rownames_to_column(MetaCyc_pathway_map %>% filter(Superclass2 == "Carbohydrate Degradation"), var = "feature"), by = "feature")
# 
# Amino_Acid_Degradation_F <- F_metacyc_with_p_0.05 %>%
#   inner_join(tibble::rownames_to_column(MetaCyc_pathway_map %>% filter(Superclass2 == "Amino Acid Degradation"), var = "feature"), by = "feature")
# 
# Fermentation_F <- F_metacyc_with_p_0.05 %>%
#   inner_join(tibble::rownames_to_column(fermentation_MetaCyc, var = "feature"), by = "feature")
# 
# 
# Fermentation_M <- M_metacyc_daa_annotated_results_dfs %>%
#   inner_join(tibble::rownames_to_column(fermentation_MetaCyc, var = "feature"), by = "feature")
# 
# 
# Respiration_F <- F_metacyc_with_p_0.05 %>%
#   inner_join(tibble::rownames_to_column(MetaCyc_pathway_map %>% filter(Superclass2 == "Respiration"), var = "feature"), by = "feature")
# 
# Carboxylate_Degradation_F <- F_metacyc_with_p_0.05 %>%
#   inner_join(tibble::rownames_to_column(MetaCyc_pathway_map %>% filter(Superclass2 == "Carboxylate Degradation"), var = "feature"), by = "feature")
# 
# FA_Lipid_bios_F <- F_metacyc_with_p_0.05 %>%
#   inner_join(tibble::rownames_to_column(MetaCyc_pathway_map %>% filter(Superclass2 == "Fatty Acid and Lipid Biosynthesis"), var = "feature"), by = "feature")
# 
# FA_Lipid_deg_F <- F_metacyc_with_p_0.05 %>%
#   inner_join(tibble::rownames_to_column(MetaCyc_pathway_map %>% filter(Superclass2 == "Fatty Acid and Lipid Degradation"), var = "feature"), by = "feature")
# 
# Carbohydrate_Degradation_F <- F_metacyc_daa_annotated_results_df %>%
#   inner_join(tibble::rownames_to_column(MetaCyc_pathway_map %>% filter(Superclass2 == "Carbohydrate Degradation"), var = "feature"), by = "feature")
# 
# A_CoA_F <- F_metacyc_daa_annotated_results_df %>%
#   inner_join(tibble::rownames_to_column(MetaCyc_pathway_map %>% filter(Superclass2 == "Acetyl-CoA Biosynthesis"), var = "feature"), by = "feature")
# 
# 
# 
# 
# Cofactor_F <- F_metacyc_with_p_0.05 %>%
#   inner_join(tibble::rownames_to_column(MetaCyc_pathway_map %>% filter(Superclass2 == "Cofactor, Carrier, and Vitamin Biosynthesis"), var = "feature"), by = "feature")
# 
# E_chains <- F_metacyc_with_p_0.05 %>%
#   inner_join(tibble::rownames_to_column(MetaCyc_pathway_map %>% filter(Superclass2 == "Electron Transfer Chains"), var = "feature"), by = "feature")
# 
# alcohol_deg <- F_metacyc_with_p_0.05 %>%
#   inner_join(tibble::rownames_to_column(MetaCyc_pathway_map %>% filter(Superclass2 == "Alcohol Degradation"), var = "feature"), by = "feature")
# 
# detox <- F_metacyc_with_p_0.05 %>%
#   inner_join(tibble::rownames_to_column(MetaCyc_pathway_map %>% filter(Superclass1 == "Detoxification"), var = "feature"), by = "feature")
# 
# Glycolysis <- F_metacyc_with_p_0.05 %>%
#   inner_join(tibble::rownames_to_column(MetaCyc_pathway_map %>% filter(Superclass2 == "Glycolysis"), var = "feature"), by = "feature")
# 
# TCA <- F_metacyc_daa_annotated_results_df %>%
#   inner_join(tibble::rownames_to_column(MetaCyc_pathway_map %>% filter(Superclass2 == "TCA cycle"), var = "feature"), by = "feature")



#######################
## MODULES ANNOTATION
#######################

library(tidyr)
library(dplyr)
library(readr)  # For exporting to CSV
# If exporting to Excel, use the writexl package: install.packages("writexl")
library(writexl)  # For exporting to Excel




### Read AD database
RWAD_path <- read.csv("RWAD_path.csv")

RWAD_path_expanded <- RWAD_path %>%
  separate_rows(pathway, sep = " // ") %>%
  distinct()  # Ensure unique pathway-process pairs


# Read path coverage from picrust2
file_path_path_coverage <- "~/Projects/BETO/Microbial Analysis/2024_Interactions/Picrust/picrust2_py/picrust2_out_pipeline_per_seq_con/pathways_out/path_cov_predictions.tsv"
path_cov_predictions <- read.delim(file_path_path_coverage, header = TRUE, sep = "\t", check.names = FALSE)

#read taxonomy
taxonomy_table<-read.table("~/Projects/BETO/Microbial Analysis/2024_Interactions/16S/taxonomy_table_2.txt", sep='\t', header=T, comment="")

#read module data

fw_node_table <- read.csv("~/Projects/BETO/Microbial Analysis/2024_Interactions/16S/fw_node_table_no55.csv", row.names = 1)
m_node_table <- read.csv("~/Projects/BETO/Microbial Analysis/2024_Interactions/16S/m_node_table_no55.csv", row.names = 1)


#rename first column to merge three datasets if needed
colnames(path_cov_predictions)[1] <- "feature"
colnames(taxonomy_table)[1] <- "feature"
colnames(fw_node_table)[1] <- "feature"
colnames(m_node_table)[1] <- "feature"

# Merge taxonomy and path coverage data frames based on the 'feature' column

path_cov_taxonomy <- merge(path_cov_predictions, taxonomy_table, by = "feature", all.x = TRUE)

# transform into long format to group by genus later
path_cov_tax_long <- pivot_longer(
  path_cov_taxonomy,
  cols = -c(feature, Kingdom, Phylum, Class, Family, Order, Genus, Species),  # All columns except the specified ones
  names_to = "pathway",  # Name of new column for pathways
  values_to = "value"    # Name of new column for pathway values
)


#Sum pathway coverage by genus and pathway (Get mean path coverage by genus)
genus_pathway_mean <- path_cov_tax_long %>%
  group_by(Phylum,Order,Family,Genus, pathway) %>%
  summarise(mean_coverage = mean(value, na.rm = TRUE)) %>%
  ungroup() 



#wide version if neeeded
# genus_pathway_wide <- genus_pathway_mean %>%
#   pivot_wider(names_from = pathway, values_from = mean_coverage)

# Merge node (template) and genus_path data data frames based on the 'Genus' column
node_pathways <- merge(fw_node_table, genus_pathway_mean, by = "Genus", all.x = TRUE)

m_node_pathways <- merge(m_node_table, genus_pathway_mean, by = "Genus", all.x = TRUE)



# Calculate weighted_coverage by multiplying Abundance and mean_coverage
node_pathways <- node_pathways %>%
  mutate(weighted_coverage = Abundance * mean_coverage)

m_node_pathways <- m_node_pathways %>%
  mutate(weighted_coverage = Abundance * mean_coverage)


# Group pathways by processes given by the RWAD database (sum)

#Annotated RWAD network by module
network_module_path_RWAD_weighted <- node_pathways %>%
  inner_join(RWAD_path_expanded, by = "pathway",relationship = "many-to-many") %>%
  group_by(process,module) %>%
  summarise(across(weighted_coverage, sum, na.rm = TRUE), .groups = "drop")

m_network_module_path_RWAD_weighted <- m_node_pathways %>%
  inner_join(RWAD_path_expanded, by = "pathway",relationship = "many-to-many") %>%
  group_by(process,module) %>%
  summarise(across(weighted_coverage, sum, na.rm = TRUE), .groups = "drop")





network_module_path_RWAD_log <- network_module_path_RWAD_weighted %>% 
  group_by(process) %>% 
  mutate(log_transformed = log(weighted_coverage + 1)) %>%
  ungroup()%>% 
  select(-weighted_coverage)  # This removes the 'weighted_coverage' column

m_network_module_path_RWAD_log <- m_network_module_path_RWAD_weighted %>% 
  group_by(process) %>% 
  mutate(log_transformed = log(weighted_coverage + 1)) %>%
  ungroup()%>% 
  select(-weighted_coverage)  # This removes the 'weighted_coverage' column



network_module_path_RWAD_log_wide <- network_module_path_RWAD_log %>%
  pivot_wider(names_from = process, values_from = log_transformed)

m_network_module_path_RWAD_log_wide <- m_network_module_path_RWAD_log %>%
  pivot_wider(names_from = process, values_from = log_transformed)


# Genomic potential per taxa per module (weighted)
network_genus_path_RWAD_weighted <- node_pathways %>%
  inner_join(RWAD_path_expanded, by = "pathway",relationship = "many-to-many") %>%
  group_by(process,module,Genus,degree, betweenness_centrality,closeness_centrality,eigenvector_centrality,Abundance,z,p,taxa_roles) %>%
  summarise(across(weighted_coverage, sum, na.rm = TRUE), .groups = "drop")


network_genus_path_RWAD_log <- network_genus_path_RWAD_weighted %>%
  group_by(process,module,Genus,degree, betweenness_centrality,closeness_centrality,eigenvector_centrality,Abundance,z,p,taxa_roles) %>%
  mutate(log_transformed = log(weighted_coverage + 1)) %>%
  ungroup()%>% 
  select(-weighted_coverage)  # This removes the 'weighted_coverage' column


m_network_genus_path_RWAD_weighted <- m_node_pathways %>%
  inner_join(RWAD_path_expanded, by = "pathway",relationship = "many-to-many") %>%
  group_by(process,module,Genus,degree, betweenness_centrality,closeness_centrality,eigenvector_centrality,Abundance,z,p,taxa_roles) %>%
  summarise(across(weighted_coverage, sum, na.rm = TRUE), .groups = "drop")


m_network_genus_path_RWAD_log <- m_network_genus_path_RWAD_weighted %>%
  group_by(process,module,Genus,degree, betweenness_centrality,closeness_centrality,eigenvector_centrality,Abundance,z,p,taxa_roles) %>%
  mutate(log_transformed = log(weighted_coverage + 1)) %>%
  ungroup()%>% 
  select(-weighted_coverage)  # This removes the 'weighted_coverage' column



network_genus_path_RWAD_log_wide <- network_genus_path_RWAD_log %>%
  pivot_wider(names_from = process, values_from = log_transformed)

m_network_genus_path_RWAD_log_wide <- m_network_genus_path_RWAD_log %>%
  pivot_wider(names_from = process, values_from = log_transformed)

library(writexl)

# Export to
write_xlsx(network_module_path_RWAD_log_wide, "annotated_fw_network_module_no55.xlsx")

write_xlsx(network_genus_path_RWAD_log_wide, "annotated_fw_network_genus_no55.xlsx")

# Export to CSV
write_xlsx(m_network_module_path_RWAD_log_wide, "annotated_m_network_module_no55.xlsx")

# OR, export to Excel
write_xlsx(m_network_genus_path_RWAD_log_wide, "annotated_m_network_genus_no55.xlsx")



# Load the necessary library
library(writexl)

# Create a named list of your dataframes
df_list <- list(
  "FW_Network_Module" = network_module_path_RWAD_log_wide,
  "FW_Network_Genus" = network_genus_path_RWAD_log_wide,
  "M_Network_Module" = m_network_module_path_RWAD_log_wide,
  "M_Network_Genus" = m_network_genus_path_RWAD_log_wide
)

# Export all data frames to one Excel file with multiple sheets
write_xlsx(df_list, "annotated_network_data_no55.xlsx")



#### Read Maaslin2 tables

maaslin2_M_wide <- read.csv("~/Projects/BETO/Microbial Analysis/2024_Interactions/16S/maaslin2_M_wide.csv")

maaslin2_F_wide <- read.csv("~/Projects/BETO/Microbial Analysis/2024_Interactions/16S/maaslin2_F_wide.csv")

library(dplyr)
m_network_genus_path_RWAD_log_wide <- m_network_genus_path_RWAD_log_wide %>%
  rename(genus = Genus)

network_genus_path_RWAD_log_wide <- network_genus_path_RWAD_log_wide %>%
  rename(genus = Genus)

head(network_genus_path_RWAD_log_wide)

head(maaslin2_F_wide)

library(dplyr)
library(stringr)

# Remove leading/trailing whitespace and any double quotes
network_genus_path_RWAD_log_wide <- network_genus_path_RWAD_log_wide %>%
  mutate(genus = str_trim(gsub('"', '', genus)))

maaslin2_F_wide <- maaslin2_F_wide %>%
  mutate(genus = str_trim(gsub('"', '', genus)))


m_network_genus_path_RWAD_log_wide <- m_network_genus_path_RWAD_log_wide %>%
  mutate(genus = str_trim(gsub('"', '', genus)))

maaslin2_M_wide <- maaslin2_M_wide %>%
  mutate(genus = str_trim(gsub('"', '', genus)))


fw_maaslin2_picrust <- merge(network_genus_path_RWAD_log_wide, maaslin2_F_wide, by = "genus", all.x = TRUE)
head(fw_maaslin2_picrust)


m_maaslin2_picrust <- merge(m_network_genus_path_RWAD_log_wide, maaslin2_M_wide, by = "genus", all.x = TRUE)


maaslin2_list <- list(
  "fw_maaslin2_picrust" = fw_maaslin2_picrust,
  "m_maaslin2_picrust" = m_maaslin2_picrust
)

# Export all data frames to one Excel file with multiple sheets
write_xlsx(maaslin2_list, "maaslin2_picrust2_network.xlsx")




install.packages("corrplot")

#
# Subset numeric columns from columns 3 to 43
subset_data <- fw_maaslin2_picrust[, 3:42]
numeric_subset_data <- subset_data[, sapply(subset_data, is.numeric)]

# Calculate the covariance matrix
covariance_matrix <- cov(numeric_subset_data, use = "complete.obs")
covariance_matrix

library(corrplot)

# Plot the covariance matrix
corrplot(covariance_matrix, method = "color", type = "upper",
         tl.col = "black", tl.srt = 45,
         col = colorRampPalette(c("blue", "white", "red"))(200),
         title = "Covariance Matrix Plot", mar = c(0, 0, 2, 0))



# Install and load the necessary package if not already installed
if (!require("corrplot")) install.packages("corrplot")
library(corrplot)

# Subset numeric columns from columns 3 to 43
subset_data <- fw_maaslin2_picrust[, 3:42]
numeric_subset_data <- subset_data[, sapply(subset_data, is.numeric)]

# Calculate the covariance matrix
covariance_matrix <- cov(numeric_subset_data, use = "complete.obs")

# Convert the covariance matrix to a correlation matrix
correlation_matrix <- cov2cor(covariance_matrix)


# Shorten the column and row names to the first 5 characters
colnames(correlation_matrix) <- substr(colnames(correlation_matrix), 1, 15)
rownames(correlation_matrix) <- substr(rownames(correlation_matrix), 1, 15)

# Plot the correlation matrix with shortened labels
corrplot(correlation_matrix, method = "color", type = "upper",
         tl.col = "black", tl.srt = 45, tl.cex = 0.5,  # Adjust text size if needed
         col = colorRampPalette(c("blue", "white", "red"))(200),
         title = "Correlation Matrix Plot", mar = c(0, 0, 2, 0))



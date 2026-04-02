#Microbiome analysis interactions study, let's do it!

setwd("~/Projects/BETO/Microbial Analysis/2024_Interactions/R_Picrust")


library(ggh4x)
library(dplyr)
library(tibble)
library(readr)
library(ggpicrust2)
library(tibble)
library(tidyverse)
library(ggprism)
library(patchwork)
library(tidyr)
library(dplyr)
library(dplyr)
library(tibble)
rm(list = ls())

library(microeco)
#load combined data
library(ggplot2)

library(dplyr)
# data_acids_0<-read.csv("data_acids_A_B.csv", comment="") #Combined data
# 
# data_acids_f1=filter(data_acids_0,rep !="D")
# data_acids_f2=filter(data_acids_f1,rep !="INOC")
# data_acids_f3=filter(data_acids_f2,day != 9.7)%>%
#   filter(between(day, 0, 13))  %>% filter(T!=55)
# 
# 
# data_acids <- pivot_wider(data_acids_f3, names_from = acid, values_from = conc_gL)
# data_acids[is.na(data_acids)] = 0
# data_acids$total <- rowSums(data_acids[,15:23])
# data_acids$day <- as.numeric(data_acids$day)
# 
# write.csv(data_acids, "data_acids_expanded.csv", row.names = FALSE)
# metadata <- read.table("Delong_metadata_final.txt", sep = "\t", header =T)
# metadata<-read.table("metadata.txt", sep='\t', header=T, row.names=1, comment="")
# data_acids_0<-read.csv("data_acids_A_B.csv", comment="") #Combined data


################################################################################
##### READ AND CLEAN TABLES
#################################################################################

# sample_table: rownames of sample_table must be sample names used
# otu_table: rownames must be feature names; colnames must be sample names
# microtable class: creating microtable object requires at least one file input (otu_table)
# tidy_taxonomy(): necessary to make taxonomic table have unified format
# tidy_dataset(): necessary to trim files in microtable object
# add_rownames2taxonomy(): add the rownames of tax_table as the last column of tax_table
# cal_abund(): powerful and flexible to cope with complex cases in tax_table, see the parameters
# taxa_abund: taxa_abund is a list stored in microtable object and have several data frame
# beta_diversity: beta_diversity is a list stored in microtable object and have several distance matri

metadata_0<-read.table("../16S/metadata_interactions.txt", sep='\t', header=T, row.names=1, comment="")
# taxonomy_table<-read.table("taxonomy_table_2.txt", sep='\t', header=T, row.names=1, comment="")
load("MetaCyc_pathway_map.RData") #Metacyc Pathways ontology
RWAD_Major <- read.csv("RWAD_path_clean_2.csv")
RWAD_Major_expanded_0 <- RWAD_Major %>%separate_rows(pathway, sep = " // ")
colnames(RWAD_Major_expanded_0)[colnames(RWAD_Major_expanded_0) == "pathway"] <- "path_id"


MetaCyc_pathway <- MetaCyc_pathway_map %>% tibble::rownames_to_column(var = "path_id")  

path_abun <- RWAD_Major_expanded_0 %>% left_join(MetaCyc_pathway, by = "path_id") 

# Convert to data frame if not already
RWAD_Major_expanded <- as.data.frame(path_abun)

# Create unique row names without modifying special characters
unique_pathways <- RWAD_Major_expanded$path_id

# Handle duplicates by appending numbers manually
unique_pathways <- make.unique(as.character(unique_pathways))

# Assign the modified unique names while keeping dashes
rownames(RWAD_Major_expanded) <- unique_pathways

# Remove the 'pathway' column to avoid redundancy
# RWAD_Major_expanded$pathway <- NULL

# View the result
head(RWAD_Major_expanded)





taxonomy_table<-RWAD_Major_expanded

taxonomy_table$Superclass1 <- NULL
taxonomy_table$Superclass2 <- NULL
taxonomy_table$path_id <- NULL
taxonomy_table$process <- NULL
taxonomy_table$type <- NULL


acids_env<-read.csv("../16S/data_acids_with_sid.csv", header=T, row.names=1, comment="") #Combined data

abundance_file <- "../picrust2_updated/picrust2_out_pipeline_per_seq_con/pathways_out/path_abun_unstrat_new_sid.txt" #path
# Load necessary library
library(readr)

asv_table <- read_delim(abundance_file, delim = "\t", col_names = TRUE)
asv_table <- as.data.frame(asv_table)
rownames(asv_table) <- asv_table[, 1]
asv_table <- asv_table[, -1]



acids_env<-as.data.frame(acids_env)
acids_env <- acids_env[ , !names(acids_env) %in% c("ISOC6","ISOC5","ISOC4")]


taxonomy_table<-as.data.frame(taxonomy_table)
asv_table<-as.data.frame(asv_table)
metadata_0<-as.data.frame(metadata_0)
samples_out <- rownames(filter(metadata_0,keep == "n"))

# Get the sample IDs from the otu_table
asv_table_ids <- colnames(asv_table)
metadata <-metadata_0[rownames(metadata_0) %in% asv_table_ids, ]
nrow(metadata)
length(asv_table_ids)

## clean taxonomy table
# taxonomy_table %<>% tidy_taxonomy
head(taxonomy_table)

# taxonomy_table<-read.table("taxonomy_table.txt", sep='\t', header=T, row.names=1, comment="")
# taxonomy_table %<>% tidy_taxonomy

# use the example data rep_phylo.tre in file2meco package https://chiliubio.github.io/microeco_tutorial/file2meco-package.html#qiime
# phylo_file_path <- system.file("extdata", "rep_phylo.tre", package="file2meco")
# tree <- ape::read.tree("rooted-tree/tree.nwk")

#create microeco environment

metadata_clean <- filter(metadata_0,keep == "y")
metadata_exp <- filter(metadata_clean,type == "experiment")

dataset_int <- microtable$new(sample_table = metadata_clean, otu_table = asv_table, tax_table = taxonomy_table)
dataset_int_e <- microtable$new(sample_table = metadata_exp, otu_table = asv_table, tax_table = taxonomy_table)
dataset_int
#remove contaminat features
# dataset_int$filter_pollution(taxa = c("mitochondria","uncultured", "chloroplast","metagenome","eukaryote"))
# dataset_int_e$filter_pollution(taxa = c("mitochondria","uncultured", "chloroplast","metagenome","eukaryote"))


# It is better to have a backup before filtering features
dataset_filter <- clone(dataset_int)
dataset_filter_e <- clone(dataset_int_e)

# In this example, mean relative abundance threshold 0.0001
# occurrence frequency 0.1; 10% samples have the target features
# dataset_filter$filter_taxa(rel_abund = 0.0001, freq = 0.01)
# dataset_filter_e$filter_taxa(rel_abund = 0.0001, freq = 0.01)

dataset_filter$sample_sums() %>% range


#make asv table and taxonomy consistent
dataset_filter$tidy_dataset()
dataset_filter_e$tidy_dataset()


dataset_rar<- clone(dataset_filter)
dataset_rar_e<- clone(dataset_filter_e)

# As an example, use 10000 sequences in each sample
# dataset_rar$rarefy_samples(sample.size =4143);dataset_rar
# dataset_rar_e$rarefy_samples(sample.size =4143);dataset_rar_e



dataset_rar$sample_sums() %>% range

# use default parameters
dataset_rar$cal_abund()
dataset_rar_e$cal_abund()
dataset_rar_e$cal_betadiv(unifrac = FALSE)
dataset_rar_e$cal_alphadiv()
# View(dataset_rar_e$taxa_abund$Genus)
# View(dataset_rar_e$taxa_abund$Species)

# tab-delimited, i.e. mpa format
# dataset_rar$save_abund(merge_all = TRUE, sep = "\t", quote = FALSE)

# remove those unclassified
dataset_rar$save_abund(merge_all = TRUE, sep = "\t", rm_un = TRUE, rm_pattern = "__$|Sedis$", quote = FALSE)

# View(dataset_rar$otu_table)
# View(dataset_rar$tax_table)
# test=dataset_rar$taxa_abund$Genus


##Add acid_data to Metadata

acids_env <- acids_env[ , !names(acids_env) %in% c("ISOC6","ISOC5","ISOC4")]
dataset_rar_e$sample_table <- data.frame(dataset_rar_e$sample_table, acids_env[rownames(dataset_rar_e$sample_table), ])
View(dataset_rar_e$sample_table)

## CREATE DATA SETS PER FEEDSTOCK

# View(dataset_rar_e_F$sample_table)

dataset_rar_e_F <- clone(dataset_rar_e)
dataset_rar_e_M <- clone(dataset_rar_e)
dataset_rar_e_F$sample_table <- subset(dataset_rar_e$sample_table,feedstock == "F");dataset_rar_e_F
dataset_rar_e_M$sample_table <- subset(dataset_rar_e$sample_table,feedstock == "M");dataset_rar_e_M

dataset_rar_e_F$tidy_dataset();dataset_rar_e_F
dataset_rar_e_M$tidy_dataset();dataset_rar_e_M


################################################################################
### MAASLIN2
################################################################################
library(Maaslin2)
library(tidyr)
library(dplyr)


# Load the necessary libraries
# 
# pathway_df$taxonomy <- rownames(pathway_df)
# 
# pathway_df_2 <- pathway_df %>%
#   separate(taxonomy, 
#            into = c("Superclass1", "Superclass2", "pathway"), 
#            sep = "\\|", fill = "right")
# 
# pathway_df_2 <- pathway_df_2 %>%
#   mutate(Superclass2_pathway = paste(Superclass2, pathway, sep = "|"))
# 
# pathway_formated <- pathway_df_2 %>%
#   select(starts_with("J")) %>%
#   mutate(Superclass2_pathway = pathway_df_2$Superclass2_pathway)
# 
# pathway_formated$Superclass2_pathway <- gsub(" ", "", pathway_formated$Superclass2_pathway)
# 
# 
# # Step 3: Assign Phylum_Genus as row names
# rownames(pathway_formated) <- pathway_formated$Superclass2_pathway
# 
# # Step 4: Remove the Phylum_Genus column as it's now in the row names
# genus_formated <- genus_formated %>%
#   select(-Family_Genus)
# 
pathway_df=dataset_rar_e$taxa_abund$process_type
pathway_df


## MAASLIN2 FOR EACH FEEDSTOCK



fit_data_fixed_effects_F = Maaslin2(input_data     = pathway_df, 
                                    input_metadata = metadata_0[metadata_0$feedstock == "F", ], 
                                    normalization  = "TSS",
                                    output         = "F_Maaslin2_fixed_effects_process_type_0.05", 
                                    fixed_effects  = c("pH","temperature","inoculum"),
                                    reference      = c("pH,7","temperature,45","inoculum,A"),
                                    random_effects = c("reactor","run"),
                                    min_prevalence = 0,
                                    min_abundance  = 0.05)


fit_data_fixed_effects_M = Maaslin2(input_data     = pathway_df, 
                                    input_metadata = metadata_0[metadata_0$feedstock == "M", ], 
                                    normalization  = "TSS",
                                    output         = "M_Maaslin2_fixed_effects_process_type_0.05", 
                                    fixed_effects  = c("pH","temperature","inoculum"),
                                    reference      = c("pH,7","temperature,45","inoculum,A"),
                                    random_effects = c("reactor","run"),
                                    min_prevalence = 0,
                                    min_abundance  = 0.05)




fit_data_int_effects_F = Maaslin2(input_data     = pathway_df, 
                                    input_metadata = metadata_0[metadata_0$feedstock == "F", ], 
                                    normalization  = "TSS",
                                    output         = "F_Maaslin2_int_effects_process_type", 
                                    fixed_effects  = c("inoculum_temperature_feedstock_pH"),
                                    reference      = c("inoculum_temperature_feedstock_pH,A_35_F_7"),
                                    random_effects = c("reactor","run"),
                                    min_prevalence = 0,
                                    min_abundance  = 0.05)

fit_data_int_effects_M = Maaslin2(input_data     = pathway_df, 
                                  input_metadata = metadata_0[metadata_0$feedstock == "M", ], 
                                  normalization  = "TSS",
                                  output         = "M_Maaslin2_int_effects_process_type", 
                                  fixed_effects  = c("inoculum_temperature_feedstock_pH"),
                                  reference      = c("inoculum_temperature_feedstock_pH,A_35_M_7"),
                                  random_effects = c("reactor","run"),
                                  min_prevalence = 0,
                                  min_abundance  = 0.05)






pathway_df=dataset_rar_e$taxa_abund$pathway


fit_data_fixed_effects_F = Maaslin2(input_data     = pathway_df, 
                                    input_metadata = metadata_0[metadata_0$feedstock == "F", ], 
                                    normalization  = "TSS",
                                    output         = "F_Maaslin2_fixed_effects_pathway_0.05", 
                                    fixed_effects  = c("pH","temperature","inoculum"),
                                    reference      = c("pH,7","temperature,45","inoculum,A"),
                                    random_effects = c("reactor","run"),
                                    min_prevalence = 0,
                                    min_abundance  = 0.05)


fit_data_fixed_effects_M = Maaslin2(input_data     = pathway_df, 
                                    input_metadata = metadata_0[metadata_0$feedstock == "M", ], 
                                    normalization  = "TSS",
                                    output         = "M_Maaslin2_fixed_effects_pathway_0.05", 
                                    fixed_effects  = c("pH","temperature","inoculum"),
                                    reference      = c("pH,7","temperature,45","inoculum,A"),
                                    random_effects = c("reactor","run"),
                                    min_prevalence = 0,
                                    min_abundance  = 0.05)




fit_data_int_effects_F = Maaslin2(input_data     = pathway_df, 
                                  input_metadata = metadata_0[metadata_0$feedstock == "F", ], 
                                  normalization  = "TSS",
                                  output         = "F_Maaslin2_int_effects_pathway", 
                                  fixed_effects  = c("inoculum_temperature_feedstock_pH"),
                                  reference      = c("inoculum_temperature_feedstock_pH,A_35_F_7"),
                                  random_effects = c("reactor","run"),
                                  min_prevalence = 0,
                                  min_abundance  = 0.05)

fit_data_int_effects_M = Maaslin2(input_data     = pathway_df, 
                                  input_metadata = metadata_0[metadata_0$feedstock == "M", ], 
                                  normalization  = "TSS",
                                  output         = "M_Maaslin2_int_effects_pathway", 
                                  fixed_effects  = c("inoculum_temperature_feedstock_pH"),
                                  reference      = c("inoculum_temperature_feedstock_pH,A_35_M_7"),
                                  random_effects = c("reactor","run"),
                                  min_prevalence = 0,
                                  min_abundance  = 0.05)



####

############
### PEARSON CORRELATIONS GENUS
##########

t1 <- trans_env$new(dataset = dataset_rar_e, env_cols = 25:31)

t1$cal_cor(use_data = "process_type", p_adjust_method = "fdr",by_group = "feedstock", p_adjust_type = "Env")
t1$cal_cor(use_data = "pathway", p_adjust_method = "fdr", p_adjust_type = "Env")
t1$cal_cor(use_data = "pathway", p_adjust_method = "fdr",by_group = "feedstock", p_adjust_type = "Env")

View(t1$res_cor)
t1$plot_cor()
t1$plot_cor(filter_feature = c("***"))

png("FigSX_Pearson_feedstock_process_type.png", width = 10, height = 20, units = "in", res = 300)
t1$plot_cor(filter_feature = c("***"))
dev.off()





##install.packages("paletteer")
# examples
paletteer_d("RColorBrewer::Spectral")
paletteer_d("RColorBrewer::Set3")
paletteer_d("RColorBrewer::Paired")
paletteer_d("ggsci::nrc_npg")
paletteer_d("ggsci::default_aaas")
paletteer_d("ggsci::lanonc_lancet")
paletteer_d("ggsci::default_nejm")
paletteer_d("ggsci::category10_d3")

paletteer_d("ggthemes::Classic_10_Light")
paletteer_d("ggthemes::Classic_10_Medium")
paletteer_d("ggthemes::Classic_Cyclic")
# 20 color values

paletteer_d("ggsci::category20c_d3")
paletteer_d("ggsci::category20_d3")

paletteer_d("ggthemes::Classic_20")
RColorBrewer::brewer.pal(15, "Dark2")

color_values = RColorBrewer:ggsci::category10_d3(15, "Dark2")
RColorBrewer::brewer.pal(12, "Dark2")
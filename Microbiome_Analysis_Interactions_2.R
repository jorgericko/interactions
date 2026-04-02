#Microbiome analysis interactions study, let's do it!

setwd("~/Projects/BETO/Microbial Analysis/2024_Interactions/16S")


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

metadata_0<-read.table("metadata_interactions.txt", sep='\t', header=T, row.names=1, comment="")
# taxonomy_table<-read.table("taxonomy_table_2.txt", sep='\t', header=T, row.names=1, comment="")
taxonomy_table<-read.table("taxonomy_table_3.txt", sep='\t', header=T, row.names=1, comment="") ## editing unclass genus

acids_env<-read.csv("data_acids_with_sid.csv", header=T, row.names=1, comment="") #Combined data
asv_table<-read.table("feature_table.txt", sep='\t', header=T, row.names=1, comment="")


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
phylo_file_path <- system.file("extdata", "rep_phylo.tre", package="file2meco")
tree <- ape::read.tree("rooted-tree/tree.nwk")

#create microeco environment

metadata_clean <- filter(metadata_0,keep == "y")
metadata_exp <- filter(metadata_clean,type == "experiment")


dataset_int <- microtable$new(sample_table = metadata_clean, otu_table = asv_table, tax_table = taxonomy_table,phylo_tree = tree)
dataset_int_e <- microtable$new(sample_table = metadata_exp, otu_table = asv_table, tax_table = taxonomy_table,phylo_tree = tree)


#remove contaminat features
dataset_int$filter_pollution(taxa = c("mitochondria","uncultured", "chloroplast","metagenome","eukaryote"))
dataset_int_e$filter_pollution(taxa = c("mitochondria","uncultured", "chloroplast","metagenome","eukaryote"))


# It is better to have a backup before filtering features
dataset_filter <- clone(dataset_int)
dataset_filter_e <- clone(dataset_int_e)

# In this example, mean relative abundance threshold 0.0001
# occurrence frequency 0.1; 10% samples have the target features
dataset_filter$filter_taxa(rel_abund = 0.0001, freq = 0.01)
dataset_filter_e$filter_taxa(rel_abund = 0.0001, freq = 0.01)

dataset_filter$sample_sums() %>% range

#make asv table and taxonomy consistent
dataset_filter$tidy_dataset()
dataset_filter_e$tidy_dataset()

View(dataset_filter$sample_table)

library(microeco)
library(mecodev)

# set.seed is used to fix the random number generation to make the results repeatable
set.seed(123)
# trans_rarefy class
t1 <- trans_rarefy$new(dataset_filter, alphadiv = "Observed", depth = c(0,10,100,600, 2000, 4143, 5000,10000,20000))
t1$plot_rarefy(color = "feedstock", show_point = TRUE, add_fitting = FALSE,point_size = 2)
t2 <- trans_rarefy$new(dataset_filter, alphadiv = "Coverage", depth = c(0,10,100,600, 2000, 4143, 5000,10000,20000))

View(t1)

png("Figures_Manuscript/FigSX_rarefaction_curves.png", width = 6, height = 6, units = "in", res = 300)
t1$plot_rarefy(color = "feedstock", show_point = TRUE, add_fitting = FALSE,point_size = 2)
dev.off()

png("Figures_Manuscript/FigSX_rarefaction_curves_coverage.png", width = 6, height = 6, units = "in", res = 300)
t2$plot_rarefy(color = "feedstock", show_point = TRUE, add_fitting = FALSE,point_size = 2)
dev.off()

dataset_rar<- clone(dataset_filter)
dataset_rar_e<- clone(dataset_filter_e)


# Calculate coverage at your chosen depth (4143)
t1_coverage <- trans_rarefy$new(dataset_filter, alphadiv = "Coverage", depth = 4143)
View(t1_coverage$res_rarefy)

# View the mean coverage across all samples
mean_coverage <- mean(t1_coverage$res_rarefy$mean)
mean_coverage

# To see the range (min and max) of coverage across samples
range(t1_coverage$res_rarefy$mean)



# Calculate mean coverage (convert to numeric)
mean_coverage <- mean(as.numeric(t1_coverage$res_rarefy$Coverage))
mean_coverage * 100  # to see as percentage

# Range
min(as.numeric(t1_coverage$res_rarefy$Coverage)) * 100
max(as.numeric(t1_coverage$res_rarefy$Coverage)) * 100



# As an example, use 10000 sequences in each sample
dataset_rar$rarefy_samples(sample.size =4143);dataset_rar
dataset_rar_e$rarefy_samples(sample.size =4143);dataset_rar_e



dataset_rar$sample_sums() %>% range

# use default parameters
dataset_rar$cal_abund()
dataset_rar_e$cal_abund()
dataset_rar_e$cal_betadiv(unifrac = TRUE)
dataset_rar_e$cal_alphadiv()
# View(dataset_rar_e$taxa_abund$Genus)
# View(dataset_rar_e$taxa_abund$Species)

# tab-delimited, i.e. mpa format
dataset_rar$save_abund(merge_all = TRUE, sep = "\t", quote = FALSE)

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

genus_df=dataset_rar_e$taxa_abund$Genus

# Load the necessary libraries

genus_df_2$taxonomy <- rownames(genus_df)

genus_df_2 <- genus_df_2 %>%
  separate(taxonomy, 
           into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus"), 
           sep = "\\|", fill = "right")

genus_df_2 <- genus_df_2 %>%
  mutate(Family_Genus = paste(Family, Genus, sep = "|"))

genus_formated <- genus_df_2 %>%
  select(starts_with("J")) %>%
  mutate(Family_Genus = genus_df_2$Family_Genus)

genus_formated$Family_Genus <- gsub(" ", "", genus_formated$Family_Genus)


# Step 3: Assign Phylum_Genus as row names
rownames(genus_formated) <- genus_formated$Family_Genus

# Step 4: Remove the Phylum_Genus column as it's now in the row names
genus_formated <- genus_formated %>%
  select(-Family_Genus)


fit_data_fixed_effects_all = Maaslin2(input_data     = genus_formated, 
                                    input_metadata = metadata_all, 
                                    normalization  = "TSS",
                                    output         = "Maaslin2_fixed_effects_genus_trt_0.05", 
                                    fixed_effects  = c("feedstock","pH","temperature","inoculum"),
                                    reference      = c("feedstock,M","pH,7","temperature,35","inoculum,A"),
                                    random_effects = c("reactor","run"),
                                    min_prevalence = 0,
                                    min_abundance  = 0.05)


fit_data_int_effects_all = Maaslin2(input_data     = genus_formated, 
                                      input_metadata = metadata_all, 
                                      normalization  = "TSS",
                                      output         = "Maaslin2_int_effects_genus_trt_0.05", 
                                      fixed_effects  = c("inoculum_temperature_feedstock_pH"),
                                      reference      = c("inoculum_temperature_feedstock_pH,A_35_M_7"),
                                      random_effects = c("reactor","run"),
                                      min_prevalence = 0,
                                      min_abundance  = 0.05)




genus_df_F=dataset_rar_e_F$taxa_abund$Genus
genus_df_M=dataset_rar_e_M$taxa_abund$Genus


## MAASLIN2 FOR EACH FEEDSTOCK FOR NETWORKS

fit_data_fixed_effects_F = Maaslin2(input_data     = genus_df_F, 
                           input_metadata = metadata_F, 
                           normalization  = "TSS",
                           output         = "F_Maasli2_fixed_effects_genus_trt_0.0001", 
                           fixed_effects  = c("pH","temperature","inoculum"),
                           reference      = c("pH,7","temperature,45","inoculum,A"),
                           random_effects = c("reactor","run"),
                           min_prevalence = 0,
                           min_abundance  = 0.0001)


fit_data_fixed_effects_M = Maaslin2(input_data     = genus_df_M, 
                                    input_metadata = metadata_M, 
                                    normalization  = "TSS",
                                    output         = "M_Maasli2_fixed_effects_genus_trt_0.0001", 
                                    fixed_effects  = c("pH","temperature","inoculum"),
                                    reference      = c("pH,7","temperature,45","inoculum,A"),
                                    random_effects = c("reactor","run"),
                                    min_prevalence = 0,
                                    min_abundance  = 0.0001)

# head(maaslin2_F)
maaslin2_F=fit_data_fixed_effects_F$results
maaslin2_M=fit_data_fixed_effects_M$results


library(dplyr)

# Assuming your dataframe is called df
maaslin2_F <- maaslin2_F %>%
  mutate(feature = gsub("\\.\\.", "|", feature))

maaslin2_M <- maaslin2_M %>%
  mutate(feature = gsub("\\.\\.", "|", feature))



# Load necessary libraries
library(dplyr)
library(tidyr)
library(stringr)


maaslin2_F <- maaslin2_F %>%
  mutate(genus = str_extract(feature, "[^|]+$"))%>% 
  mutate(genus = gsub("\\.\\.", "|", genus))


maaslin2_M <- maaslin2_M %>%
  mutate(genus = str_extract(feature, "[^|]+$"))



maaslin2_F <- maaslin2_F %>%
  mutate(genus = str_replace_all(genus, c("\\.g" = "|g", "\\.f" = "|f", "\\.o" = "|o", "\\.c" = "|c")))



maaslin2_M <- maaslin2_M %>%
  mutate(genus = str_replace_all(genus, c("\\.g" = "|g", "\\.f" = "|f", "\\.o" = "|o", "\\.c" = "|c")))



# Step 2: Annotate significance based on p-value
maaslin2_F <- maaslin2_F %>%
  mutate(significance = case_when(
    pval < 0.001 ~ "***",
    pval < 0.01 ~ "**",
    pval < 0.05 ~ "*",
    pval < 0.1 ~ ".",
    TRUE ~ "ns"
  ))

maaslin2_M <- maaslin2_M %>%
  mutate(significance = case_when(
    pval < 0.001 ~ "***",
    pval < 0.01 ~ "**",
    pval < 0.05 ~ "*",
    pval < 0.1 ~ ".",
    TRUE ~ "ns"
  ))


maaslin2_F <- maaslin2_F %>%
  mutate(metadata_value = paste(metadata, value, sep = "_"))

maaslin2_M <- maaslin2_M %>%
  mutate(metadata_value = paste(metadata, value, sep = "_"))


maaslin2_F_wide <- maaslin2_F %>%
  select(genus, metadata_value, coef, pval, significance) %>%
  pivot_wider(names_from = metadata_value, 
              values_from = c(coef, pval, significance),
              values_fn = list) %>%
  unnest(cols = everything())  # Flatten list columns to regular numeric columns


maaslin2_M_wide <- maaslin2_M %>%
  select(genus, metadata_value, coef, pval, significance) %>%
  pivot_wider(names_from = metadata_value, 
              values_from = c(coef, pval, significance),
              values_fn = list) %>%
  unnest(cols = everything())  # Flatten list columns to regular numeric columns



write.csv(maaslin2_M_wide, "maaslin2_M_wide.csv", row.names = FALSE)

write.csv(maaslin2_F_wide, "maaslin2_F_wide.csv", row.names = FALSE)



################################################################################
##### VENN DIAGRAMS
#################################################################################


t1_venn_F <- trans_venn$new(dataset_rar$merge_samples("feedstock"), ratio = "numratio")

t1_venn_I$plot_venn()

t1_venn_I <- trans_venn$new(dataset_rar$merge_samples("inoculum"), ratio = "numratio")
# 
# t1_venn_F$plot_venn()
# t1_venn_I$plot_venn()

grid.arrange(t1_venn_F$plot_venn(),t1_venn_I$plot_venn(), ncol = 2)

# Set the filename and dimensions (in inches)
png("Figures_Manuscript/Fig6_Venn_diagrams.png", width = 10, height = 4, units = "in", res = 300)
grid.arrange(t1_venn_F$plot_venn(),t1_venn_I$plot_venn(), ncol = 2)
dev.off()



################################################################################
#### BETA DIVERSITY - PCoA
################################################################################

library(microeco)
library(magrittr)
library(ggplot2)
library(aplot)
library(ggpubr)
theme_set(theme_pubr())


# create trans_beta object
dataset_rar$cal_betadiv(unifrac = TRUE)
dataset_rar_e$cal_betadiv(unifrac = TRUE)


# measure parameter should be either one of names(dataset$beta_diversity) or a customized symmetric matrix
t1 <- trans_beta$new(dataset = dataset_rar, group = "feedstock_pH", measure = "unwei_unifrac")
# PCoA, PCA, DCA and NMDS are available
t1$cal_ordination(method = "PCoA")
# t1$res_ordination is the ordination result list
class(t1$res_ordination)
# save dataset$beta_diversity to a directory
dataset_rar$save_betadiv(dirpath = "beta_diversity")
# extract the axis scores
tmp <- t1$res_ordination$scores
# differential test with trans_env class
t2 <- trans_env$new(dataset = dataset_rar, add_data = tmp[, 1:2])
# 'KW_dunn' for non-parametric test
t2$cal_diff(group = "feedstock_pH", method = "anova")

colors_PCoA = c("purple","#CD534CFF","#EFC000FF", "#0073C2FF","darkgreen","#A73030FF","#8F7700FF","#003C67FF")

p1 <- t1$plot_ordination(plot_color = "feedstock_pH",
                         color_values = colors_PCoA,
                         shape_values = c(16, 15, 8, 7, 18, 17),
                         plot_shape = "inoculum_temperature",
                         plot_type = c("point", "ellipse"))
p1

# groups order in p2 is same with p1; use legend.position = "none" to remove redundant legend
p2 <- t2$plot_diff(measure = "PCo1", add_sig = T,color_values = colors_PCoA) + theme_pubr() + coord_flip() + 
  theme(legend.position = "none", axis.title.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())

p3 <- t2$plot_diff(measure = "PCo2", add_sig = T,color_values = colors_PCoA) + theme_pubr() + 
  theme(legend.position = "none", axis.title.y = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())

# height of the upper figure and width of the right-hand figure are both 0.2-fold of the main figure
PCoA_all <- p1 %>% insert_top(p2, height = 0.2) %>% insert_right(p3, width = 0.2)

# Set the filename and dimensions (in inches)
png("Figures_Manuscript/Fig5_PCoA_all.png", width = 9, height = 9, units = "in", res = 300)
PCoA_all
dev.off()






###########
### FATTY ACID DATA PLOTS
#################

library(ggplot2)


t1 <- trans_env$new(dataset = dataset_rar_e, add_data = acids_env)
t1 <- trans_env$new(dataset = dataset_rar_e, env_cols = 25:31)

t1$cal_diff(method = "anova", group = "inoculum_temperature_feedstock_pH")

color_palette<-c("red",     "gold", "blue", 
                 "#A50F15", "yellow4", "darkblue", 
                 "#FCAE91" ,"khaki2", "#BDD7E7", 
                 "#FB6A4A", "yellow","#6BAED6",
                 "red",     "gold", "blue", 
                 "#A50F15", "yellow4", "darkblue", 
                 "#FCAE91" ,"khaki2", "#BDD7E7", 
                 "#FB6A4A", "yellow","#6BAED6")

# place all the plots into a list
tmp <- list()
for(i in colnames(t1$data_env)){
  tmp[[i]] <- t1$plot_diff(measure = i, add_sig_text_size = 0.1, xtext_size = 12, color_values=color_palette) + theme(plot.margin = unit(c(0.1, 0, 0, 1), "cm"))
}

plot(gridExtra::arrangeGrob(grobs = tmp, ncol = 2))


### CORRELATION BETWEEN FA DATA
t1$cal_autocor(group = "feedstock")
t1$cal_autocor()


####

###########
### RDA AND DBRDA
#################



color_palette<-c("darkred", "yellow4", "darkblue", 
                 "firebrick1", "darkgoldenrod1","#6BAED6")

t_F <- trans_env$new(dataset = dataset_rar_e_F, add_data = acids_env)
t_M <- trans_env$new(dataset = dataset_rar_e_M, add_data = acids_env)


#### dbRDA PLOTS

t_F$cal_ordination(method = "dbRDA", use_measure = "unwei_unifrac")
t_M$cal_ordination(method = "dbRDA", use_measure = "unwei_unifrac")

t_F$trans_ordination(adjust_arrow_length = TRUE, max_perc_env = 1.5)
t_M$trans_ordination(adjust_arrow_length = TRUE, max_perc_env = 1.5)

t_F_plot=t_F$plot_ordination(plot_color = "inoculum_pH", plot_shape = "temperature", plot_type = c("point", "ellipse"),color_values = color_palette,
                             ellipse_level = .95,
                             ellipse_chull_alpha = 0.04,
                             point_alpha = 0.6,
                             point_size =1.5,
                             env_text_size = 4.5, # Slightly smaller for clarity
                             env_text_color = "darkred",
                             env_arrow_color = "darkred")+ 
  theme(
    legend.position = "none",
    axis.title = element_text(size = 15), # Bolds axis titles
    axis.text = element_text(size = 15),
    # This bolds the actual taxa/env labels via ggplot override if the function didn't
  )

t_M_plot=t_M$plot_ordination(plot_color = "inoculum_pH", plot_shape = "temperature", plot_type = c("point", "ellipse"),color_values = color_palette,
                             ellipse_level = .95,
                             ellipse_chull_alpha = 0.04,
                             point_alpha = 0.6,
                             point_size =1.5,
                             env_text_size = 4.5, # Slightly smaller for clarity
                             env_text_color = "darkred",
                             env_arrow_color = "darkred")+ 
  theme(
    legend.position = "none",
    axis.title = element_text(size = 15), # Bolds axis titles
    axis.text = element_text(size = 15),
    # This bolds the actual taxa/env labels via ggplot override if the function didn't
  )


# t_F$cal_ordination_anova()
# t_F$cal_ordination_envfit()
# t_M$cal_ordination_anova()
# t_M$cal_ordination_envfit()
# 
# #PRINT RESULT TABLES COPY PASTE INTO SUPPLEMENTAL INFO dbRDA
# 
# t_F$res_ordination_terms
# t_F$res_ordination_envfit
# t_F$res_ordination_R2
# t_M$res_ordination_terms
# t_M$res_ordination_envfit
# t_M$res_ordination_R2




#RDA
t_F$cal_ordination(method = "RDA", taxa_level = "Genus")
t_M$cal_ordination(method = "RDA", taxa_level = "Genus")


t_F$trans_ordination(show_taxa = 10, adjust_arrow_length = TRUE, max_perc_env = 1.5, max_perc_tax = 1.5, min_perc_env = 0.2, min_perc_tax = 0.2)
t_M$trans_ordination(show_taxa = 10, adjust_arrow_length = TRUE, max_perc_env = 1.5, max_perc_tax = 1.5, min_perc_env = 0.2, min_perc_tax = 0.2)

#t_F_plot_tax=t_F$plot_ordination(plot_color = "inoculum_pH", plot_shape = "temperature",color_values = color_palette)

#t_M_plot_tax=t_M$plot_ordination(plot_color = "inoculum_pH", plot_shape = "temperature",color_values = color_palette)


t_F_plot_tax=t_F$plot_ordination(
  plot_color = "inoculum_pH", 
  plot_shape = "temperature",
  plot_type = c("point"),
  ellipse_level = .95,
  ellipse_chull_alpha = 0.05,
  color_values = color_palette,
  point_alpha = 0.6,
  point_size = 1.5,
  taxa_text_size = 3,
  # Change: Use taxa_text_style for bold/italic formatting
  env_text_size = 4.5,
  taxa_text_color = "black",
  taxa_arrow_color = "grey35",
  env_text_color = "darkred",
  env_arrow_color = "darkred"
) + 
  theme(
    legend.position = "none",
    axis.title = element_text(size = 15), # Bolds axis titles
    axis.text = element_text(size = 15),
    # This bolds the actual taxa/env labels via ggplot override if the function didn't
  )


t_M_plot_tax = t_M$plot_ordination(
  plot_color = "inoculum_pH", 
  plot_shape = "temperature",
  plot_type = c("point"),
  ellipse_level = .95,
  ellipse_chull_alpha = 0.05,
  color_values = color_palette,
  point_alpha = 0.6,
  point_size = 1.5,
  taxa_text_size = 3,
  # Change: Use taxa_text_style for bold/italic formatting
  env_text_size = 4.5,
  taxa_text_color = "black",
  taxa_arrow_color = "grey35",
  env_text_color = "darkred",
  env_arrow_color = "darkred"
) + 
  theme(
    legend.position = "none",
    axis.title = element_text(size = 15), # Bolds axis titles
    axis.text = element_text(size = 15),
  )



# 1. Create the full plot object
t_M_plot_legend = t_M$plot_ordination(
  plot_color = "inoculum_pH", 
  plot_shape = "temperature", 
  plot_type = c("point", "ellipse"),
  color_values = color_palette,
  ellipse_level = .95,
  ellipse_chull_alpha = 0.04,
  point_alpha = 0.6,
  point_size = 1.5,
  env_text_size = 4.5,
  env_text_color = "darkred",
  env_arrow_color = "darkred"
) +  
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 12, face = "bold"),   # Bolds legend text
    legend.title = element_text(size = 13, face = "bold")  # Bolds legend title
  )


t_M_plot_legend = t_M_plot_legend + 
  theme(
    # 1. Position the legend (Change from "none" to see it!)
    legend.position = "bottom",
    legend.box = "vertical",           # Stack color and shape legends
    
    # 2. Bold and increase the Title (e.g., "inoculum_pH")
    legend.title = element_text(size = 14, face = "bold"), 
    
    # 3. Bold and increase the Labels (e.g., "pH 5", "pH 7")
    legend.text = element_text(size = 12, face = "bold"),
    
    # 4. Increase the size of the colored symbols/icons themselves
    legend.key.size = unit(1.2, "cm"), 
    
    # Keep your axis settings
    axis.title = element_text(size = 15, face = "bold"),
    axis.text = element_text(size = 15, face = "bold")
  ) +
  # This helper ensures the points inside the legend are larger/clearer
  guides(
    color = guide_legend(override.aes = list(size = 5)),
    shape = guide_legend(override.aes = list(size = 6))
  )


t_M_plot_legend


# dbRDA and RDA plots
library(gridExtra)
grid.arrange(t_F_plot,t_F_plot_tax, t_M_plot,t_M_plot_tax, ncol = 2)

grid.arrange(t_F_plot_tax,t_M_plot_tax, ncol = 2)



t_M_plot_legend

png("Figures_Manuscript/Fig5_dBRDA_RDA_diagrams_top10.png", width = 9, height = 8, units = "in", res = 300)
grid.arrange(t_F_plot,t_F_plot_tax, t_M_plot,t_M_plot_tax, ncol = 2)
dev.off()

png("Figures_Manuscript/Fig5_dBRDA_RDA_diagrams_top10_legend.png", width = 9, height = 8, units = "in", res = 300)
t_M_plot_legend
dev.off()




#RDA analysis

t_F$cal_ordination_anova()
t_F$cal_ordination_envfit()
t_M$cal_ordination_anova()
t_M$cal_ordination_envfit()

#PRINT RESULT TABLES COPY PASTE INTO SUPPLEMENTAL INFO dbRDA

t_F$res_ordination_terms
t_F$res_ordination_envfit
t_F$res_ordination_R2
t_M$res_ordination_terms
t_M$res_ordination_envfit
t_M$res_ordination_R2
t_F_plot_tax
t_M_plot
t_M_plot_tax

############
### PEARSON CORRELATIONS GENUS
##########

t1 <- trans_env$new(dataset = dataset_rar_e, env_cols = 25:31)

t1$cal_cor(use_data = "Genus", p_adjust_method = "fdr",by_group = "feedstock", p_adjust_type = "Env")
View(t1$res_cor)
t1$plot_cor()
t1$plot_cor(filter_feature = c("", "*", "**"))

png("Figures_Manuscript/FigSX_Pearson_genus.png", width = 7, height = 11, units = "in", res = 300)
t1$plot_cor(filter_feature = c("", "*", "**"))
dev.off()


############
### PEARSON ALPHA DIVERSITY
##########

tf <- trans_env$new(dataset = dataset_rar_e_F, add_data = acids_env)
tm <- trans_env$new(dataset = dataset_rar_e_M, add_data = acids_env)


tf$cal_cor(add_abund_table = dataset_rar_e_F$alpha_diversity)
tm$cal_cor(add_abund_table = dataset_rar_e_M$alpha_diversity)
# try to use ggplot2 with clustering plot
# require ggtree and aplot packages to be installed (https://chiliubio.github.io/microeco_tutorial/intro.html#dependence)

png("Figures_Manuscript/FigSX_alpha_correlation_f.png", width = 3, height = 5, units = "in", res = 300)
tf$plot_cor(cluster_ggplot = "both")
dev.off()

png("Figures_Manuscript/FigSX_alpha_correlation_m.png", width = 3, height = 5, units = "in", res = 300)
tm$plot_cor(cluster_ggplot = "both")
dev.off()

tf$plot_cor(cluster_ggplot = "both")


################################################################################
#### BETA DIVERSITY - PERMANOVA - TABLE 1 - TABLE S
################################################################################

#PERMANOVA WITH INTERACTIONS EFFECTS
t3 <- trans_beta$new(dataset = dataset_rar_e, measure = "unwei_unifrac")
t3$cal_manova(manova_set = "feedstock*pH*temperature*inoculum") 
t3$res_manova #PRINT AND COPY INTO EXCEL TABLE 1

color_palette<-c("red",     "gold", "blue", 
                 "#A50F15", "yellow4", "darkblue", 
                 "#FCAE91" ,"khaki2", "#BDD7E7", 
                 "#FB6A4A", "yellow","#6BAED6")

t3$plot_clustering(group = "inoculum_feedstock_pH", replace_name = c("inoculum_temperature_feedstock_pH"),color_values = color_palette)

# inoculum_temperature_feedstock_pH


# the group parameter is not necessary when it is provided in creating the object


#PERMANOVA GROUP COMPARISON 
t4 <- trans_beta$new(dataset = dataset_rar_e, group = "feedstock_pH", measure = "unwei_unifrac")
t4$cal_group_distance(within_group = FALSE)
t4$cal_group_distance_diff(method = "anova")

t4$plot_group_distance()

#PERMANOVA GROUP COMPARISON pH feedstock
t4 <- trans_beta$new(dataset = dataset_rar_e, group = "feedstock_pH", measure = "unwei_unifrac")
t4$cal_manova(manova_all = FALSE, group = "pH", by_group = "feedstock")
t4$res_manova #PRINT AND COPY INTO EXCEL TABLE S5


# parameters in plot_group_distance function will be passed to the plot_alpha function of trans_alpha class
t4$plot_group_distance(plot_type = "ggviolin", add = "mean_se")
t4$plot_group_distance(add = "mean")


# DISTANCE PLOT FOR FEEDSTOCK PH INTERACTION
distance_plot <- t4$plot_group_distance(jitter = FALSE)  # Check if 'jitter' argument is available


t4 <- trans_beta$new(dataset = dataset_rar_e, group = "inoculum_temperature_pH", measure = "unwei_unifrac")
t4$cal_manova(manova_all = FALSE, group = "inoculum_temperature_pH", by_group = "feedstock")
permanova_it_fpH=t4$res_manova ### SAVE AS DATA FRAME

t4$res_manova #PRINT AND COPY INTO EXCEL TABLE S5




t4 <- trans_beta$new(dataset = dataset_rar_e, group = "inoculum_temperature_feedstock_pH", measure = "unwei_unifrac")
t4$cal_manova(manova_all = FALSE, group = "inoculum", by_group = "temperature_feedstock_pH")
# permanova_it_fpH=t4$res_manova ### SAVE AS DATA FRAME

t4$res_manova #PRINT AND COPY INTO EXCEL TABLE S5


################################################################################
#### TAXONOMIC ANALYSIS
################################################################################
# Inspect the sample_table to understand its structure
str(dataset_rar$sample_table)

dataset_rar$sample_table

t_bp <- trans_abund$new(dataset = dataset_rar, taxrank = "Genus", ntaxa = 24)
t_bp_sp <- trans_abund$new(dataset = dataset_rar, taxrank = "Species", ntaxa = 24)

view(t_bp_sp$data_abund)

# Order the samples by 'day.2', using row.names as the sample identifier
sample_order <- row.names(dataset_rar$sample_table)[order(dataset_rar$sample_table$day)]

# Now use this ordered list of sample names in your plot function
taxa_bp_gen = t_bp$plot_bar(
  others_color = "grey90", 
  color_values = colors_24,
  facet = c("feedstock", "pH", "temperature", "inoculum"), 
  order_x = sample_order,  # Use ordered row names
  xtext_keep = FALSE, 
  legend_text_italic = TRUE, 
  barwidth = 1)


taxa_bp_gen

str(dataset_rar)
str(t_bp)


# t_bp_fam <- trans_abund$new(dataset = dataset_rar, taxrank = "Family", ntaxa = 24)
# 
# taxa_bp_p=t_bp_fam$plot_bar(
#   others_color = "grey90", 
#   color_values =colors_24,
#   facet = c("feedstock", "pH", "temperature", "inoculum"), 
#   order_x = sample_order,  # Use the ordered sample names
#   xtext_keep = FALSE, 
#   legend_text_italic = TRUE, 
#   barwidth = 1)
# 
# taxa_bp_p


png("Figures_Manuscript/Fig6_Genus_diagrams.png", width = 13, height = 6, units = "in", res = 300)
taxa_bp_gen
dev.off()

t1 <- trans_diff$new(dataset = dataset_rar, method = "lefse", group = "inoculum_feedstock_pH", alpha = 0.05, lefse_subgroup = NULL)

t1$plot_diff_cladogram(use_taxa_num = 200, use_feature_num = 50, clade_label_level = 3)








################################################################################
#### TAXONOMIC ANALYSIS - ANCOM-BC
################################################################################

# ANCOMBC2 method
# Assuming the 'pH' column exists and it's a factor, relevel it to set pH 7 as reference
# Check the current levels of the pH variable in the sample_table
levels(factor(dataset_rar_e$sample_table$pH))
levels(factor(dataset_rar_e$sample_table$temperature))
levels(factor(dataset_rar_e$sample_table$feedstock))
levels(factor(dataset_rar_e$sample_table$inoculum))

# Relevel pH to set pH 7 as the reference
dataset_rar_e$sample_table$pH <- relevel(factor(dataset_rar_e$sample_table$pH), ref = "7")
dataset_rar_e$sample_table$temperature <- relevel(factor(dataset_rar_e$sample_table$temperature), ref = "45")
dataset_rar_e$sample_table$feedstock <- relevel(factor(dataset_rar_e$sample_table$feedstock), ref = "M")
dataset_rar_e$sample_table$inoculum <- relevel(factor(dataset_rar_e$sample_table$inoculum), ref = "B")

dataset_rar_e$sample_table

# Proceed with your analysis using trans_diff
t_ABC_fe <- trans_diff$new(dataset = dataset_rar_e, 
                     method = "ancombc2", 
                     fix_formula = "inoculum + feedstock + pH + temperature",
                     rand_formula = "(1 | run)+(1|reactor)",
                     taxa_level = "Genus", 
                     filter_thres = 0.005)



dataset_rar_e_F <- clone(dataset_rar_e)
dataset_rar_e_M <- clone(dataset_rar_e)
dataset_rar_e_F$sample_table <- subset(dataset_rar_e$sample_table,feedstock == "F")
dataset_rar_e_M$sample_table <- subset(dataset_rar_e$sample_table,feedstock == "M")

dataset_rar_e_F$tidy_dataset()
dataset_rar_e_M$tidy_dataset()

View(dataset_rar_e_F$sample_table)

View(dataset_rar_e_F$sample_table)


t_ABC_fe_F <- trans_diff$new(dataset = dataset_rar_e_F, 
                           method = "ancombc2", 
                           fix_formula = "inoculum_temperature_pH",
                           rand_formula = "(1|reactor)",
                           taxa_level = "Genus", 
                           filter_thres = 0.005)

View(t_ABC_fe_F$res_diff_raw)
View(t_ABC_fe_F$res_diff)





# view(t_ABC_fe$res_diff)

# View the structure of the res_diff data
View(t_ABC_fe$res_diff)

# Filter out the intercept from res_diff
filtered_diff <- t_ABC_fe$res_diff[!grepl("(Intercept)", t_ABC_fe$res_diff$Factors), ]

# Use the filtered data to generate the plot without intercept
t_ABC_fe$res_diff <- filtered_diff

# Plot the filtered differential abundance data
t_ABC_fe$plot_diff_bar(keep_full_name = FALSE, 
                       heatmap_cell = "lfc", 
                       heatmap_sig = "Significance", 
                       heatmap_x = "Factors", 
                       heatmap_y = "Taxa",
                       heatmap_lab_fill = "lfc")




t_ABC_F <- trans_diff$new(dataset = dataset_rar_e_F, 
                          method = "ancombc2", 
                          fix_formula = "inoculum+temperature+pH",
                          rand_formula = "(1 | run)+(1|reactor)",
                          taxa_level = "Genus", 
                          filter_thres = 0.005)


t_ABC_F$plot_diff_bar(keep_full_name = FALSE, 
                      heatmap_cell = "lfc", 
                      heatmap_sig = "Significance", 
                      heatmap_x = "Factors", 
                      heatmap_y = "Taxa",
                      heatmap_lab_fill = "lfc")





tl <- trans_diff$new(dataset = dataset_rar_e_F, method = "lefse", group = "pH", alpha = 0.01, lefse_subgroup = NULL)
# see t1$res_diff for the result
# From v0.8.0, threshold is used for the LDA score selection.
tl$plot_diff_bar(threshold = 4)
# we show 20 taxa with the highest LDA (log10)
tl$plot_diff_bar(use_number = 1:30, width = 0.3)


View(tl$res_diff)



t_anova_F <- trans_diff$new(dataset = dataset_rar_e_F, method = "anova", group = "inoculum_temperature_pH", taxa_level = "Genus", filter_thres = 0.001)
t_anova_M <- trans_diff$new(dataset = dataset_rar_e_M, method = "anova", group = "inoculum_temperature_pH", taxa_level = "Genus", filter_thres = 0.001)
View(t_anova_F$res_diff)
View(t_anova_F$res_abund)

color_palette<-c("red",     "gold", "blue", 
                 "#A50F15", "yellow4", "darkblue", 
                 "#FCAE91" ,"khaki2", "#BDD7E7", 
                 "#FB6A4A", "yellow","#6BAED6")



anova_F_plot_1_20 = t_anova_F$plot_diff_abund(
  use_number = 1:20, 
  add_sig = TRUE, 
  coord_flip = TRUE,
  color_values = color_palette,
  barwidth = 0.8,
  text_y_size = 25,
  group_order = c("A_35_5", "A_35_7", "A_35_9",
                  "A_45_5", "A_45_7", "A_45_9",
                  "B_35_5", "B_35_7", "B_35_9",
                  "B_45_5", "B_45_7", "B_45_9")) +
  labs(
    x = "Genus",   # Customize your x-axis label
    y = "Relative Abundance", 
    ) +
  theme_pubr() +  # Apply the publication-ready theme
  theme(
    # title = element_text(size = 30),
    axis.text.x = element_text(size = 30),   # Increase x-axis text size
    axis.title.x = element_text(size = 30),  # Increase x-axis title size
    axis.title.y = element_text(size = 35),  # Increase x-axis title size
    axis.text.y = element_text(size = 30), 
    legend.text = element_text(size = 30),   # Increase legend text size
    legend.title = element_text(size = 35),
    legend.position = "right",# Increase legend title size
    panel.grid.major.x = element_line(color = "gray90", size = 0.5),
    panel.grid.minor.x = element_line(color = "gray90", size = 0.5),# Add gray major grid lines along x-axis
  )


anova_F_plot_21_40 = t_anova_F$plot_diff_abund(
  use_number = 21:40, 
  add_sig = TRUE, 
  coord_flip = TRUE,
  color_values = color_palette,
  barwidth = 0.8,
  text_y_size = 25,
  group_order = c("A_35_5", "A_35_7", "A_35_9",
                  "A_45_5", "A_45_7", "A_45_9",
                  "B_35_5", "B_35_7", "B_35_9",
                  "B_45_5", "B_45_7", "B_45_9")) +
  labs(
    x = "Genus",   # Customize your x-axis label
    y = "Relative Abundance", 
  ) +
  theme_pubr() +  # Apply the publication-ready theme
  theme(
    # title = element_text(size = 30),
    axis.text.x = element_text(size = 30),   # Increase x-axis text size
    axis.title.x = element_text(size = 30),  # Increase x-axis title size
    axis.title.y = element_text(size = 35),  # Increase x-axis title size
    axis.text.y = element_text(size = 30), 
    legend.text = element_text(size = 30),   # Increase legend text size
    legend.title = element_text(size = 35),
    legend.position = "right",# Increase legend title size
    panel.grid.major.x = element_line(color = "gray90", size = 0.5),
    panel.grid.minor.x = element_line(color = "gray90", size = 0.5),# Add gray major grid lines along x-axis
  )


anova_M_plot_1_20 = t_anova_M$plot_diff_abund(
  use_number = 1:20, 
  add_sig = TRUE, 
  coord_flip = TRUE,
  color_values = color_palette,
  barwidth = 0.8,
  text_y_size = 25,
  group_order = c("A_35_5", "A_35_7", "A_35_9",
                  "A_45_5", "A_45_7", "A_45_9",
                  "B_35_5", "B_35_7", "B_35_9",
                  "B_45_5", "B_45_7", "B_45_9")) +
  labs(
    x = "Genus",   # Customize your x-axis label
    y = "Relative Abundance", 
  ) +
  theme_pubr() +  # Apply the publication-ready theme
  theme(
    # title = element_text(size = 30),
    axis.text.x = element_text(size = 30),   # Increase x-axis text size
    axis.title.x = element_text(size = 30),  # Increase x-axis title size
    axis.title.y = element_text(size = 35),  # Increase x-axis title size
    axis.text.y = element_text(size = 30), 
    legend.text = element_text(size = 30),   # Increase legend text size
    legend.title = element_text(size = 35),
    legend.position = "right",# Increase legend title size
    panel.grid.major.x = element_line(color = "gray90", size = 0.5),
    panel.grid.minor.x = element_line(color = "gray90", size = 0.5),# Add gray major grid lines along x-axis
  )


anova_M_plot_21_40 = t_anova_M$plot_diff_abund(
  use_number = 21:40, 
  add_sig = TRUE, 
  coord_flip = TRUE,
  color_values = color_palette,
  barwidth = 0.8,
  text_y_size = 25,
  group_order = c("A_35_5", "A_35_7", "A_35_9",
                  "A_45_5", "A_45_7", "A_45_9",
                  "B_35_5", "B_35_7", "B_35_9",
                  "B_45_5", "B_45_7", "B_45_9")) +
  labs(
    x = "Genus",   # Customize your x-axis label
    y = "Relative Abundance", 
  ) +
  theme_pubr() +  # Apply the publication-ready theme
  theme(
    # title = element_text(size = 30),
    axis.text.x = element_text(size = 30),   # Increase x-axis text size
    axis.title.x = element_text(size = 30),  # Increase x-axis title size
    axis.title.y = element_text(size = 35),  # Increase x-axis title size
    axis.text.y = element_text(size = 30), 
    legend.text = element_text(size = 30),   # Increase legend text size
    legend.title = element_text(size = 35),
    legend.position = "right",# Increase legend title size
    panel.grid.major.x = element_line(color = "gray90", size = 0.5),
    panel.grid.minor.x = element_line(color = "gray90", size = 0.5),# Add gray major grid lines along x-axis
  )



png("Figures_Manuscript/FigS8_anova_F_1_20.png", width = 15, height = 40, units = "in", res = 300)
print(anova_F_plot_1_20)
dev.off()


png("Figures_Manuscript/FigS8_anova_F_21_40.png", width = 15, height = 40, units = "in", res = 300)
print(anova_F_plot_21_40)
dev.off()

png("Figures_Manuscript/FigS8_anova_M_1_20.png", width = 15, height = 40, units = "in", res = 300)
print(anova_M_plot_1_20)
dev.off()

png("Figures_Manuscript/FigS8_anova_M_21_40.png", width = 15, height = 40, units = "in", res = 300)
print(anova_M_plot_21_40)
dev.off()






################################################################################

################################################################################
#### FOOD WASTE NETWORK
################################################################################


## CREATE DATA SETS PER FEEDSTOCK

# View(dataset_rar_e_F$sample_table)

dataset_rar_e_F <- clone(dataset_rar_e)
dataset_rar_e_M <- clone(dataset_rar_e)
dataset_rar_e_F$sample_table <- subset(dataset_rar_e$sample_table,feedstock == "F");dataset_rar_e_F
dataset_rar_e_M$sample_table <- subset(dataset_rar_e$sample_table,feedstock == "M");dataset_rar_e_M

dataset_rar_e_F$tidy_dataset();dataset_rar_e_F
dataset_rar_e_M$tidy_dataset();dataset_rar_e_M


# if(!require("WGCNA")) install.packages("WGCNA", repos = BiocManager::repositories())
t_F <- trans_network$new(dataset = dataset_rar_e_F, cor_method = "spearman",taxa_level = "Genus", use_WGCNA_pearson_spearman = TRUE, filter_thres = 0.0001)

# network <- trans_network$new(dataset = dataset_rar_e_F, cal_cor = "base", taxa_level = "Genus", filter_thres = 0.0001, cor_method = "spearman")

# construct network; require igraph package
t_F$cal_network(COR_p_thres = 0.01, COR_optimization = TRUE, COR_cut = 0.7)



t_F$cal_module(method = "cluster_fast_greedy")

# get node properties
t_F$get_node_table(node_roles = TRUE)

t_F$cal_eigen()
# return t1$res_eigen

t_F$cal_network_attr()
View(t_F$res_network_attr)
t_F$get_edge_table()
t_F$plot_taxa_roles(use_type = 1)


View(t_F$res_edge_table)

t_F$get_node_table(node_roles = TRUE)

fw_node_table<-t_F$res_node_table

write.csv(fw_node_table,'fw_node_table_no55.csv') 



# plot node roles with phylum information
t_F$plot_taxa_roles(use_type = 2) 
# View(t1)

t_F$get_adjacency_matrix()
View(t_F$res_adjacency_matrix)

# View(t1$plot_taxa_roles)

t_F$cal_eigen()

t_F$res_eigen_expla


# create trans_env object
t_F_cor <- trans_env$new(dataset = dataset_rar_e_F, add_data = acids_env)
# calculate correlations
t_F_cor$cal_cor(add_abund_table = t1$res_eigen)
# plot the correlation heatmap
t_F_cor$plot_cor()

# t1$res_node_table


t_F$save_network(filepath = "network_fw_genus_spearman_no55.gexf")


# return t1$res_eigen



#########
### MANURE  NETWORK
#####


## CREATE DATA SETS PER FEEDSTOCK

dataset_rar_e_M <- clone(dataset_rar_e)
dataset_rar_e_M$sample_table <- subset(dataset_rar_e$sample_table,feedstock == "M");dataset_rar_e_M
dataset_rar_e_M$tidy_dataset();dataset_rar_e_M

# if(!require("WGCNA")) install.packages("WGCNA", repos = BiocManager::repositories())
t_M <- trans_network$new(dataset = dataset_rar_e_M, cor_method = "spearman",taxa_level = "Genus", use_WGCNA_pearson_spearman = TRUE, filter_thres = 0.0001)

# network <- trans_network$new(dataset = dataset_rar_e_F, cal_cor = "base", taxa_level = "Genus", filter_thres = 0.0001, cor_method = "spearman")

# construct network; require igraph package
t_M$cal_network(COR_p_thres = 0.01, COR_optimization = TRUE, COR_cut = 0.7)



t_M$cal_module(method = "cluster_fast_greedy")

# get node properties
t_M$get_node_table(node_roles = TRUE)

t_M$cal_eigen()
# return t1$res_eigen

t_M$cal_network_attr()
t_M$res_network_attr

t_M$get_node_table(node_roles = TRUE)

m_node_table<-t_M$res_node_table

write.csv(m_node_table,'m_node_table_no55.csv') 


# plot node roles with phylum information
t_M$plot_taxa_roles(use_type = 1)

t_M$get_edge_table()

t_M$cal_eigen()

t_M$res_eigen_expla

t_M$cal_network_attr()
View(t_M$res_network_attr)
# create trans_env object
t_M_cor <- trans_env$new(dataset = dataset_rar_e_M, add_data = acids_env)
# calculate correlations
t_M_cor$cal_cor(add_abund_table = t2$res_eigen)
# plot the correlation heatmap
t_M_cor$plot_cor()

t_M$res_node_table


t2$save_network(filepath = "network_M_genus_spearman_no55.gexf")

#### IDENTIFY IMPORTANT NODES AND COMPARE TAXA ROLES FOR EACH NETWORK

t_F$plot_taxa_roles(use_type = 1)
t_M$plot_taxa_roles(use_type = 1)

t_M$get_edge_table()
t_F$get_edge_table()

fw_edge_table=t_F$res_edge_table
m_edge_table=t_M$res_edge_table


# List of data frames to save into different sheets
edgenode_list <- list(
  m_node_table = m_node_table,
  fw_node_table = fw_node_table,
  fw_edge_table = fw_edge_table,
  m_edge_table = m_edge_table
)

library(writexl)

# Write to Excel with different sheets
write_xlsx(edgenode_list, "edge_node_table_fw_m_networks.xlsx")




#########
## CORRELATION MODULES PER FEEDSTOCK
#########

library(tidyr)

m_corplot = t_M$plot_cor()[["data"]]
f_corplot = t_F$plot_cor()[["data"]]

plot_cor

# Filter out unwanted columns and create the correlation matrix
m_corplot_wide_cor <- m_corplot %>%
  select(Taxa, Env, Correlation) %>%  # Keep only relevant columns
  pivot_wider(names_from = Env, values_from = Correlation)

m_corplot_wide_sig <- m_corplot %>%
  select(Taxa, Env, Significance) %>%  # Keep only relevant columns
  pivot_wider(names_from = Env, values_from = Significance)

m_corplot_wide_p <- m_corplot %>%
  select(Taxa, Env, AdjPvalue) %>%  # Keep only relevant columns
  pivot_wider(names_from = Env, values_from = AdjPvalue)

# Filter out unwanted columns and create the correlation matrix
f_corplot_wide_cor <- f_corplot %>%
  select(Taxa, Env, Correlation) %>%  # Keep only relevant columns
  pivot_wider(names_from = Env, values_from = Correlation)

f_corplot_wide_sig <- f_corplot %>%
  select(Taxa, Env, Significance) %>%  # Keep only relevant columns
  pivot_wider(names_from = Env, values_from = Significance)

f_corplot_wide_p <- f_corplot %>%
  select(Taxa, Env, AdjPvalue) %>%  # Keep only relevant columns
  pivot_wider(names_from = Env, values_from = AdjPvalue)




# List of data frames to save into different sheets
df_list <- list(
  m_corplot_wide_cor = m_corplot_wide_cor,
  m_corplot_wide_sig = m_corplot_wide_sig,
  m_corplot_wide_p = m_corplot_wide_p,
  f_corplot_wide_cor = f_corplot_wide_cor,
  f_corplot_wide_sig = f_corplot_wide_sig,
  f_corplot_wide_p = f_corplot_wide_p
)
library(writexl)

# Write to Excel with different sheets
write_xlsx(df_list, "correlation_data_no55.xlsx")






###########
#### ALPHA DIVERSITY
#########


tmp <- clone(dataset_rar_e)
tmp$cal_alphadiv()
# View(tmp$alpha_diversity)




####  APLHA DIVERSITYY FIXED EFFECTS

t1 <- trans_alpha$new(dataset = dataset_rar_e)
t1$cal_diff(method = "lme", formula = "inoculum+temperature+feedstock+pH + (1|reactor)", return_model = TRUE)
View(t1$res_diff)
t1$plot_alpha()


anova(t1$res_model$Shannon)


library(sjPlot)
library(sjlabelled)
library(sjmisc)
library(ggplot2)
theme_set(theme_sjplot())
library(ggplot2)  # Ensure ggplot2 is loaded

shannon_model=t1$res_model$Shannon
# Generate the plot and add a title
p <- plot_model(shannon_model, vline.color = "red") +
  ggtitle("Fixed effects on Shannon Diversity")  # Replace with your desired title

# Display the plot
print(p)



shannon_model

### FIXED EFFECTS PLOT FOR SHANNON DIVERSITY 


t1 <- trans_alpha$new(dataset = dataset_rar_e, group = "feedstock")
t2 <- trans_alpha$new(dataset = dataset_rar_e, group = "pH")
t3 <- trans_alpha$new(dataset = dataset_rar_e, group = "temperature")
t4 <- trans_alpha$new(dataset = dataset_rar_e, group = "inoculum")

t1$cal_diff(method = "wilcox")
t2$cal_diff(method = "wilcox")
t3$cal_diff(method = "wilcox")
t4$cal_diff(method = "wilcox")

library(gridExtra)

# Generate the alpha diversity plots
p1 <- t1$plot_alpha(measure = "Shannon")
p2 <- t2$plot_alpha(measure = "Shannon", color_values = c("red3","yellow3","blue3"))
p3 <- t3$plot_alpha(measure = "Shannon",color_values = c("black","purple"))
p4 <- t4$plot_alpha(measure = "Shannon",color_values = c("darkgreen","darkorange"))


# Arrange the plots in a single panel ### 
grid.arrange(p1, p2, p3, p4, ncol = 2)


#### ANOVA TABLE APLHA DIVERSITYY INTERACTION EFFECTS

t1 <- trans_alpha$new(dataset = dataset_rar_e, group = "inoculum_temperature_pH", by_group ="feedstock")
t1$cal_diff(method = "anova")
t1$plot_alpha(measure = "Shannon")


t1 <- trans_alpha$new(dataset = dataset_rar_e, group = "temperature", by_group ="pH")
t1$cal_diff(method = "anova")
t1$plot_alpha(measure = "Shannon")



################################################################################
#### NETWORK ANALYSIS - NETWORK COMPARISON
################################################################################


library(microeco)
library(meconetcomp)
# use pipe operator in magrittr package
library(magrittr)
library(igraph)
library(ggplot2)
theme_set(theme_bw())

all_networks <- list()


# Initialize an empty list to store the networks
all_networks <- list() #feedstock

# Get the unique values from the column 'inoculum_temperature_feedstock_pH'
# unique_values <- unique(dataset_rar_e$sample_table$inoculum_temperature_feedstock_pH)
unique_values <- unique(dataset_rar_e$sample_table$feedstock)



# all_networks_trt <- list()
# all_networks <- list()
# dataset_rar_e
# Iterate over each unique value
for (value in unique_values) {
  
  # Clone the dataset to create a deep copy
  tmp <- clone(dataset_rar_e)
  
  # Subset the sample table based on the current value
  # tmp$sample_table %<>% subset(inoculum_temperature_feedstock_pH == value)
  tmp$sample_table %<>% subset(feedstock == value)
  
  
  # Tidy the dataset to remove unnecessary features based on the subset
  tmp$tidy_dataset()
  
  # Create a new trans_network object with the filtered dataset
  tmp <- trans_network$new(dataset = tmp, cor_method = "spearman", filter_thres = 0.0001)
  
  # Calculate the network with correlation thresholds and optimization
  tmp$cal_network(COR_p_thres = 0.01, COR_optimization = TRUE, COR_cut = 0.7,
                  add_taxa_name = "Genus")
  
  tmp$cal_module(method = "cluster_fast_greedy")
  tmp$get_edge_table()
  tmp$cal_eigen()
  tmp$get_node_table()
  tmp$get_edge_table()
  
  
  
  
  # Create a filename based on the network name
  filename <- paste0("network_", value, "_sub_genus.gexf")
  
  # Save the network using the save_network method
  tmp$save_network(filepath = filename)
  
  # Store the network in the list with a name corresponding to the current value
  # all_networks_trt[[value]] <- tmp
  all_networks[[value]] <- tmp
  
}


View(all_networks_trt$B_35_F_5$res_node_table)


tmp_feedstock<- cal_network_attr(all_networks)
write.csv(tmp_feedstock,'network_properties_feedstock.csv')

tmp_all_trt<- cal_network_attr(all_networks_trt)
write.csv(tmp_feedstock,'network_properties_all_trt.csv')



# write.csv(tmp_property,'network_properties.csv') 



# HEATMAP FOR NETWORK PROPERTIES AND CORRELATION WITH FATTY ACIDS

# tmp <- subnet_property(all_networks)
tmp <- subnet_property(all_networks_trt)


rownames(tmp) <- tmp[, 2]
tmp <- tmp[, -c(1:2)]
tmp1 <- trans_env$new(dataset = dataset_rar_e, add_data = acids_env)
# View(tmp1)
tmp1$cal_cor(use_data = "other", by_group = "inoculum_temperature_feedstock_pH", add_abund_table = tmp, cor_method = "spearman")
g1 <- tmp1$plot_cor()
g1


                         
                         
## NETWORK CONNECTIVITY  FEEDSTOCK EFFECT

tmp <- robustness$new(all_networks, remove_strategy = c("edge_rand", "edge_strong", "node_rand", "node_degree_high"), 
                      remove_ratio = seq(0, 0.99, 0.1), measure = c("Eff", "Eigen", "Pcr"), run = 10)

tmp1 <- tmp$res_table %>% .[.$remove_strategy == "node_rand" & .$measure == "Eigen", ]
t1 <- trans_env$new(dataset = NULL, add_data = tmp1)
t1$dataset$sample_table <- t1$data_env
t1$plot_scatterfit(x = "remove_ratio", y = "value",line_se = FALSE, type = "cor", group = "Network") + 
  xlab("Ratio of randomly removed nodes") + ylab("Network connectivity") + theme(axis.title = element_text(size = 15))+
   scale_color_manual(values = c('#D95F02','#1B9E77'))



## NETWORK CONNECTIVITY ALL OPERATING CONDITIONS EFFECT

tmp <- robustness$new(all_networks_trt, remove_strategy = c("edge_rand", "edge_strong", "node_rand", "node_degree_high"), 
                      remove_ratio = seq(0, 0.99, 0.1), measure = c("Eff", "Eigen", "Pcr"), run = 10)

tmp1 <- tmp$res_table %>% .[.$remove_strategy == "node_rand" & .$measure == "Eigen", ]


t1 <- trans_env$new(dataset = NULL, add_data = tmp1)
t1$dataset$sample_table <- t1$data_env

# Use dplyr to create a new column that represents the interaction of feedstock and pH
t1$dataset$sample_table <- t1$dataset$sample_table %>%
  separate(Network, into = c("inoc", "T", "feedstock", "pH"), sep = "_", remove = FALSE) %>%
  mutate(inoc_T = paste(inoc, T, sep = "_"))

    

# Ensure ggplot2 is loaded for facet plotting
library(ggplot2)


library(ggplot2)
library(ggpubr)


## Figure for manuscript
ggplot(t1$dataset$sample_table, aes(x = remove_ratio, y = value, color = interaction(T, inoc))) + 
  geom_point() + 
  geom_smooth(method = "lm", se = TRUE, aes(fill = interaction(T, inoc))) +  # Add fill aesthetic for SE color
  xlab("Ratio of randomly removed nodes") + 
  ylab("Network connectivity") + 
  theme(axis.title = element_text(size = 15)) + 
  scale_color_manual(values = c( "green3","darkgreen",  "orange","orange4")) + 
  scale_fill_manual(values = c("green3","darkgreen",  "orange","orange4")) +  # Apply same colors for SE fill
  facet_wrap(~ feedstock + pH) +  # Faceting by 'feedstock' and 'pH'
  theme_pubr()  # Apply the 'theme_pubr' for clean publication-ready plots






############
#CIRCLE plots for correlations between taxa
#############

library("paletteer")

all_networks_trt$A_35_F_5$cal_sum_links(taxa_level = "Genus")  
all_networks_trt$B_35_F_5$cal_sum_links(taxa_level = "Genus")  
all_networks_trt$A_45_F_5$cal_sum_links(taxa_level = "Genus")  
all_networks_trt$B_45_F_5$cal_sum_links(taxa_level = "Genus")  


#print plots and save as PDF individual
all_networks_trt$A_35_F_5$plot_sum_links(method = "circlize", transparency = 0.5, plot_num = 20, annotationTrackHeight = circlize::mm_h(c(2, 2)),color_values = paletteer_d("ggthemes::Classic_20"))
all_networks_trt$B_35_F_5$plot_sum_links(method = "circlize", transparency = 0.5, plot_num = 20, annotationTrackHeight = circlize::mm_h(c(2, 2)),color_values = paletteer_d("ggthemes::Classic_20"))
all_networks_trt$A_45_F_5$plot_sum_links(method = "circlize", transparency = 0.5, plot_num = 20, annotationTrackHeight = circlize::mm_h(c(2, 2)),color_values = paletteer_d("ggthemes::Classic_20"))
all_networks_trt$B_45_F_5$plot_sum_links(method = "circlize", transparency = 0.5, plot_num = 20, annotationTrackHeight = circlize::mm_h(c(2, 2)),color_values = paletteer_d("ggthemes::Classic_20"))





# Open a PDF device for a single-page grid layout
pdf("FigureSX_F_35_5_Circle_plot_sum_links_grid_2x2.pdf", width = 20, height = 20) # Adjust width and height as needed

# Set up a 2x2 grid layout for the plots
par(mfrow = c(2, 2))

# Plot each network in the grid layout
all_networks_trt$A_35_F_5$plot_sum_links(
  method = "circlize",
  transparency = 0.5,
  plot_num = 20,
  annotationTrackHeight = circlize::mm_h(c(4, 4)),
  color_values = paletteer_d("ggthemes::Classic_20")
)

all_networks_trt$B_35_F_5$plot_sum_links(
  method = "circlize",
  transparency = 0.5,
  plot_num = 20,
  annotationTrackHeight = circlize::mm_h(c(4, 4)),
  color_values = paletteer_d("ggthemes::Classic_20")
)

all_networks_trt$A_45_F_5$plot_sum_links(
  method = "circlize",
  transparency = 0.5,
  plot_num = 20,
  annotationTrackHeight = circlize::mm_h(c(4, 4)),
  color_values = paletteer_d("ggthemes::Classic_20")
)

all_networks_trt$B_45_F_5$plot_sum_links(
  method = "circlize",
  transparency = 0.5,
  plot_num = 20,
  annotationTrackHeight = circlize::mm_h(c(4, 4)),
  color_values = paletteer_d("ggthemes::Classic_20")
)

# Close the PDF device
dev.off()



all_networks_trt$A_35_F_5$plot_sum_links(plot_pos = TRUE, plot_num = 20, color_values = RColorBrewer::brewer.pal(20, "Paired"))
all_networks_trt$B_35_F_5$plot_sum_links(plot_pos = TRUE, plot_num = 20, color_values =paletteer_d("ggthemes::Classic_20"))
all_networks_trt$A_45_F_5$plot_sum_links(plot_pos = TRUE, plot_num = 20, color_values = RColorBrewer::brewer.pal(20, "Paired"))

all_networks_trt$B_45_F_5$plot_sum_links(plot_pos = TRUE, plot_num = 20, color_values =  paletteer_d("ggthemes::Classic_20"))




install.packages("paletteer")
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

color_values = RColorBrewer::brewer.pal(8, "Dark2")
RColorBrewer::brewer.pal(12, "Dark2")

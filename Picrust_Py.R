#Per sequence contribution picrust


setwd("~/Projects/BETO/Microbial Analysis/2024_Interactions/R_Picrust")

source("functions.R")
source("libraries.R")


# Load the necessary libraries
#library(dplyr)
#library(tidyr)

#Read and format tables
load("MetaCyc_pathway_map.RData") #Metacyc Pathways ontology

metadata_0<-read.table("../16S/metadata_interactions.txt", sep='\t', header=T, comment="")

metadata<-metadata_0
# metadata<-read.table("DeLong_metadata_final.txt", sep='\t', header=T, comment="")
# metadata <- metadata_0[-1, ]

names(metadata)[names(metadata) == "sample.id"] <- "sample_id"
abundance_file_strat <- "../picrust2_updated/picrust2_out_pipeline_per_seq_con/pathways_out/path_abun_contrib.tsv" #path

path_abun_contrib <- read.delim(abundance_file_strat, header = TRUE)
# path_abun_contrib_psc <- read.delim("path_abun_contrib_psc.tsv", header = TRUE)


colnames(path_abun_contrib)[1] <- "sample_id"
# colnames(path_abun_contrib_psc)[1] <- "sample_id"
colnames(path_abun_contrib)[2] <- "path_id"
# colnames(path_abun_contrib_psc)[2] <- "path_id"

MetaCyc_pathway <- MetaCyc_pathway_map %>% tibble::rownames_to_column(var = "path_id")  
path_abun <- path_abun_contrib %>% left_join(MetaCyc_pathway, by = "path_id")  # Make sure 'function' matches in both data frames
# path_abun_psc <- path_abun_contrib_psc %>% left_join(MetaCyc_pathway, by = "path_id")  # Make sure 'function' matches in both data frames

##### Format master table taxonomy names
##### 
library(readr)

master_table_clean_rel <- read_csv("vault/master.table.clean_rel.csv",col_types = cols(...1 = col_integer())) 
names(master_table_clean_rel)[names(master_table_clean_rel) == "...1"] <- "number"
names(master_table_clean_rel)[names(master_table_clean_rel) == "FeatureID"] <- "taxon"
master_table_clean_rel <- as.data.frame(master_table_clean_rel)

master_table_clean_rel$int.slv.lws.lvl[master_table_clean_rel$int.slv.lws.lvl == "Species"] <- "(s)"
master_table_clean_rel$int.slv.lws.lvl[master_table_clean_rel$int.slv.lws.lvl == "Genus"] <- "(g)"
master_table_clean_rel$int.slv.lws.lvl[master_table_clean_rel$int.slv.lws.lvl == "Family"] <- "(f)"
master_table_clean_rel$int.slv.lws.lvl[master_table_clean_rel$int.slv.lws.lvl == "Order"] <- "(o)"
master_table_clean_rel$int.slv.lws.lvl[master_table_clean_rel$int.slv.lws.lvl == "Class"] <- "(c)"
master_table_clean_rel$int.slv.lws.lvl[master_table_clean_rel$int.slv.lws.lvl == "Phylum"] <- "(p)"
master_table_clean_rel$int.slv.lws.lvl[master_table_clean_rel$int.slv.lws.lvl == "Kingdom"] <- "(d)"

master_table_clean_rel$int.ggs.lws.lvl[master_table_clean_rel$int.ggs.lws.lvl == "Species"] <- "(s)"
master_table_clean_rel$int.ggs.lws.lvl[master_table_clean_rel$int.ggs.lws.lvl == "Genus"] <- "(g)"
master_table_clean_rel$int.ggs.lws.lvl[master_table_clean_rel$int.ggs.lws.lvl == "Family"] <- "(f)"
master_table_clean_rel$int.ggs.lws.lvl[master_table_clean_rel$int.ggs.lws.lvl == "Order"] <- "(o)"
master_table_clean_rel$int.ggs.lws.lvl[master_table_clean_rel$int.ggs.lws.lvl == "Class"] <- "(c)"
master_table_clean_rel$int.ggs.lws.lvl[master_table_clean_rel$int.ggs.lws.lvl == "Phylum"] <- "(p)"
master_table_clean_rel$int.ggs.lws.lvl[master_table_clean_rel$int.ggs.lws.lvl == "Kingdom"] <- "(d)"

master_table_clean_rel$ASV <- sapply(master_table_clean_rel$number, function(x) paste('ASV',x))
mastertable<-master_table_clean_rel

mastertable$gg.lws.taxon <- paste(mastertable$ASV,
                                  mastertable$int.ggs.lws.txn,
                                  mastertable$int.ggs.lws.lvl,
                                  sep = " ")

mastertable$slv.lws.taxon <- paste(mastertable$ASV,
                                  mastertable$int.slv.lws.txn,
                                  mastertable$int.slv.lws.lvl,
                                  sep = " ")

##### Input data to plot
######

metadata_0 <- metadata %>% filter(temperature != "55")

specific_functions <- c(
  "GLYCOCAT-PWY", 
  "PWY-6737", 
  "GALACTUROCAT-PWY", 
  "FUCCAT-PWY", 
  "PWY-6901", 
  "PWY-6507", 
  "PWY-5677",
  "PWY-5676", 
  "PWY-5100", 
  "ANAEROFRUCAT-PWY", 
  "P124-PWY", 
  "FASYN-ELONG-PWY", 
  "FASYN-INITIAL-PWY")

metadata_F <- metadata_0 %>% filter(feedstock == "F")
metadata_M <- metadata_0 %>% filter(feedstock == "M")

## Metadata reference detaframe
metadata_df <- as.data.frame(metadata_F) ## change here as reference to choose what data you want to plot



x=10 #change to desired top X for visualization
library(paletteer)

colors<-paletteer_d("ggsci::category20_d3")
path_abun_table=path_abun ## choose reference table here
# install.packages("paletteer")


# Initialize an empty list to store plots, need to change colors as well
plots_list <- list()
missing_path_ids <- c()

head(fun_tax_df)

for (specific_function in specific_functions) {
  
  if(!specific_function %in% path_abun_table$path_id) {
    # If the specific_function is not in path_abun, skip this iteration
    missing_path_ids <- c(missing_path_ids, specific_function)
    next
  }
  
  
  description <- MetaCyc_pathway %>%
    filter(path_id == specific_function) %>%
    pull(pathway) %>% unique()  # Assure you get a single description
  
  # If specific_description is NULL or length 0, provide a fallback
  if (is.null(description) || length(description) == 0) {
    description <- specific_function
  } else {
    description <- description[1]  # In case of multiple, take the first
  }
  
  # Filter for the specific function
  filtered_df <- path_abun_table %>% filter(path_id == specific_function)
  #filtered_df_psc <- path_abun_contrib_psc %>% filter(function. == specific_function)
  
  merged_df <- filtered_df %>% inner_join(metadata_df, by = "sample_id")
  
  fun_tax_df <- merged_df %>% left_join(mastertable %>% select(taxon, slv.lws.taxon,gg.lws.taxon), by = "taxon")
  
  top_x_taxa<- data_sum(fun_tax_df,varname="taxon_rel_function_abun",
                        groupnames=c("taxon","slv.lws.taxon","gg.lws.taxon")) %>% arrange(desc(taxon_rel_function_abun)) %>% head(n = x)
  
  top_taxa_gg <- top_x_taxa$gg.lws.taxon
  top_taxa_slv <- top_x_taxa$slv.lws.taxon
  
  # Mark taxa not in the top 10 as "Other"
  fun_tax_df$gg_taxon <- ifelse(fun_tax_df$gg.lws.taxon %in% top_taxa_gg, fun_tax_df$gg.lws.taxon, "Other")
  fun_tax_df$silva_taxon <- ifelse(fun_tax_df$slv.lws.taxon %in% top_taxa_slv, fun_tax_df$slv.lws.taxon, "Other")
  
  
  # Step 2: Sum taxon_rel_function_abun for each sample_id and taxon, including "Other"
  fun_tax_df_sum<- data_sum(fun_tax_df,
                            varname="taxon_rel_function_abun",
                            groupnames=c("sample_id",
                                         "day.2",
                                         "gg_taxon",
                                         "silva_taxon",
                                         "inoculum",
                                         "pH", 
                                         "day", 
                                         "temperature", 
                                         "feedstock", 
                                         "path_id",
                                         "pathway"))

  
  fun_tax_df_summary <- data_summary(fun_tax_df_sum, 
                                     varname="taxon_rel_function_abun",
                                     groupnames=c("gg_taxon",
                                                  "silva_taxon",
                                                  "inoculum",
                                                  "pH", 
                                                  "day",
                                                  "day.2",
                                                  "temperature", 
                                                  "feedstock",
                                                  "path_id",
                                                  "pathway"))
  
  
  fun_tax_df_summary$sample <- paste(fun_tax_df_summary$feedstock, 
                                     fun_tax_df_summary$temperature,
                                     fun_tax_df_summary$pH,
                                     fun_tax_df_summary$inoculum,
                                     fun_tax_df_summary$day,
                                     sep = "_")
  
  fun_tax_df_summary$pH_I_D <- paste(fun_tax_df_summary$pH,
                                     fun_tax_df_summary$inoculum,
                                     fun_tax_df_summary$day,
                                     sep = "_")
  
  fun_tax_df_summary$T_pH_I <- paste(fun_tax_df_summary$temperature,
                                     fun_tax_df_summary$pH,
                                     fun_tax_df_summary$inoculum,
                                     sep = "_")
  
  plot_df <- fun_tax_df_summary
  
  library(RColorBrewer)
  
  #colors
  #colors <-c( "aquamarine", "olivedrab", "firebrick", "gold", "purple", "steelblue2", "red", "blue4", "yellow3","coral","grey")
  #colors <- c("#386cb9","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33","blue4","grey")
  
  # Generate the plot
  
  plot <- ggplot(plot_df, aes(x = day, y = taxon_rel_function_abun_mean, fill = silva_taxon)) + 
    geom_bar(stat = "identity", position = "stack") +
    theme_Publication() +  
    labs(y = "tax_rel_path_ab", title = description) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          legend.position = "right")+
    scale_fill_manual(values = colors, name = " ") +
    guides(fill = guide_legend(ncol = 1)) +
    facet_grid(~ T_pH_I)
  
  # Store the plot in the list with the function name as the key
  plots_list[[specific_function]] <- plot

  
}

##### display plots
######

# Print missing pathway IDs note
if (length(missing_path_ids) > 0) {
  cat("Note: The following pathway IDs were not found:", paste(missing_path_ids, collapse = ", "), "\n")
}






  F_plot_list<-plots_list
 # M_plot_list<-plots_list

F_total=wrap_plots(F_plot_list[1:3], ncol = 1)
F_C4=wrap_plots(F_plot_list[7:8], ncol = 1)
F_MCFA_1=wrap_plots(F_plot_list[9:11], ncol = 1)
F_MCFA_2=wrap_plots(F_plot_list[12:13], ncol = 1)
M_total=wrap_plots(M_plot_list[4:6], ncol = 1)



 png("FigS_F_total_taxa_pwy.png", width = 10.5, height = 7, units = "in", res = 300)
 F_total
 dev.off()

 png("FigS_F_C4_taxa_pwy.png", width = 10.5, height = 7, units = "in", res = 300)
 F_C4
 dev.off()
 
 png("FigS_F_MCFA1_taxa_pwy.png", width = 10.5, height = 7, units = "in", res = 300)
 F_MCFA_1
 dev.off()
 
 png("FigS_F_MCFA2_taxa_pwy.png", width = 10.5, height = 7, units = "in", res = 300)
 F_MCFA_2
 dev.off()
 
 png("FigS_M_total_taxa_pwy.png", width = 10.5, height = 7, units = "in", res = 300)
 M_total
 dev.off()
 
 
 














# 
# # Vector of specific function names
# #specific_functions <- c("CENTFERM-PWY" , "PWY-8190" , "PWY-5022" , "P163-PWY" , "PWY-5677" , "PWY-5676" , "GLUDEG-II-PWY" , "P162-PWY")  # Update with your function names
# specific_functions <- c("PWY-5494" , "PROPFERM-PWY" , "PWY-7013" , "P108-PWY" , "PWY-5088" , "PWY-8188" , "PWY-8086")  # Update with your function names
# #specific_functions <- c("P124-PWY" , "P122-PWY" , "PWY-5481" , "P461-PWY" , "PWY-8274" , "ANAEROFRUCAT-PWY" , "PWY-5100" , "P41-PWY")  # Update with your function names
# #specific_functions <- c("PWY66-398","TCA","PWY-6969","REDCITCYC","PWY-7254","P105-PWY","PWY-5913","PWY-5690","TCA-1")
# specific_functions <- c("METH-ACETATE-PWY","CO2FORM-PWY","METHANOGENESIS-PWY","PWY-6830","METHFORM-PWY")
# 
# specific_functions <- c("PWY-5100", #acetate
#                         "P108-PWY", #Prop
#                         "PWY-5676", ## Butyrate
#                         "P122-PWY", ## Lactate 
#                         "FASYN-ELONG-PWY",
#                         "METHANOGENESIS-PWY")
# 
# 
# 
# 
# specific_functions <- c("PWY-6901")
# specific_functions <- c("PWY-5431")
# 
# specific_functions <- c("PWY-6901",
#                         "PWY-5431",
#                         "FUCCAT-PWY",
#                         "LEU-DEG2-PWY",
#                         "PWY-6338",
#                         "PWY-6590",
#                         "FERMENTATION-PWY") ### Degradation



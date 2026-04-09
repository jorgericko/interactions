

#SET WORKING DIRECTORY
setwd("~/Projects/BETO/Experiments_Inoculum/A_B_F_M")
####### If you are running the Data processing step from here (source(file="Data_Processing.R) below)
#### then running the libraries.R (source(file="libraries.R")) is redundant but doesn't hurt
source(file="libraries.R")
#### The following loads some functions we use to plot and look at our data (description is here and in the "functions.R" file)
source(file="functions.R")
library(dplyr)

library(ggpubr)


####################################
###  LINE PLOTS 
#################################






#load combined data
data_acids_0<-read.csv("data_acids_A_B.csv", comment="") #Combined data


# Define the conversion factors (gCOD per g of acid)
# Ensure the names match the labels in your 'acid' column (e.g., C2, C3, etc.)
cod_factors <- c(
  "C2" = 1.07,
  "C3" = 1.51,
  "ISOC4" = 1.82, # Isobutyric
  "C4" = 1.82,  # n-Butyric
  "ISOC5" = 2.04, # Isovaleric
  "C5" = 2.04,  # n-Valeric
  "ISOC6" = 2.21, # Isocaproic
  "C6" = 2.21,  # n-Caproic
  "ISOC7" = 2.34, # Isoheptanoic
  "C7" = 2.34   # Heptanoic 
) 

data_acids_00 <- data_acids_0 %>%
  mutate(conc_gCOD_L = conc_gL * cod_factors[acid])


data_acids_f1=filter(data_acids_00,rep !="D")
# data_acids_f1=data_acids_0
data_acids_f2=filter(data_acids_f1,rep !="INOC")
data_acids_f3=filter(data_acids_f2,T!=55)
data_acids_f4=filter(data_acids_f3,day != 9.7)%>%
  filter(between(day, 0, 13))  %>% filter(T!=55)



#############################
#### Data in gCOD/L ########
################################
data_acids_f5 <- data_acids_f4 %>%
  select(-conc_gL)

data_acids_COD <- pivot_wider(data_acids_f5, names_from = acid, values_from = conc_gCOD_L)

data_acids_COD[is.na(data_acids_COD)] = 0
data_acids_COD$total <- rowSums(data_acids_COD[,15:23])
data_acids_COD$day <- as.numeric(data_acids_COD$day)
data_acids_COD=as.data.frame(data_acids_COD)
data_acids_COD$pH <- as.factor(data_acids_COD$pH)


# Change here if you want to plot concentrations
all_data_acids_stacked_COD<-(data_acids_COD %>% pivot_longer(
  cols = C2:total,
  names_to = c("acid"),
  values_to = "conc_gCOD_L"))


# all_data_acids_stacked<-(data_acids %>% pivot_longer(
#   cols = C2:total,
#   names_to = c("acid"),
#   values_to = "yield"))
#  
# filtered_sd <- all_data_acids_stacked[!all_data_acids_stacked$acid %in% c("ISOC4", "ISOC5", "ISOC6", "C7"), ]
filtered_sd_COD <- all_data_acids_stacked_COD[!all_data_acids_stacked_COD$acid %in% c("ISOC4", "ISOC5", "ISOC6"), ]
filtered_sd_COD$T_with_degree <- paste0(filtered_sd_COD$T, "°C")
filtered_sd_COD$F_T <- interaction(filtered_sd_COD$T_with_degree, filtered_sd_COD$feedstock,sep = "_")
filtered_sd_COD$pH_plot <- paste0("pH ",filtered_sd_COD$pH)
filtered_sd_COD$pH_Inoculum <- interaction(filtered_sd_COD$pH_plot, filtered_sd_COD$inoc,sep = "_")
str(filtered_sd_COD)

head(filtered_sd_COD)

# 
# feedstock_load <- 15 #gCOD/L of feed
# 
# filtered_sd <- filtered_sd %>%
#   mutate(yield = conc_gL / feedstock_load)
# 


library(dplyr)

filtered_sd_COD <- filtered_sd_COD %>%
  mutate(yield = case_when(
    feedstock == "F" & inoc == "A" ~ conc_gCOD_L / 15.8,  # Specific case: F and A
    TRUE                           ~ conc_gCOD_L / 15   # Everything else
  ))



# Create the ggline plot
FA_plots_COD = ggline(
  filtered_sd_COD , 
  x = "day", 
  y = "conc_gCOD_L", 
  numeric.x.axis = TRUE, 
  xlab = "day", 
  ylab = "Fatty Acid Concentration  [gCOD/L]", 
  facet.by = c("acid","F_T"),
  scales = "free_y",
  error.plot = "errorbar", 
  add = c("mean_sd", "jitter"), 
  add.params = list(size = 0.3), 
  color = "pH_inoc", 
  shape = "pH_inoc", 
  point.size = 1, 
  palette = c("#A73030FF","#CD534CFF","#8F7700FF", "#EFC000FF","#003C67FF","#0073C2FF")
) + 
  theme_pubr() +
  labs_pubr()+  
  labs(color = "pH_Inoculum", shape = "pH_Inoculum") + labs_pubr(base_size = 12)+  
  theme(
    legend.position = "right",
    legend.key = element_rect(fill = "white", color = "black"),  # Rectangle legend key
    legend.key.size = unit(1.5, "lines"),  # Adjust size of legend keys
    axis.text.x = element_text(size = rel(0.8)),  # X-axis text size (default size * 0.8)
    axis.text.y = element_text(size = rel(0.8))   # Y-axis text size (default size * 0.8)
  )+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 4), labels = scales::number_format(accuracy = 1))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 3), labels = scales::number_format(accuracy = 0.2))



# Create the ggline plot
FA_plots_yield_COD = ggline(
  filtered_sd_COD , 
  x = "day", 
  y = "yield", 
  numeric.x.axis = TRUE, 
  xlab = "day", 
  ylab = "Fatty Acid Yield  [gCOD/gCODfed]", 
  facet.by = c("acid","F_T"),
  scales = "free_y",
  error.plot = "errorbar", 
  add = c("mean_sd", "jitter"), 
  add.params = list(size = 0.3), 
  color = "pH_inoc", 
  shape = "pH_inoc", 
  point.size = 1, 
  palette = c("#A73030FF","#CD534CFF","#8F7700FF", "#EFC000FF","#003C67FF","#0073C2FF")
) + 
  theme_pubr() +
  labs_pubr()+  
  labs(color = "pH_Inoculum", shape = "pH_Inoculum") + labs_pubr(base_size = 12)+  
  theme(
    legend.position = "right",
    legend.key = element_rect(fill = "white", color = "black"),  # Rectangle legend key
    legend.key.size = unit(1.5, "lines"),  # Adjust size of legend keys
    axis.text.x = element_text(size = rel(0.8)),  # X-axis text size (default size * 0.8)
    axis.text.y = element_text(size = rel(0.8))   # Y-axis text size (default size * 0.8)
  )+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 4), labels = scales::number_format(accuracy = 1))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4), labels = scales::number_format(accuracy = 0.01))



filtered_sd_total_COD <- filtered_sd_COD[filtered_sd_COD$acid == "total", ]

# Create the ggline plot for total
FA_plots_yield_total_COD = ggline(
  filtered_sd_total_COD , 
  x = "day", 
  y = "yield", 
  numeric.x.axis = TRUE, 
  xlab = "day", 
  ylab = "Fatty Acid Yield  [gCOD/gCODfed]", 
  facet.by = c("acid","F_T"),
  scales = "free_y",
  error.plot = "errorbar", 
  add = c("mean_sd", "jitter"), 
  add.params = list(size = 0.3), 
  color = "pH_inoc", 
  shape = "pH_inoc", 
  point.size = 1, 
  palette = c("#A73030FF","#CD534CFF","#8F7700FF", "#EFC000FF","#003C67FF","#0073C2FF")
) + 
  theme_pubr() +
  labs_pubr()+  
  labs(color = "pH_Inoculum", shape = "pH_Inoculum") + labs_pubr(base_size = 12)+  
  theme(
    legend.position = "right",
    legend.key = element_rect(fill = "white", color = "black"),  # Rectangle legend key
    legend.key.size = unit(1.5, "lines"),  # Adjust size of legend keys
    axis.text.x = element_text(size = rel(0.8)),  # X-axis text size (default size * 0.8)
    axis.text.y = element_text(size = rel(0.8))   # Y-axis text size (default size * 0.8)
  )+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 4), labels = scales::number_format(accuracy = 1))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4), labels = scales::number_format(accuracy = 0.01))



print(FA_plots_yield_total_COD)


print(FA_plots_yield_COD)


print(FA_plots_COD)


png("Final_Figures/fig_3_FA_YIELD_corrected_total_2_COD.png", width = 7, height = 3, units = "in", res = 300)
print(FA_plots_yield_total_COD)
dev.off()


png("Final_Figures/fig_3_FA_YIELD_corrected_COD.png", width = 7, height = 8, units = "in", res = 300)
print(FA_plots_yield_COD)
dev.off()



png("Final_Figures/fig_SX_FA_all_COD.png", width = 14, height = 7, units = "in", res = 300)
print(FA_plots_COD)
dev.off()










#############################
#### Data in g/L ########
################################

data_acids_f6 <- data_acids_f4 %>%
  select(-conc_gCOD_L)

data_acids <- pivot_wider(data_acids_f6, names_from = acid, values_from = conc_gL)
data_acids[is.na(data_acids)] = 0
data_acids$total <- rowSums(data_acids[,15:23])
data_acids$day <- as.numeric(data_acids$day)
data_acids=as.data.frame(data_acids)
data_acids$pH <- as.factor(data_acids$pH)


# Change here if you want to plot concentrations
all_data_acids_stacked<-(data_acids %>% pivot_longer(
  cols = C2:total,
  names_to = c("acid"),
  values_to = "conc_gL"))


# all_data_acids_stacked<-(data_acids %>% pivot_longer(
#   cols = C2:total,
#   names_to = c("acid"),
#   values_to = "yield"))
#  
# filtered_sd <- all_data_acids_stacked[!all_data_acids_stacked$acid %in% c("ISOC4", "ISOC5", "ISOC6", "C7"), ]
filtered_sd <- all_data_acids_stacked[!all_data_acids_stacked$acid %in% c("ISOC4", "ISOC5", "ISOC6"), ]
filtered_sd$T_with_degree <- paste0(filtered_sd$T, "°C")
filtered_sd$F_T <- interaction(filtered_sd$T_with_degree, filtered_sd$feedstock,sep = "_")
filtered_sd$pH_plot <- paste0("pH ",filtered_sd$pH)
filtered_sd$pH_Inoculum <- interaction(filtered_sd$pH_plot, filtered_sd$inoc,sep = "_")
str(filtered_sd)

head(filtered_sd)

# 
# feedstock_load <- 15 #gCOD/L of feed
# 
# filtered_sd <- filtered_sd %>%
#   mutate(yield = conc_gL / feedstock_load)
# 


library(dplyr)

filtered_sd <- filtered_sd %>%
  mutate(yield = case_when(
    feedstock == "F" & inoc == "A" ~ conc_gL / 15.8,  # Specific case: F and A
    TRUE                           ~ conc_gL / 15   # Everything else
  ))



# Create the ggline plot
FA_plots = ggline(
  filtered_sd , 
  x = "day", 
  y = "conc_gL", 
  numeric.x.axis = TRUE, 
  xlab = "day", 
  ylab = "Fatty Acid Concentration  [g/L]", 
  facet.by = c("acid","F_T"),
  scales = "free_y",
  error.plot = "errorbar", 
  add = c("mean_sd", "jitter"), 
  add.params = list(size = 0.3), 
  color = "pH_inoc", 
  shape = "pH_inoc", 
  point.size = 1, 
  palette = c("#A73030FF","#CD534CFF","#8F7700FF", "#EFC000FF","#003C67FF","#0073C2FF")
) + 
  theme_pubr() +
  labs_pubr()+  
  labs(color = "pH_Inoculum", shape = "pH_Inoculum") + labs_pubr(base_size = 12)+  
  theme(
    legend.position = "right",
    legend.key = element_rect(fill = "white", color = "black"),  # Rectangle legend key
    legend.key.size = unit(1.5, "lines"),  # Adjust size of legend keys
    axis.text.x = element_text(size = rel(0.8)),  # X-axis text size (default size * 0.8)
    axis.text.y = element_text(size = rel(0.8))   # Y-axis text size (default size * 0.8)
  )+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 4), labels = scales::number_format(accuracy = 1))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 3), labels = scales::number_format(accuracy = 0.2))



# Create the ggline plot
FA_plots_yield = ggline(
  filtered_sd , 
  x = "day", 
  y = "yield", 
  numeric.x.axis = TRUE, 
  xlab = "day", 
  ylab = "Fatty Acid Yield  [g/gCOD]", 
  facet.by = c("acid","F_T"),
  scales = "free_y",
  error.plot = "errorbar", 
  add = c("mean_sd", "jitter"), 
  add.params = list(size = 0.3), 
  color = "pH_inoc", 
  shape = "pH_inoc", 
  point.size = 1, 
  palette = c("#A73030FF","#CD534CFF","#8F7700FF", "#EFC000FF","#003C67FF","#0073C2FF")
) + 
  theme_pubr() +
  labs_pubr()+  
  labs(color = "pH_Inoculum", shape = "pH_Inoculum") + labs_pubr(base_size = 12)+  
  theme(
    legend.position = "right",
    legend.key = element_rect(fill = "white", color = "black"),  # Rectangle legend key
    legend.key.size = unit(1.5, "lines"),  # Adjust size of legend keys
    axis.text.x = element_text(size = rel(0.8)),  # X-axis text size (default size * 0.8)
    axis.text.y = element_text(size = rel(0.8))   # Y-axis text size (default size * 0.8)
  )+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 4), labels = scales::number_format(accuracy = 1))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4), labels = scales::number_format(accuracy = 0.01))



filtered_sd_total <- filtered_sd[filtered_sd$acid == "total", ]

# Create the ggline plot for total
FA_plots_yield_total = ggline(
  filtered_sd_total , 
  x = "day", 
  y = "yield", 
  numeric.x.axis = TRUE, 
  xlab = "day", 
  ylab = "Fatty Acid Yield  [g/gCOD]", 
  facet.by = c("acid","F_T"),
  scales = "free_y",
  error.plot = "errorbar", 
  add = c("mean_sd", "jitter"), 
  add.params = list(size = 0.3), 
  color = "pH_inoc", 
  shape = "pH_inoc", 
  point.size = 1, 
  palette = c("#A73030FF","#CD534CFF","#8F7700FF", "#EFC000FF","#003C67FF","#0073C2FF")
) + 
  theme_pubr() +
  labs_pubr()+  
  labs(color = "pH_Inoculum", shape = "pH_Inoculum") + labs_pubr(base_size = 12)+  
  theme(
    legend.position = "right",
    legend.key = element_rect(fill = "white", color = "black"),  # Rectangle legend key
    legend.key.size = unit(1.5, "lines"),  # Adjust size of legend keys
    axis.text.x = element_text(size = rel(0.8)),  # X-axis text size (default size * 0.8)
    axis.text.y = element_text(size = rel(0.8))   # Y-axis text size (default size * 0.8)
  )+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 4), labels = scales::number_format(accuracy = 1))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4), labels = scales::number_format(accuracy = 0.01))



print(FA_plots_yield_total)


print(FA_plots_yield)


print(FA_plots)


png("Final_Figures/fig_3_FA_YIELD_corrected_total_2.png", width = 7, height = 3, units = "in", res = 300)
print(FA_plots_yield_total)
dev.off()


png("Final_Figures/fig_3_FA_YIELD_corrected.png", width = 7, height = 8, units = "in", res = 300)
print(FA_plots_yield)
dev.off()


### nature width 7 in
png("Final_Figures/fig_3_FA.png", width = 7, height = 8, units = "in", res = 300)
print(FA_plots)
dev.off()

png("Final_Figures/fig_SX_FA_all.png", width = 14, height = 7, units = "in", res = 300)
print(FA_plots)
dev.off()

#################################

########################################
### GAS PRODUCTION
########################################
gas_data<-read.csv("gas_data_all.csv", comment="") #Combined data
gas_data$pH <- as.factor(gas_data$pH)
gas_data$I <- as.factor(gas_data$I)
gas_data$T <- as.factor(gas_data$T)
gas_data$F <- as.factor(gas_data$F)

gas_data$reactor <- paste(gas_data$I, gas_data$T,gas_data$F, gas_data$pH,gas_data$rep,sep = "_")
gas_data$reactor <- as.factor(gas_data$reactor)


gas_stacked<-(gas_data %>% pivot_longer(
  cols = CH4:total,
  names_to = c("product"),
  values_to = "mL"))

detach("package:dplyr", unload=TRUE) ### always deatach for data summary
gas_data_summary<- data_summary(gas_stacked, varname="mL",groupnames=c("I", "T", "F","pH","product"))
gas_data_summary$T_with_degree <- paste0(gas_data_summary$T, "°C")

gas_data_summary$F_T <- interaction(gas_data_summary$T_with_degree, gas_data_summary$F,sep = "_")
gas_data_summary_filter <- gas_data_summary[!gas_data_summary$product %in% c("H2", "CO2"),]
gas_data_summary_filter <- gas_data_summary_filter[!gas_data_summary_filter$T %in% c("55"),]




gas_p<- ggplot(gas_data_summary_filter, aes(x=I, y=mL, fill=pH)) + 
  geom_bar(stat="identity", color="black", position="dodge") + 
  geom_errorbar(aes(ymin=mL-sd, ymax=mL+sd), width=0.3,  position=position_dodge(0.9)) + 
  facet_grid(product ~ F_T)+
  xlab("Inoculum") + ylab("Cumulative Gas [mL]") +
  scale_fill_manual(breaks = c("5","7","9"), values=c("#CD534CFF", "#EFC000FF","#0073C2FF"))+ 
  theme_pubr(legend = c("top"))+ labs_pubr(base_size = 11)

gas_p

png("Final_Figures/fig_2_gas.png", width = 3.3, height = 4.5, units = "in", res = 300)
print(gas_p)
dev.off()

library(broom)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(lme4)
library(emmeans)
library(lmerTest)
library(pbkrtest)
library(ggpubr)

library(sjPlot)
library(sjlabelled)
library(sjmisc)
library(ggplot2)


gas_data=filter(gas_data,T!=55)

# Fit the model
gas_data$pH <- relevel(gas_data$pH, ref = "7")
model_CH4 <- lm(CH4~ pH+I+T+F, data = gas_data)
model_total <- lm(total~ pH+I+T+F, data = gas_data)

summary(model_CH4)
summary(model_total)

anova(model_CH4)
anova(model_total)

small_text_theme <- theme(
  plot.title = element_text(size = 10), # Adjust title size here
  axis.title = element_text(size = 10),
  axis.text = element_text(size = 1),
  legend.title = element_text(size = 8),
  legend.text = element_text(size = 7)
)


plot_gas_models<-plot_models(model_CH4, model_total, #grid.breaks = TRUE,
                         axis.labels = c("Feedstock = M","T = 45°C","Inoculum = B", "pH = 9","pH = 5"), 
                         p.shape = TRUE,
                         ci.lvl = 0.95,
                         p.threshold = c(0.05, 0.01, 0.001),
                         spacing = 0.7,
                         dot.size = 2,legend.title = "Dependent Variables",
                         colors = c("blue", "red"),
                         vline.color = "#061423")+ theme_pubr()+ labs_pubr()+
  theme(legend.position = "right",panel.grid.minor = element_line(color = "grey90"))+ 
  labs(y = "Coefficient Effect Size") 

plot_gas_models 

png("Final_Figures/fig_SX_fe_coef_gas.png", width = 6, height = 4, units = "in", res = 300)
print(plot_gas_models)
dev.off()




#################################

########################################
### DATA SUMMARIES
########################################
###################




####################
### DATA SUMMARIES YIELD ###
######################


head(filtered_sd_COD)


detach("package:dplyr", unload=TRUE) ### always deatach for data summary

all_data_acids_sum_yield_COD <- data_summary(filtered_sd_COD, varname="yield",groupnames=c("inoc","pH", "day","acid", "T", "feedstock"))

library(writexl)
write_xlsx(all_data_acids_sum_yield_COD, path = "~/Projects/BETO/Experiments_Inoculum/A_B_F_M/Final_Tables/table_SX_fa_summary_yield_COD.xlsx")







####################
### DATA SUMMARIES CONCENTRATIONS [G/L]  ###
######################



detach("package:dplyr", unload=TRUE) ### always deatach for data summary

all_data_acids_sum <- data_summary(filtered_sd, varname="conc_gL",groupnames=c("inoc","pH", "day","acid", "T", "feedstock"))

write.csv(all_data_acids_sum, "~/Projects/BETO/Experiments_Inoculum/A_B_F_M/Final_Tables/all_data_acids_summary.csv", row.names=FALSE)
library(writexl)

write_xlsx(all_data_acids_sum, path = "~/Projects/BETO/Experiments_Inoculum/A_B_F_M/Final_Tables/table_SX_fa_summary.xlsx")






#################################



########################################
### STATISTICAL ANALYSIS
########################################



# # Install jtools if not already installed
# if (!requireNamespace("jtools", quietly = TRUE)) {
#   install.packages("jtools")
# }
# # Install broom.mixed if not already installed
# if (!requireNamespace("broom.mixed", quietly = TRUE)) {
#   install.packages("broom.mixed")
# }
# 
# # Load jtools package
# library(jtools)
# Load required libraries

library(broom)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(lme4)
library(emmeans)
library(lmerTest)
library(pbkrtest)
library(ggpubr)




#################################

############################################################
###  FIXED EFFECTS MODEL CONSTRUCTION WITH RANDOM EFFECTS AND PLOT QQ - RESIDUALS FIG S2f
############################################################
library(broom)
library(dplyr)
library(ggplot2)
library(gridExtra)

library(lme4)
library(emmeans)
library(lmerTest)
library(pbkrtest)
library(ggpubr)

library(performance)
library(gridExtra)
library(grid)
library(ggplot2)


# Define the response variables

data_acids_stats =data_acids
data_acids_stats$reactor <- paste(data_acids_stats$inoc, data_acids_stats$T,data_acids_stats$feedstock, data_acids_stats$pH,data_acids_stats$rep,sep = "_")
data_acids_stats=as.data.frame(data_acids_stats)
data_acids_stats$reactor <- paste(data_acids_stats$inoc, data_acids_stats$T,data_acids_stats$feedstock, data_acids_stats$pH,data_acids_stats$rep,sep = "_")


data_acids_stats_end <- data_acids_stats %>% filter(between(day, 0.1, 12))
data_acids_stats_end$pH <- as.factor(data_acids_stats_end$pH)
data_acids_stats_end$T <- as.factor(data_acids_stats_end$T)
data_acids_stats_end$inoc <- as.factor(data_acids_stats_end$inoc)
data_acids_stats_end$feedstock <- as.factor(data_acids_stats_end$feedstock)
data_acids_stats_end$reactor <- factor(data_acids_stats_end$reactor)

data_acids_stats_end$pH <- relevel(data_acids_stats_end$pH, ref = "7")
data_acids_stats_end$inoc <- relevel(data_acids_stats_end$inoc, ref = "B")
data_acids_stats_end$T <- relevel(data_acids_stats_end$T, ref = "45")
data_acids_stats_end$feedstock <- relevel(data_acids_stats_end$feedstock, ref = "M")
#day is a continous variable
response_vars <- c("C2", "C3", "C4", "C5", "C6","C7", "total")



str(data_acids_stats_end)


# Loop through each response variable
fe_models <- list()
fe_anovas <- list()
fe_models_summary<- list()


# Apply effect coding
# contrasts(data_acids_stats_end$pH) <- contr.sum
# contrasts(data_acids_stats_end$inoc) <- contr.sum
# contrasts(data_acids_stats_end$T) <- contr.sum
# contrasts(data_acids_stats_end$feedstock) <- contr.sum


for (response in response_vars) {
  
  
  formula <- as.formula(paste(response, "~ pH + inoc + T + feedstock+ (1 | reactor)"))
  
  # Fit the model
  model <- lmer(formula, data = data_acids_stats_end)
  
  # Summarize the model
  summary_model <- summary(model)
  anova_model <- anova(model)
  # Store the model and summary
  fe_models[[response]] <- model
  fe_models_summary[[response]] <- summary_model
  fe_anovas[[response]] <- anova_model
}


# model <- nlmer(C2 ~ pH * inoc * T * feedstock + (1 | reactor), data = data_acids_stats_end,)
library(ggResidpanel)

# Function to create a residual plot with title
create_residual_plot_with_title <- function(model, title) {
  residuals_plot <- resid_panel(model)
  grid.arrange(residuals_plot, top = textGrob(title, gp = gpar(fontsize = 14, fontface = "bold")))
}

# Create and capture residual plots for each model
feplot_C2 <- create_residual_plot_with_title(fe_models$C2, "C2")
feplot_C3 <- create_residual_plot_with_title(fe_models$C3, "C3")
feplot_C4 <- create_residual_plot_with_title(fe_models$C4, "C4")
feplot_C5 <- create_residual_plot_with_title(fe_models$C5, "C5")
feplot_C6 <- create_residual_plot_with_title(fe_models$C6, "C6")
feplot_C7 <- create_residual_plot_with_title(fe_models$C7, "C7")
feplot_total <- create_residual_plot_with_title(fe_models$total, "Total")


png("Final_Figures/fig_SX_resid_fe_with_randomeffects.png", width = 9, height = 9, units = "in", res = 300)
resid_fe_combined_plot <- grid.arrange(feplot_C2, feplot_C3, feplot_C4,feplot_C5, feplot_C6,feplot_C7, feplot_total,  ncol = 3)
dev.off()





############################################################
###   MODEL SUMMARY TABLES FOR FIXED EFFECTS WITH RANDOM EFFECTS ANOVA TABLE S1 and SUMMARY MODEL TABLES2
############################################################
library(writexl)
library(dplyr)
library(writexl)
library(dplyr)
library(broom.mixed)
library(sjPlot)


# Function to add significance levels
add_significance <- function(df) {
  df %>%
    mutate(significance = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01 ~ "**",
      p.value < 0.05 ~ "*",
      p.value < 0.1 ~ ".",
      TRUE ~ ""
    ))
}

# Function to process a list of models or ANOVAs
process_list <- function(model_list) {
  tidy_list <- lapply(model_list, function(model) {
    tidy_result <- tidy(model)
    tidy_result <- add_significance(tidy_result)
    return(tidy_result)
  })

  
  combined_df <- do.call(rbind, lapply(seq_along(tidy_list), function(i) {
    data <- tidy_list[[i]]
    data$model <- names(model_list)[i]
    return(data)
  }))
  
  combined_df <- combined_df %>%
    mutate(across(where(is.numeric), ~ signif(., 3)))
  
  return(combined_df)
}

# Remove NULL elements from the fe_anovas list
fe_anovas <- fe_anovas[!sapply(fe_anovas, is.null)]

# Process the cleaned fe_anovas list
anovas_processed <- process_list(fe_anovas)
#fe_models_summary <- process_list(fe_models_summary)


# Create a list of data frames to export to different sheets
data_to_export <- list("table_SX_fe_ANOVAs" = anovas_processed)


# table_SX_fe_ANOVAs
write_xlsx(data_to_export, path = "Final_Tables/table_SX_fe_ANOVAs.xlsx")
# load required package

#print and paste in excel_see word file created manually: table_SX_model_summaries_fe_and_int.docx
tab_model(fe_models$C2,fe_models$C3,fe_models$C4) ## Figure will be from this table table_SX_fe_model_summary
tab_model(fe_models$C5,fe_models$C6,fe_models$C7,fe_models$total)

#################################


#################################

############################################################
###  FIG 4 FIXED EFFECT MODEL PLOTS WITH RANDOM EFFECTS (NO RANDOM EFFECTS COMMENTED)
############################################################
library(sjPlot)
library(sjlabelled)
library(sjmisc)
library(ggplot2)

small_text_theme <- theme(
  plot.title = element_text(size = 10), # Adjust title size here
  axis.title = element_text(size = 10),
  axis.text = element_text(size = 10),
  legend.title = element_text(size = 8),
  legend.text = element_text(size = 7)
)


plot_all_fe<-plot_models(fe_models$C2, fe_models$C3, fe_models$C4,
                         fe_models$C5, fe_models$C6,fe_models$C7, fe_models$total, grid.breaks = TRUE,
            axis.labels = c("Feedstock = F","T = 35°C","Inoculum = A", "pH = 9","pH = 5"), 
            p.shape = TRUE,
            ci.lvl = 0.95,
            p.threshold = c(0.05, 0.01, 0.001),
            spacing = 0.7,
            dot.size = 2,legend.title = "Dependent Variables",
            colors = c("black", "blue", "darkgreen", "darkorange", "purple","gray", "red"),
            vline.color = "#061423")+ theme_pubr()+ labs_pubr()+
            theme(legend.position = "right",panel.grid.minor = element_line(color = "grey90"))+ 
            labs(y = "Coefficient Effect Size") 
plot_all_fe 

png("Final_Figures/fig_4_fe_random_coef.png", width = 6, height = 4, units = "in", res = 300)
print(plot_all_fe)
dev.off()




############################################################

############################################################
###  INTERACTION EFFECTS MODEL CONSTRUCTION ONLY AND STATS
############################################################
library(broom)
library(dplyr)
library(ggplot2)
library(gridExtra)

library(lme4)
library(emmeans)
library(lmerTest)
library(pbkrtest)
library(ggpubr)

library(performance)
library(gridExtra)
library(grid)
library(ggplot2)


# Define the response variables

data_acids_stats_end <- data_acids_stats %>% filter(between(day, 0.1, 12))
head(data_acids_stats_end)
data_acids_stats_end$pH <- as.factor(data_acids_stats_end$pH)
data_acids_stats_end$T <- as.factor(data_acids_stats_end$T)
data_acids_stats_end$inoc <- as.factor(data_acids_stats_end$inoc)
data_acids_stats_end$feedstock <- as.factor(data_acids_stats_end$feedstock)

data_acids_stats_end$pH <- relevel(data_acids_stats_end$pH, ref = "7")
data_acids_stats_end$inoc <- relevel(data_acids_stats_end$inoc, ref = "B")
data_acids_stats_end$T <- relevel(data_acids_stats_end$T, ref = "45")
data_acids_stats_end$feedstock <- relevel(data_acids_stats_end$feedstock, ref = "M")
#day is a continous variable
response_vars <- c("C2", "C3", "C4", "C5", "C6","C7", "total")

# Loop through each response variable
int_models <- list()
int_anovas <- list()
int_models_summary<- list()

l_emm_F_T <- list()
l_con_F_T <- list()
l_emm_F_pH<- list()
l_con_F_pH<- list()
l_emm_F_inoc <- list()
l_con_F_inoc <- list()
l_con_all_simple <- list()
l_con_all_expanded <- list()



# Apply effect coding
# contrasts(data_acids_stats_end$pH) <- contr.sum
# contrasts(data_acids_stats_end$inoc) <- contr.sum
# contrasts(data_acids_stats_end$T) <- contr.sum
# contrasts(data_acids_stats_end$feedstock) <- contr.sum


for (response in response_vars) {


  formula <- as.formula(paste(response, "~ pH * inoc * T * feedstock  + (1 | reactor)"))
  
  # Fit the model
  model <- lmer(formula, data = data_acids_stats_end)
  
  # Summarize the model
  summary_model <- summary(model)
  
  #statistical anaysis for each model
  anova_int_model <- anova(model)
  
  emm_F_T <- emmeans(model, ~ feedstock|T)
  contrast_F_T <- contrast(emm_F_T,
                        method = "pairwise",
                        simple = "each",
                        combine = TRUE,
                        adjust = "none") %>% summary(infer = TRUE)
  
  emm_F_pH <- emmeans(model, ~ feedstock|pH)
  contrast_F_pH <- contrast(emm_F_pH,
                                 method = "pairwise",
                                 simple = "each",
                                 combine = TRUE,
                                 adjust = "none") %>% summary(infer = TRUE)
  
  emm_F_inoc <- emmeans(model, ~ feedstock|inoc)
  contrast_F_inoc <- contrast(emm_F_inoc,
                                 method = "pairwise",
                                 simple = "each",
                                 combine = TRUE,
                                 adjust = "none") %>% summary(infer = TRUE)
  
  
  contrast_simple <- contrast(emmeans(model, ~ inoc*pH*T|feedstock),
                              method = "pairwise",
                              simple = "each",
                              combine = TRUE,
                              adjust = "none") %>% summary(infer = TRUE)
  
  contrast_expanded <- contrast(emmeans(model, ~ inoc*pH*T|feedstock),
                              method = "pairwise",
                              #simple = "each",
                              combine = TRUE,
                              adjust = "none") %>% summary(infer = TRUE)
  
  

  # Store the model and summary
  int_models[[response]] <- model
  int_models_summary[[response]] <- summary_model
  int_anovas[[response]] <- anova_int_model
  l_emm_F_T[[response]] <- emm_F_T
  l_con_F_T[[response]] <- contrast_F_T
  l_emm_F_pH[[response]] <- emm_F_pH
  l_con_F_pH[[response]] <- contrast_F_pH
  l_emm_F_inoc[[response]] <- emm_F_inoc
  l_con_F_inoc[[response]] <- contrast_F_inoc
  l_con_all_simple[[response]] <- contrast_simple
  l_con_all_expanded[[response]] <- contrast_expanded
}

#########
#### RESIDUAL PLOTS FIG S5
###########


# model <- nlmer(C2 ~ pH * inoc * T * feedstock + (1 | reactor), data = data_acids_stats_end,)

# Function to create a residual plot with title
create_residual_plot_with_title <- function(model, title) {
  residuals_plot <- resid_panel(model)
  grid.arrange(residuals_plot, top = textGrob(title, gp = gpar(fontsize = 14, fontface = "bold")))
}

# Create and capture residual plots for each model
iplot_C2 <- create_residual_plot_with_title(int_models$C2, "C2")
iplot_C3 <- create_residual_plot_with_title(int_models$C3, "C3")
iplot_C4 <- create_residual_plot_with_title(int_models$C4, "C4")
iplot_C5 <- create_residual_plot_with_title(int_models$C5, "C5")
iplot_C6 <- create_residual_plot_with_title(int_models$C6, "C6")
iplot_C7 <- create_residual_plot_with_title(int_models$C7, "C7")
iplot_total <- create_residual_plot_with_title(int_models$total, "Total")


png("Final_Figures/figS_X_resid_int_with_randomeffects_2.png", width = 10, height = 9, units = "in", res = 300)
resid_int_combined_plot <- grid.arrange(iplot_C2, iplot_C3, iplot_C4,iplot_C5, iplot_C6,iplot_C7, iplot_total,  ncol = 3)
dev.off()




############################################################



############################################################

############################################################
###   TABLES FOR INTERACTION EFFECTS TABLE S3_MODELSUM AND S4_ANOVA
############################################################
library(writexl)
library(dplyr)
library(writexl)
library(dplyr)

# Process the all_models_summary list
# int_models_summary_processed <- process_list(int_models_summary)

int_anovas
# Process the all_anovas list
int_anovas_processed <- process_list(int_anovas)


# Create a list of data frames to export to different sheets



int_anovas_processed

int_anovas_F=select(int_anovas_processed, term, model, statistic)
int_anovas_p=select(int_anovas_processed, term, model, p.value)
int_anovas_sig=select(int_anovas_processed, term, model, significance)

int_anovas_wide_F <- pivot_wider(int_anovas_F, names_from = model, values_from = statistic)
int_anovas_wide_p <- pivot_wider(int_anovas_p, names_from = model, values_from = p.value)
int_anovas_wide_sig <- pivot_wider(int_anovas_sig, names_from = model, values_from = significance)


int_tables_anovas <- list(
  "table_SX_int_ANOVAs" = int_anovas_processed,
  "F_stat_int_ANOVAs"=int_anovas_wide_F,
  "p_stat_int_ANOVAs"=int_anovas_wide_p,
  "sig_stat_int_ANOVAs"=int_anovas_wide_sig)
# Export to Excel file
write_xlsx(int_tables_anovas, path = "Final_Tables/table_SX_int_ANOVAs.xlsx") ## Figure 4 in paper will be from this table 






#print and paste in excel_see word file created manually: table_SX_model_summaries_fe_and_int.docx


tab_model(int_models$C2,int_models$C3,int_models$C4) ## Figure will be from this table table_SX_fe_model_summary
tab_model(int_models$C5,int_models$C6,int_models$C7,int_models$total)


#################################

###############################
#### STATISTICAL ANALYSIS INTERACTION EMMEANS TABLES ###
##############################
library(openxlsx)
library(dplyr)

# Create a new Excel workbook
wb <- createWorkbook()

# Function to combine contrast lists into a single data frame
combine_contrasts <- function(contrast_list) {
  bind_rows(lapply(names(contrast_list), function(response) {
    mutate(contrast_list[[response]], Acid = response)
  }))
}

# List of combined data frames and their corresponding sheet names
combined_contrasts <- list(
  F_T = combine_contrasts(l_con_F_T),
  F_pH = combine_contrasts(l_con_F_pH),
  F_inoc = combine_contrasts(l_con_F_inoc),
  All_Simple = combine_contrasts(l_con_all_simple),
  All_Expanded = combine_contrasts(l_con_all_expanded)
)

# Loop through each combined data frame and add it to the workbook
for (sheet_name in names(combined_contrasts)) {
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet = sheet_name, combined_contrasts[[sheet_name]])
}

# Save the workbook to a file
saveWorkbook(wb, "Final_Tables/table_SX_emm_contrast_results.xlsx", overwrite = TRUE)









#################################




#####################################
#### ADDITIONAL CODE NOT USE FROM HERE
#################################



########################################
#STACKED BARPLOTS FOR SUPPLEMENTAL
########################################
#########

library(dplyr) ##need to load again

data_filter_last <- all_data_acids_sum %>%
  filter(!grepl("total", acid)) %>%
  filter(between(day, 10.7, 13))  %>% filter(T!=55)

data_filter_last$inoc_T <- interaction(data_filter_last$inoc, data_filter_last$pH,sep = "°C_pH ")
data_filter_last$T_pH <- interaction(data_filter_last$T, data_filter_last$pH,sep = "°C_pH ")

# Plot
sbp_D4 = ggplot(data_filter_last, aes(x = inoc, y = conc_gL, fill = acid)) + 
  geom_bar(stat = "identity", position = "stack", colour = "black") +
  theme_minimal() + 
  facet_grid(feedstock ~ T_pH) +
  labs(x = "Inoculum", y = "Total FA [g/L]") +
  scale_fill_manual(values = c("red", "blue", "green", "orange", "purple", "yellow", "magenta", "darkgreen", "darkred")) + 
  theme_pubr() + labs_pubr()+theme(legend.position = "right")


png("Figures Manuscript/figS1_FA_stacked_day12.png", width = 7, height = 8, units = "in", res = 300)
print(sbp_D4)
dev.off()



#################################




############################################################
###  EFFECTS PLOTS FOR COLLECTIVE MODEL ON FIXED  PARAMETERS SUPPLEMENTAL FIGURE S4
############################################################

library(effects)



png("Figures Manuscript/FigS4_C2_fixed_effects.png", width = 6, height = 6, units = "in", res = 300)
plot(allEffects(fe_models$C2))
dev.off()
png("Figures Manuscript/FigS4_C3_fixed_effects.png", width = 6, height = 6, units = "in", res = 300)
plot(allEffects(fe_models$C3))
dev.off()
png("Figures Manuscript/FigS4_C4_fixed_effects.png", width = 6, height = 6, units = "in", res = 300)
plot(allEffects(fe_models$C4))
dev.off()
png("Figures Manuscript/FigS4_C5_fixed_effects.png", width = 6, height = 6, units = "in", res = 300)
plot(allEffects(fe_models$C5))
dev.off()
png("Figures Manuscript/FigS4_C6_fixed_effects.png", width = 6, height = 6, units = "in", res = 300)
plot(allEffects(fe_models$C6))
dev.off()
png("Figures Manuscript/FigS4_total_fixed_effects.png", width = 6, height = 6, units = "in", res = 300)
plot(allEffects(fe_models$total))
dev.off()


# Plot the effects with customized colors for lines and error bars



############################################################


############################################################
###   FIXED EFFECCT BOXPLOTS FIG S4B 
############################################################

my_comparisons <- list( c("A","B"),c("5","7"),c("5","9"),c("7","9"))

data_acids_stats_end$pH_inoc<- as.factor(data_acids_stats_end$pH_inoc)
data_acids_stats_end$T_pH<- as.factor(data_acids_stats_end$T_pH)
data_acids_stats_end$feedstock_T<- as.factor(data_acids_stats_end$feedstock_T)


str(filtered_sd)

filtered_sd$pH <- as.factor(filtered_sd$pH)
filtered_sd$T <- as.factor(filtered_sd$T)
filtered_sd$inoc <- as.factor(filtered_sd$inoc)
filtered_sd$feedstock <- as.factor(filtered_sd$feedstock)
filtered_sd$acid <- as.factor(filtered_sd$acid)


filtered_sd_stacked_end <- filtered_sd %>% filter(between(day, 4.8, 12))
#response_vars <- c("C2", "C3", "C4", "C5", "C6", "total")


library(ggplot2)
library(ggpubr)

my_comparisons <- list(c(1, 2))


boxplot_fe_inoc<-ggboxplot(filtered_sd_stacked_end, x = "inoc", y = "conc_gL",
                           color = "inoc",
                           shape = "inoc",
                           scales = "free_y", 
                           palette = c("darkgreen","#E69F00"), add = "jitter", facet.by = c("acid"), 
                           add.params = list(size = 0.4)) + 
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", size = 3) + 
  #stat_compare_means(method = "wilcox.test", label.y = max(filtered_sd_stacked_end$conc_gL), size = 2.5) +
  ylab("Concentration (g/L)") +
  xlab("Inoculum") + 
  labs(color = "Inoculum", shape = "Inoculum") +
  theme_pubr() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 10, face = "bold"),legend.position = "none")

boxplot_fe_T<-ggboxplot(filtered_sd_stacked_end, x = "T", y = "conc_gL",
                        color = "T",
                        shape = "T",
                        scales = "free_y", 
                        palette = c("black","purple"), add = "jitter", facet.by = c("acid"), 
                        add.params = list(size = 0.4)) + 
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", size = 3) + 
  #stat_compare_means(method = "wilcox.test", label.y = max(filtered_sd_stacked_end$conc_gL), size = 2.5) +
  ylab("Concentration (g/L)") +
  xlab("Temperature") + 
  labs(color = "Temperature", shape = "Temperature") +
  theme_pubr() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 10, face = "bold"),legend.position = "none")



boxplot_fe_F<-ggboxplot(filtered_sd_stacked_end, x = "feedstock", y = "conc_gL",
                        color = "feedstock",
                        shape = "feedstock",
                        scales = "free_y", 
                        palette = c("blue","red"), add = "jitter", facet.by = c("acid"), 
                        add.params = list(size = 0.4)) + 
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", size = 3) + 
  #stat_compare_means(method = "wilcox.test", label.y = max(filtered_sd_stacked_end$conc_gL), size = 2.5) +
  ylab("Concentration (g/L)") +
  xlab("Feedstock") + 
  labs(color = "Feedstock", shape = "Feedstock") +
  theme_pubr() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 10, face = "bold"),legend.position = "none")

my_comparisons <- list(c(1, 2),c(2, 3),c(1, 3))


boxplot_fe_pH<-ggboxplot(filtered_sd_stacked_end, x = "pH", y = "conc_gL",
                         color = "pH",
                         shape = "pH",
                         scales = "free_y", 
                         palette = c("#A73030FF","#8F7700FF","#003C67FF"), add = "jitter", facet.by = c("acid"), 
                         add.params = list(size = 0.4)) + 
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", size = 3) + 
  #stat_compare_means(size = 2) +
  ylab("Concentration (g/L)") +
  xlab("pH") + 
  labs(color = "pH", shape = "pH") +
  theme_pubr() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 10, face = "bold"), legend.position = "none")


grid.arrange(arrangeGrob(boxplot_fe_pH,boxplot_fe_inoc, boxplot_fe_T, boxplot_fe_F, ncol = 2))


png("Figures Manuscript/FigS4b_boxplots_fixed_effects.png", width = 8, height = 12, units = "in", res = 300)
grid.arrange(arrangeGrob(boxplot_fe_pH,boxplot_fe_inoc, boxplot_fe_T, boxplot_fe_F, ncol = 2))
dev.off()




############################################################



############################################################
###   INTERACTION PLOTS 2-WAY
############################################################


# Assuming int_models is a list of interaction models and titles is a list of titles for the plots

F_inoc_plots <- list()
F_pH_plots <- list()
F_T_plots <- list()


pal_okabe_ito <- c(
  "#E69F00",
  "#56B4E9",
  "#009E73",
  "#F0E442",
  "#0072B2",
  "#D55E00",
  "#CC79A7"
)

response_vars <- c("C2", "C3", "C4", "C5", "C6", "total")

for (response in response_vars) {
  
  plot_F_inoc <- ggplot_the_model(
    fit = int_models[[response]],  # Assuming int_models is a list or named vector
    fit_emm = l_emm_F_inoc[[response]],  # Assuming l_emm_F_inoc is a list or named vector
    fit_pairs = l_con_F_inoc[[response]],  # Assuming l_con_F_inoc is a list or named vector
    palette = c("darkgreen","#E69F00"),
    y_label = bquote(.(response)~"(g/L)"),  # Corrected line
    g_label = "none",
    effect_label = "inoculum:feedstock effects")
  
  plot_F_pH <- ggplot_the_model(
    fit = int_models[[response]],  # Assuming int_models is a list or named vector
    fit_emm = l_emm_F_pH[[response]],  # Assuming l_emm_F_inoc is a list or named vector
    fit_pairs = l_con_F_pH[[response]],  # Assuming l_con_F_inoc is a list or named vector
    palette = c("#D55E10","#EFC000FF","#0072B2"),
    y_label = bquote(.(response)~"(g/L)"),  # Corrected line
    g_label = "none",
    effect_label = "pH:feedstock effects")
  
  plot_F_T <- ggplot_the_model(
    fit = int_models[[response]],  # Assuming int_models is a list or named vector
    fit_emm = l_emm_F_T[[response]],  # Assuming l_emm_F_inoc is a list or named vector
    fit_pairs = l_con_F_T[[response]],  # Assuming l_con_F_inoc is a list or named vector
    palette = c("black","purple"),
    y_label = bquote(.(response)~"(g/L)"),  # Corrected line
    g_label = "none",
    effect_label = "T:feedstock effects")
  
  F_inoc_plots[[response]] <- plot_F_inoc
  F_pH_plots[[response]] <- plot_F_pH
  F_T_plots[[response]] <- plot_F_T
  
}
#plot_F_T



# Load the gridExtra package
library(gridExtra)



png("Figures Manuscript/Fig4a_F_pH_int.png", width = 9, height = 9, units = "in", res = 300)
grid.arrange(grobs = F_pH_plots, ncol = 3)
dev.off()

png("Figures Manuscript/Fig4b_F_T_int.png", width = 9, height = 9, units = "in", res = 300)
grid.arrange(grobs = F_T_plots, ncol = 3)
dev.off()

png("Figures Manuscript/Fig4c_F_inoc_int.png", width = 9, height = 9, units = "in", res = 300)
grid.arrange(grobs = F_inoc_plots, ncol = 3)
dev.off()


# Assuming your plots are stored in F_inoc_plots, F_pH_plots, and F_T_plots

# # Create a list of plots in the order you want to display them
# F_pH_plots <- list()
# 
# 
# response_vars <- c("C2", "C3", "C4", "C5", "C6", "total")
# 
# 
# for (response in response_vars) {
#   # all_plots <- c(all_plots, list(F_inoc_plots[[response]], F_pH_plots[[response]], F_T_plots[[response]]))
#   F_pH_plots <- c(all_plots, list(F_pH_plots[[response]]))
#   F_T_plots <- c(all_plots, list(F_T_plots[[response]]))
#   F_inoc_plots <- c(all_plots, list(F_inoc_plots[[response]]))
# }




emm_F_ipT <- emmeans(model, ~ inoc*pH*T|feedstock)

contrast_ex <- contrast(emm_F_ipT,
                        method = "pairwise",
                        simple = "each",
                        combine = TRUE,
                        adjust = "none") %>% summary(infer = TRUE)


contrast_all <- contrast(emm_F_ipT,
                         method = "pairwise",
                         #simple = "each",
                         combine = TRUE,
                         adjust = "none") %>% summary(infer = TRUE)




plot(emm_F_ipT, comparisons = TRUE)

contrast_F_ipT<- contrast(emm_F_ipT,
                          method = "pairwise",
                          
                          combine = TRUE,
                          adjust = "none") %>% summary(infer = TRUE)


# Plot the pairs with the specified model
p <- plot(pairs(emmeans(int_models$C2, ~T|feedstock|pH), adjust="none"))

# Add vertical line, theme, and labs, while removing facet labels
p <- p + 
  geom_vline(xintercept = 0, size = 1, color = "red") + 
  theme_pubr(base_size = 7) +
  labs_pubr() +
  #theme(strip.text = element_blank()) +  # Remove facet labels
  ggtitle("C2") +
  theme(
    legend.position = "right",
    panel.grid.minor = element_line(color = "grey90"),
    panel.grid.major = element_line(color = "grey90"))

# Display the plot
print(p)

plot_C2



############################################################


############################################################
###   INTERACTION EFFECCT PLOTS FIG S7 3 AND 4 - WAY FOG 5 FIGS7
############################################################

library(ggplot2)
library(gridExtra)
library(emmeans)
library(ggpubr)

# List of models and response variables
models <- int_models[c("C2", "C3", "C4", "C5", "C6", "total")]
response_titles <- response_vars


# Function to create a plot with customized colors and a title
create_custom_plot <- function(model, formula, colors, title) {
  emmeans::emmip(model, formula, CIs = TRUE, dotarg = list()) +
    theme_pubr(legend = c("right")) +
    labs_pubr() +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors) +
    ggtitle(title) + 
    theme(plot.title = element_text(hjust = 0.5))
}

# Define the formulas and color schemes for each type of plot
plot_settings <- list(
  list(formula = T ~ pH | feedstock, colors = c("black", "purple")),
  list(formula = T ~ pH |inoc| feedstock, colors = c("black", "purple")),
  list(formula = inoc ~ pH | feedstock, colors = c("darkgreen", "#E69F00")),
  list(formula = inoc ~ pH |T| feedstock, colors = c("darkgreen", "#E69F00"))

  )

# Generate the plots
plots_TpHF <- list()
plots_TpHiF <- list()

plots_ipHF <- list()
plots_ipHTF <- list()


for (i in seq_along(models)) {
  model <- models[[i]]
  title <- response_titles[i]
  
  plots_TpHF[[i]] <- create_custom_plot(model, plot_settings[[1]]$formula, plot_settings[[1]]$colors, title)
  plots_TpHiF[[i]] <- create_custom_plot(model, plot_settings[[2]]$formula, plot_settings[[2]]$colors, title)
  plots_ipHF[[i]] <- create_custom_plot(model, plot_settings[[3]]$formula, plot_settings[[3]]$colors, title)
  plots_ipHTF[[i]] <- create_custom_plot(model, plot_settings[[4]]$formula, plot_settings[[4]]$colors, title)
  
  }

# Combine and display the plots in grids


ipHF_combined_plot_all <- grid.arrange(grobs = plots_ipHF, ncol = 2)
ipHTF_combined_plot_all <- grid.arrange(grobs = plots_ipHTF, ncol = 2)
TpHF_combined_plot_all <- grid.arrange(grobs = plots_TpHF, ncol = 2)
TpHiF_combined_plot_all <- grid.arrange(grobs = plots_TpHiF, ncol = 2)


ipHF_TpHF_combined_plot_all=grid.arrange(grobs = c(rbind(plots_ipHF, plots_TpHF)), ncol = 4)




# To save the combined plot

ggsave("Figures Manuscript/FigS7a_plots_ipHTF.png.png", ipHTF_combined_plot_all, width =8, height = 12)
ggsave("Figures Manuscript/FigS7b_plots_TpHiF.png.png", TpHiF_combined_plot_all, width = 8, height = 12)
ggsave("Figures Manuscript/Fig5_plots_ipHF_TpHF_3way.png.png", ipHF_TpHF_combined_plot_all, width = 16, height = 8)






# 
# 
# 
# emmeans::emmip(int_models$C2, inoc ~ pH ~ T | feedstock, CIs = TRUE,dotarg = list()) +
#   theme_pubr(legend = c("right")) +
#   labs_pubr() +
#   scale_color_manual(values = c("#A73030FF","#CD534CFF","#8F7700FF", "#EFC000FF","#003C67FF","#0073C2FF")) +
#   scale_fill_manual(values = c("#A73030FF","#CD534CFF","#8F7700FF", "#EFC000FF","#003C67FF","#0073C2FF")) +
#   ggtitle("title") + 
#   theme(plot.title = element_text(hjust = 0.5))
# 
# emmeans::emmip(int_models$C2, inoc ~ pH | T | feedstock, CIs = TRUE,dotarg = list()) +
#   theme_pubr(legend = c("right")) +
#   labs_pubr() +
#   scale_color_manual(values = c("#A73030FF","#CD534CFF","#8F7700FF", "#EFC000FF","#003C67FF","#0073C2FF")) +
#   scale_fill_manual(values = c("#A73030FF","#CD534CFF","#8F7700FF", "#EFC000FF","#003C67FF","#0073C2FF")) +
#   ggtitle("title") + 
#   theme(plot.title = element_text(hjust = 0.5))
# 
# 
# 
# # Function to create a plot with customized colors and a title
# create_custom_plot <- function(model, title) {
#   emmeans::emmip(model, inoc ~ pH ~ T | feedstock, CIs = TRUE,dotarg = list()) +
#     theme_pubr(legend = c("right")) +
#     labs_pubr() +
#     scale_color_manual(values = c("#A73030FF","#CD534CFF","#8F7700FF", "#EFC000FF","#003C67FF","#0073C2FF")) +
#     scale_fill_manual(values = c("#A73030FF","#CD534CFF","#8F7700FF", "#EFC000FF","#003C67FF","#0073C2FF")) +
#     ggtitle(title) + 
#     theme(plot.title = element_text(hjust = 0.5))
# }




















############################################################
###  EXTRA CODE NOT USED IN MAIN
############################################################


############################################################
###   DATA DISTRIBUTION PLOT
############################################################

# Load necessary libraries
library(dplyr)

# List of acids to log-transform
acids <- c("C2", "C3", "C4", "C5", "C6", "ISOC4", "ISOC5", "C7", "ISOC6", "total")

# Function to create a distribution plot with title
create_distribution_plot_with_title <- function(data, acid, title) {
  p <- ggplot(data, aes_string(x = acid)) +
    geom_histogram(aes(y = ..density..), bins = 30, fill = "blue", alpha = 0.7) +
    geom_density(color = "red", size = 1) +
    theme_minimal() +
    labs(title = title, x = acid, y = "Density")
  p
}




# Create distribution plots for each acid
plots <- lapply(acids, function(acid) {
  create_distribution_plot_with_title(data_acids_stats_end, acid, acid)
})

# Combine all plots into a single grid
combined_plot <- grid.arrange(grobs = plots, ncol = 2)

# Display the combined plot


png("Figures Manuscript/FigS6_data_distribution_plot.png", width = 9, height = 9, units = "in", res = 300)
grid.newpage()
grid.draw(combined_plot)
dev.off()



############################################################
###  FIXED EFFECT MODEL CONSTRUCTION FIG S2 AND S3 NO RANDOM EFFECTS
############################################################
# 
# #List model plots
# 
# pH_models <- list()
# pH_anovas <- list()
# pH_plots_residuals <- list()
# 
# T_models <- list()
# T_anovas <- list()
# T_plots_residuals <- list()
# 
# F_models <- list()
# F_anovas <- list()
# F_plots_residuals <- list()
# 
# I_models <- list()
# I_anovas <- list()
# I_plots_residuals <- list()
# 
# all_models <- list()
# all_anovas <- list()
# all_plots_residuals <- list()
# all_models_summary<- list()
# all_plots_qq <- list()
# 
# # Define the response variables
# data_acids_stats_end <- data_acids_stats %>% filter(between(day, 10, 12))
# data_acids_stats_end$pH <- relevel(data_acids_stats_end$pH, ref = "7")
# data_acids_stats_end$inoc <- relevel(data_acids_stats_end$inoc, ref = "B")
# data_acids_stats_end$T <- relevel(data_acids_stats_end$T, ref = "45")
# data_acids_stats_end$feedstock <- relevel(data_acids_stats_end$feedstock, ref = "M")
# response_vars <- c("C2", "C3", "C4", "C5", "C6", "total")
# 
# ### pH MODELS
# # Loop through each response variable
# for (response in response_vars) {
#   # Create the formula
#   formula <- as.formula(paste(response, "~ pH"))
#   
#   # Fit the model
#   model <- lm(formula, data = data_acids_stats_end)
#   anova_lm <- anova(model)
#   pH_models[[response]] <- model
#   pH_anovas[[response]] <- anova_lm
#   
#   # Get augmented data with fitted values and residuals
#   aug_data <- augment(model)
#   
#   # Create diagnostic plots manually
#   residuals_plot <- ggplot(aug_data, aes(.fitted, .resid)) +
#     geom_point() +
#     geom_smooth(method = "lm", se = FALSE, col = "red") +
#     labs(
#       title = paste(response, "Residuals vs Fitted"),
#       subtitle = paste("Model:", deparse(formula)),
#       x = "Fitted values", 
#       y = "Residuals"
#     ) +
#     theme_minimal()
#   
#   # Store the plots in lists
#   pH_plots_residuals[[response]] <- residuals_plot
# }
# 
# ### TEMPERATURE MODELS
# # Loop through each response variable
# for (response in response_vars) {
#   # Create the formula
#   formula <- as.formula(paste(response, "~ T"))
#   
#   # Fit the model
#   model <- lm(formula, data = data_acids_stats_end)
#   anova_lm <- anova(model)
#   T_models[[response]] <- model
#   T_anovas[[response]] <- anova_lm
#   
#   # Get augmented data with fitted values and residuals
#   aug_data <- augment(model)
#   
#   # Create diagnostic plots manually
#   residuals_plot <- ggplot(aug_data, aes(.fitted, .resid)) +
#     geom_point() +
#     geom_smooth(method = "lm", se = FALSE, col = "red") +
#     labs(
#       title = paste(response, "Residuals vs Fitted"),
#       subtitle = paste("Model:", deparse(formula)),
#       x = "Fitted values", 
#       y = "Residuals"
#     ) +
#     theme_minimal()
#   
#   
#   # Store the plots in lists
#   T_plots_residuals[[response]] <- residuals_plot
# }
# 
# ### FEEDSTOCK MODELS
# # Loop through each response variable
# for (response in response_vars) {
#   # Create the formula
#   formula <- as.formula(paste(response, "~ feedstock"))
#   
#   # Fit the model
#   model <- lm(formula, data = data_acids_stats_end)
#   anova_lm <- anova(model)
#   F_models[[response]] <- model
#   F_anovas[[response]] <- anova_lm
#   
#   # Get augmented data with fitted values and residuals
#   aug_data <- augment(model)
#   
#   # Create diagnostic plots manually
#   residuals_plot <- ggplot(aug_data, aes(.fitted, .resid)) +
#     geom_point() +
#     geom_smooth(method = "lm", se = FALSE, col = "red") +
#     labs(
#       title = paste(response, "Residuals vs Fitted"),
#       subtitle = paste("Model:", deparse(formula)),
#       x = "Fitted values", 
#       y = "Residuals"
#     ) +
#     theme_minimal()
#   
#   
#   # Store the plots in lists
#   F_plots_residuals[[response]] <- residuals_plot
# }
# 
# ### INOCULUM MODELS 
# # Loop through each response variable
# for (response in response_vars) {
#   # Create the formula
#   formula <- as.formula(paste(response, "~ inoc"))
#   
#   # Fit the model
#   model <- lm(formula, data = data_acids_stats_end)
#   anova_lm <- anova(model)
#   I_models[[response]] <- model
#   I_anovas[[response]] <- anova_lm
#   
#   # Get augmented data with fitted values and residuals
#   aug_data <- augment(model)
#   
#   # Create diagnostic plots manually
#   residuals_plot <- ggplot(aug_data, aes(.fitted, .resid)) +
#     geom_point() +
#     geom_smooth(method = "lm", se = FALSE, col = "red") +
#     labs(
#       title = paste(response, "Residuals vs Fitted"),
#       subtitle = paste("Model:", deparse(formula)),
#       x = "Fitted values", 
#       y = "Residuals"
#     ) +
#     theme_minimal()
# 
#   
#   # Store the plots in lists
#   I_plots_residuals[[response]] <- residuals_plot
# }
# 
# ### ALL FIXED EFFECTS
# # Loop through each response variable
# for (response in response_vars) {
#   # Create the formula
#   formula <- as.formula(paste(response, "~ pH+T+feedstock+inoc"))
#   
#   # Fit the model
#   model <- lm(formula, data = data_acids_stats_end)
#   anova_lm <- anova(model)
#   summary_model<- summary(model)
#   all_models[[response]] <- model
#   all_anovas[[response]] <- anova_lm
#   all_models_summary[[response]] <- summary_model
#   
#   # Get augmented data with fitted values and residuals
#   aug_data <- augment(model)
#   
#   # Create diagnostic plots manually
#   residuals_plot <- ggplot(aug_data, aes(.fitted, .resid)) +
#     geom_point() +
#     geom_smooth(method = "lm", se = FALSE, col = "red") +
#     labs(
#       title = paste(response, "Residuals vs Fitted"),
#       subtitle = paste("Model:", deparse(formula)),
#       x = "Fitted values", 
#       y = "Residuals"
#     ) +
#     theme_minimal()
#   
#   qq_plot <- ggplot(aug_data, aes(sample = .std.resid)) +
#     stat_qq() +
#     stat_qq_line(col = "red") +
#     labs(
#       title = paste(response, "Normal Q-Q"),
#       subtitle = paste("Model:", deparse(formula)),
#       x = "Theoretical Quantiles", 
#       y = "Standardized Residuals"
#     ) +
#     theme_minimal()
#   
#   # Store the plots in lists
#   F_plots_residuals[[response]] <- residuals_plot
#   F_plots_qq[[response]] <- qq_plot
#   
#   
#   # Store the plots in lists
#   all_plots_residuals[[response]] <- residuals_plot
#   all_plots_qq[[response]] <- qq_plot
#   
# }
# 



#################################


############################################################
###  FIXED EFFECT MODEL RESIDUAL PLOTS FIGURE S2 NO RANDOM EFFECTS
############################################################


# 
# 
# png("Figures Manuscript/FigS2a_residualplot_all.png", width = 7, height = 7, units = "in", res = 300)
# all_all_plots_residuals <- do.call(grid.arrange, c(all_plots_residuals, ncol = 2))
# dev.off()
# 
# png("Figures Manuscript/FigS2a_qq_all.png", width = 7, height = 7, units = "in", res = 300)
# all_all_plots_qq <- do.call(grid.arrange, c(all_all_plots_qq, ncol = 2))
# dev.off()
# 
# png("Figures Manuscript/FigS2b_residualplot_pH.png", width = 7, height = 7, units = "in", res = 300)
# pH_all_plots_residuals <- do.call(grid.arrange, c(pH_plots_residuals, ncol = 2))
# dev.off()
# 
# png("Figures Manuscript/FigS2c_residualplot_T.png", width = 7, height = 7, units = "in", res = 300)
# T_all_plots_residuals <- do.call(grid.arrange, c(T_plots_residuals, ncol = 2))
# dev.off()
# 
# png("Figures Manuscript/FigS2d_residualplot_F.png", width = 7, height = 7, units = "in", res = 300)
# F_all_plots_residuals <- do.call(grid.arrange, c(F_plots_residuals, ncol = 2))
# dev.off()
# 
# png("Figures Manuscript/FigS2e_residualplot_I.png", width = 7, height = 7, units = "in", res = 300)
# I_all_plots_residuals <- do.call(grid.arrange, c(I_plots_residuals, ncol = 2))
# dev.off()






############################################################


############################################################
###  EFFECTS PLOTS FOR COLLECTIVE MODEL ON FIXED  PARAMETERS SUPPLEMENTAL FIGURE S4
############################################################

# library(effects)


# 
# plot(allEffects(fe_models$C2))
# allEffects(fe_models$C2)
# 
# plot(tidy_all_effects(fe_models$C2))
# 

# # facet by group, use pre-defined color palette
# dat <- predict_response(fe_models$C2, terms = c("pH"))
# plot(dat, colors = "hero",show_ci = TRUE)
# library(jtools)
# 
# effect_plot(fe_models$C2, pred = pH,jitter = .2, interval = TRUE, y.label = "Concentration",
#             cat.geom = "line",plot.points = TRUE)
# 
# 
# 
# library(effects)
# library(car)
# library(MASS)
# library(splines)
# library(lattice)
# 
# # plot(Effect(focal.predictors = c("pH"), 
# #             mod = fe_models$C2, latent = TRUE),
# #      rug = FALSE, axes = list(grid = TRUE), multiline=TRUE, 
# #      colors = c("red", "blue", "green"), 
# #      lattice = list(layout = c(5,1), key.args = list(space="top")))
# 
# 
# png("Figures Manuscript/FigS4_C2_fixed_effects.png", width = 6, height = 6, units = "in", res = 300)
# plot(allEffects(all_models$C2))
# dev.off()
# png("Figures Manuscript/FigS4_C3_fixed_effects.png", width = 6, height = 6, units = "in", res = 300)
# plot(allEffects(all_models$C3))
# dev.off()
# png("Figures Manuscript/FigS4_C4_fixed_effects.png", width = 6, height = 6, units = "in", res = 300)
# plot(allEffects(all_models$C4))
# dev.off()
# png("Figures Manuscript/FigS4_C5_fixed_effects.png", width = 6, height = 6, units = "in", res = 300)
# plot(allEffects(all_models$C5))
# dev.off()
# png("Figures Manuscript/FigS4_C6_fixed_effects.png", width = 6, height = 6, units = "in", res = 300)
# plot(allEffects(all_models$C6))
# dev.off()
# png("Figures Manuscript/FigS4_total_fixed_effects.png", width = 6, height = 6, units = "in", res = 300)
# plot(allEffects(all_models$total))
# dev.off()



# plot_pH_fe<-plot_models(pH_models$C2, pH_models$C3, pH_models$C4,
#                         pH_models$C5, pH_models$C6, pH_models$total, grid.breaks = TRUE, 
#                         p.shape = TRUE,
#                         title="pH effect",
#                         ci.lvl = 0.95,
#                         p.threshold = c(0.05, 0.01, 0.001),
#                         spacing = 0.7,
#                         dot.size = 2,legend.title = "Dependent Variables",
#                         colors = c("black", "blue", "darkgreen", "darkorange", "purple", "red"),
#                         vline.color = "#061423")+ theme_pubr()+ labs_pubr()+
#   theme(legend.position = "right",panel.grid.minor = element_line(color = "grey90"))+small_text_theme
# 
# 
# plot_T_fe<-plot_models(T_models$C2, T_models$C3, T_models$C4,
#                         T_models$C5, T_models$C6, T_models$total, grid.breaks = TRUE, 
#                         p.shape = TRUE,
#                         title="Temperature effect",
#                         ci.lvl = 0.95,
#                         p.threshold = c(0.05, 0.01, 0.001),
#                         spacing = 0.7,
#                         dot.size = 2,legend.title = "Dependent Variables",
#                         colors = c("black", "blue", "darkgreen", "darkorange", "purple", "red"),
#                         vline.color = "#061423")+ theme_pubr()+ labs_pubr()+
#   theme(legend.position = "right",panel.grid.minor = element_line(color = "grey90"))+small_text_theme
# 
# plot_F_fe<-plot_models(F_models$C2, F_models$C3, F_models$C4,
#                        F_models$C5, F_models$C6, F_models$total, grid.breaks = TRUE, 
#                        p.shape = TRUE,
#                        title="Feedstock effect",
#                        ci.lvl = 0.95,
#                        p.threshold = c(0.05, 0.01, 0.001),
#                        spacing = 0.7,
#                        dot.size = 2,legend.title = "Dependent Variables",
#                        colors = c("black", "blue", "darkgreen", "darkorange", "purple", "red"),
#                        vline.color = "#061423")+ theme_pubr()+ labs_pubr()+
#   theme(legend.position = "right",panel.grid.minor = element_line(color = "grey90"))+small_text_theme
# 
# 
# plot_I_fe<-plot_models(I_models$C2, I_models$C3, I_models$C4,
#                        I_models$C5, I_models$C6, I_models$total, grid.breaks = TRUE, 
#                        p.shape = TRUE,
#                        title="Inoculum effect",
#                        ci.lvl = 0.95,
#                        p.threshold = c(0.05, 0.01, 0.001),
#                        spacing = 0.7,
#                        dot.size = 2,legend.title = "Dependent Variables",
#                        colors = c("black", "blue", "darkgreen", "darkorange", "purple", "red"),
#                        vline.color = "#061423")+ theme_pubr()+ labs_pubr()+
#   theme(legend.position = "right",panel.grid.minor = element_line(color = "grey90"))+small_text_theme
# 


#NO RANDOM EFFECTS
# png("Figures Manuscript/FigS3_individual_fixed_effects.png", width = 7, height = 7, units = "in", res = 300)
# grid.arrange(arrangeGrob(plot_pH_fe, plot_T_fe, plot_F_fe, plot_I_fe, ncol = 2))
# dev.off()

############################################################


plot_all_int_mod <- plot_models( int_models$total, int_models$C6, int_models$C5, int_models$C4, 
                                 int_models$C3, int_models$C2,
                                 grid.breaks = TRUE,
                                 # axis.labels = c("Feedstock = F","T = 35°C","Inoculum = A", "pH = 9","pH = 5"), 
                                 p.shape = TRUE,
                                 grid = TRUE,
                                 ci.lvl = 0.95,
                                 p.threshold = c(0.05, 0.01, 0.001),
                                 spacing = 1, # Increased spacing value
                                 dot.size = 2,
                                 legend.title = "Dependent Variables",
                                 colors = c("red", "purple", "darkorange","darkgreen",  "blue", "black"),
                                 vline.color = "#061423") + 
  theme_pubr() + 
  labs_pubr() +
  theme(legend.position = "right",
        # panel.grid.minor = element_line(color = "grey90"),
        panel.grid.major.y = element_line(color = "grey95"))

print(plot_all_int_mod)

png("Figures Manuscript/figS6a_interaction_effects_models_coef.png", width = 12, height = 7, units = "in", res = 300)
print(plot_all_int_mod)
dev.off()


plot_all_int_mod



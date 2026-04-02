### ************************************
### Data Processing ----
### ************************************

# function takes an input .fasta file in two-line format which looks like:
# line 1: >header
# line 2: DNA sequence
# ... and creates a two column dataframe




fasta_to_df <- function(fasta_file, hdr_colname = "", seq_colname = "") {
  # read in the .fasta file line by line
  fasta <- readLines(fasta_file)
  # identify header lines by finding line numbers with the '>' character
  hlines <- grep(pattern = ">", x = fasta)
  # create a three column df containing:
  # header line numbers (hdr)
  # sequence line number beginning (beg)
  # sequence line number end (end)
  slines <- data.frame(hdr = hlines,
                       beg = hlines+1,
                       end = c((hlines-1)[-1], length(fasta)))
  # create a vector of identical length to hdr_lines to be used for
  # storing values obtained in the loop below
  storage_vec <- rep(NA, length(hlines))
  # loop to obtain sequences
  for(i in 1:length(hlines)) {
    storage_vec[i] <- paste(fasta[slines$beg[i]:slines$end[i]], collapse = "")
  }
  # create a two column df containing:
  # the header with '>' character replaced (hdr)
  # representative sequence associated with the header (seq)
  new_df_v0 <- data.frame("hdr" = gsub(pattern = ">", replacement = "",
                                       x = fasta[hlines]),
                          "seq" = storage_vec)
  # replace generic column names with names specified in the function input
  new_df_v1 <- new_df_v0
  hdr_colnum <- which(names(new_df_v1) == "hdr")
  seq_colnum <- which(names(new_df_v1) == "seq")
  names(new_df_v1)[hdr_colnum] <- hdr_colname
  names(new_df_v1)[seq_colnum] <- seq_colname
  return(new_df_v1)
}
#
# example usage:
#df_fasta <- fasta_to_df(fasta_file = ifp_fasta,
#                        hdr_colname = "FeatureID",
#                        seq_colname = "RepSeq")

## *********************************** ## *********************************** ##
## *********************************** ## *********************************** ##

# trunc_tax() function requires an input df with a column containing either:
# greengenes formatted taxonomic strings, e.g.:
# k__; p__ ; c__; o__; f__; g__; s__"
# or silva taxonomic strings, e.g.:
# "D_0__;D_1__;D_2__;D_3__;D_4__;D_5__;D_6__"
# the output df will contain a new column with:
# the desired taxonomic level for the specified database naming convention
# NOTE: if the above strings are not in consecutive order...
# ... the output may not be correct i.e.:
# if string input is missing D_2__ = "D_0__;D_1__Example;D_3__Example;D_4__"
# output for ;D_1__ = "Example;D_3__Example;D_4__" rather than just "Example"
# output for ;D_2__ = "Unassigned"
# output for ;D_3__ = "Example"
trunc_tax <- function(data = data.frame, tax_type = c("Greengenes", "SILVA"),
                      tax_lvl = c("Kingdom", "Phylum", "Class", "Order",
                                  "Family", "Genus", "Species"),
                      input_col = "", str1 = "", str2 = "", classifier = "") {
  
  # internal checks to ensure correct input classes
  if (!inherits(data, "data.frame")) {
    stop("input data must be class 'data.frame'")
  }
  if (!tax_type == "Greengenes" && !tax_type == "SILVA") {
    stop("tax_type should be either 'Greengenes' or 'SILVA'")
  }
  if (!tax_lvl == "Kingdom" && !tax_lvl == "Phylum" && !tax_lvl == "Class" &&
      !tax_lvl == "Order" && !tax_lvl == "Family" && !tax_lvl == "Genus" &&
      !tax_lvl == "Species") {
    stop("tax_lvl should be one of: 'Kingdom', 'Phylum', 'Class', 'Order',
         'Family', 'Genus', 'Species'")
  }
  if (!inherits(input_col, "character")) {
    stop("input_col must be character")
  }
  if (!inherits(str1, "character")) {
    stop("str1 must be character")
  }
  if (!inherits(str2, "character")) {
    stop("str2 must be character")
  }
  if (!inherits(classifier, "character")) {
    stop("classifier must be character")
  }
  
  # if tax_lvl is Species, pasting the string split is handled differently
  num <- ifelse(tax_lvl == "Species", yes = 2, no = 1)
  
  # store original input df and format the new df
  new_df1 <- data
  new_df1$new_col <- 1
  
  # loop through data and truncate the taxonomic strings
  for(row in 1:nrow(new_df1)) {
    # split the string in the specified input column with input of str1
    splt1 <- strsplit(new_df1[row, input_col], str1)[[1]][2]
    # goal of first ifelse test:
    # yes = tax level was not assigned
    # no = tax level was likely assigned, split the string with input of str2
    splt2 <- ifelse(is.na(splt1), yes = "Unassigned",
                    no = paste("", strsplit(splt1, str2)[[1]][num], sep = ""))
    # goal of second ifelse test:
    # yes = tax level was not assigned
    # no = tax level was actually assigned
    new_df1[row, "new_col"] <- ifelse(nchar(splt2) == 0 | "NA" %in% splt2,
                                      yes = "Unassigned", no = splt2)
    # if tax_type is Greengenes, tax_lvl is Species and row is not Unassigned,
    # append Genus in front of Species assignment
    if (tax_type == "Greengenes" & tax_lvl == "Species" &
        !new_df1[row, "new_col"] == "Unassigned") {
      new_df1[row, "new_col"] <- paste("", strsplit(splt1, str2)[[1]][1],
                                       strsplit(splt1, str2)[[1]][2], sep = " ")
    }
  } # close loop
  
  # internal warning in the event that all rows in the new column are Unassigned
  # this outcome could be correct if the specified tax lvl had no assignemnts;
  # however, if assignments are expected for the specified tax lvl, this outcome
  # is incorrect and the function inputs need to be double checked
  if (isTRUE(all(grepl(pattern = "Unassigned", x = new_df1[, "new_col"],
                       ignore.case = FALSE, fixed = TRUE)))) {
    warning("in new_col all row values in output column are 'Unassigned'")
  }
  
  # format the new df
  new_df2 <- new_df1
  new_col_num <- which(names(new_df2) == 'new_col')
  if (tax_type == "Greengenes") {
    names(new_df2)[new_col_num] <- paste(classifier, "ggs", tax_lvl, sep = ".")
  }
  if (tax_type == "SILVA") {
    names(new_df2)[new_col_num] <- paste(classifier, "slv", tax_lvl, sep = ".")
  }
  return(new_df2)
  }

# inputs are incorrect, prints warning
#new_df <- trunc_tax(data = a.df, tax_type = "silva",
#                   tax_lvl = "Species", input_col = "ggTaxon",
#                   str1 = "; g__", str2 = "; s__")

# corrected inputs, no warning
#new_df <- trunc_tax(data = a.df, tax_type = "greengenes",
#                   tax_lvl = "Species", input_col = "ggTaxon",
#                   str1 = "; g__", str2 = "; s__")


### ************************************
### Taxonomic Barplots ----
### ************************************



#' Microbiome.Barplot.R
#'
#' @description Uses ggplot2 to create a stacked barplot, for example on phylum level abundances. The most abundant features (defaults to 10, based on rowMeans) will be plotted unless user specified. Anything of over 10 features will use default coloring which may be very difficult to interpret.
#' @param FEATURES Table of feature/OTU/SV counts where Samples are columns, and IDs are row names
#' @param METADATA A Table of metadata where sample names are row names.
#' @param NTOPLOT A number of features to plot.
#' @param CATEGORY A Metadata category to block samples by (faceting via ggplot2)
#' @return Barplot
#' @usage Microbiome.Barplot(table,metadata,10, "group")
#' @export




Microbiome.Barplot<-function(FEATURES,METADATA, NTOPLOT, CATEGORY){
  
  if(missing(NTOPLOT) & nrow(FEATURES)>10){NTOPLOT=10}
  else if (missing(NTOPLOT)) {NTOPLOT=nrow(FEATURES)}
  
  FEATURES<-Make.Percent(FEATURES)
  
  FEATURES<-FEATURES[order(rowMeans(FEATURES), decreasing = T),]
  if(NTOPLOT<nrow(FEATURES)){ #if left over, is added to remainder
    Other<-colSums(FEATURES[(NTOPLOT+1):nrow(FEATURES),])
    FEATURES<-rbind(FEATURES[1:NTOPLOT,], Other)
  }
  
  forplot<-TidyConvert.ToTibble(FEATURES, "Taxa") %>% gather(-Taxa, key="SampleID", value="Abundance")
  forplot$Taxa<-factor(forplot$Taxa,levels=rev(unique(forplot$Taxa)))
  if(!missing(METADATA) & !missing(CATEGORY)){
    if(TidyConvert.WhatAmI(METADATA)=="data.frame" | TidyConvert.WhatAmI(METADATA)=="matrix") {METADATA<-TidyConvert.ToTibble(METADATA, "SampleID")}
    forplot<-left_join(forplot, METADATA, by="SampleID")
  }
  
  
  library(RColorBrewer)
  # Define the number of colors you want
  nb.cols <- NTOPLOT+1
  mycolors <- colorRampPalette(brewer.pal(20, "Paired"))(nb.cols)
  
  
  PLOT<-(ggplot(forplot, aes(x=SampleID, y=Abundance, fill=Taxa))
         + geom_bar(stat="identity")
         + theme_classic()
         + ylab("% Abundance")
         + xlab("SampleID")
         + theme(axis.text.x =element_text(angle=90, hjust=1, size=7))
         + theme(legend.text = element_text(size=8,face="italic"))
  ) + scale_fill_manual(values = mycolors)
  
  
  if(NTOPLOT<=10){
    COLORS<-rev(c(
      "blue4",
      "olivedrab",
      "firebrick",
      "gold",
      "darkorchid",
      "steelblue2",
      "chartreuse1",
      "aquamarine",
      "yellow3",
      "coral",
      
      "grey"
    ))
    
    PLOT<-(PLOT + scale_fill_manual(values=COLORS, name=" "))
  }
  
  if(!missing(CATEGORY)){
    PLOT<-( PLOT + facet_grid(~get(CATEGORY), margins=FALSE, drop=TRUE, scales = "free", space = "free") )
    
  }
  return(PLOT)
}




Microbiome.Barplot2<-function(FEATURES,METADATA, NTOPLOT, CATEGORY){
  
  if(missing(NTOPLOT) & nrow(FEATURES)>10){NTOPLOT=10}
  else if (missing(NTOPLOT)) {NTOPLOT=nrow(FEATURES)}
  
  FEATURES<-Make.Percent(FEATURES)
  
  FEATURES<-FEATURES[order(rowMeans(FEATURES), decreasing = T),]
  if(NTOPLOT<nrow(FEATURES)){ #if left over, is added to remainder
    Other<-colSums(FEATURES[(NTOPLOT+1):nrow(FEATURES),])
    FEATURES<-rbind(FEATURES[1:NTOPLOT,], Other)
  }
  
  forplot<-TidyConvert.ToTibble(FEATURES, "Taxa") %>% gather(-Taxa, key="SampleID", value="Abundance")
  forplot$Taxa<-factor(forplot$Taxa,levels=rev(unique(forplot$Taxa)))
  if(!missing(METADATA) & !missing(CATEGORY)){
    if(TidyConvert.WhatAmI(METADATA)=="data.frame" | TidyConvert.WhatAmI(METADATA)=="matrix") {METADATA<-TidyConvert.ToTibble(METADATA, "SampleID")}
    forplot<-left_join(forplot, METADATA, by="SampleID")
  }
  
  
  library(RColorBrewer)
  # Define the number of colors you want
  nb.cols <- NTOPLOT+1
  mycolors <- colorRampPalette(brewer.pal(20, "PRGn"))(nb.cols)
  
  
  PLOT<-(ggplot(forplot, aes(x=SampleID, y=Abundance, fill=Taxa))
         + geom_bar(stat="identity")
         + theme_classic()
         + ylab("% Abundance")
         + xlab("SampleID")
         + theme(axis.text.x =element_text(angle=90, hjust=1, size=7))
         + theme(legend.text = element_text(size=8,face="italic"))
  ) + scale_fill_manual(values = mycolors)
  
  
  if(NTOPLOT<=10){
    COLORS<-rev(c(
      "mediumorchid4",
      "mediumpurple",
      "mediumblue"
    ))
    
    PLOT<-(PLOT + scale_fill_manual(values=COLORS, name=" "))
  }
  
  if(!missing(CATEGORY)){
    PLOT<-( PLOT + facet_grid(~get(CATEGORY), margins=FALSE, drop=TRUE, scales = "free", space = "free") )
    
  }
  return(PLOT)
}




data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}



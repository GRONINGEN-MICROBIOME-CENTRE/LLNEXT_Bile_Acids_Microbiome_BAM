############ BA Microbiome associations infants ##########
# Author: Trishla Sinha 
# Date: December 21st, 2024
# Last update: 11 january 2026

# Load libraries 
library(tidyverse)
library(reshape2)
library(dplyr)
library(pheatmap)
library(viridis)

########## Functions ################
run_invrank_dataFrame = function(data){
  invrank= function(x) {qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))} 
  data_types = unlist(lapply(data,class))
  newdata = data
  for(i in 1:ncol(data)){
    if(data_types[i]=="numeric"|data_types[i]=="integer") {
      newdata[,i] = invrank(data[,i])
    }
  }
  newdata
}
linear_model_taxa_simple_per_timepoint <- function(metadata, ID, CLR_transformed_data, pheno_list, tp = "W2" ) {
  metadata %>% filter(Timepoint == tp) -> metadata
  
  df<-merge(metadata, CLR_transformed_data, by='row.names')
  row.names(df) <- df$Row.names
  df$Row.names <- NULL
  
  Prevalent= c(colnames(CLR_transformed_data))
  #pheno_list= phenotypes
  
  Overall_result_phenos =tibble() 
  
  for (Bug in Prevalent){
    if (! Bug %in% colnames(df)){ next }
    Bug2 = paste(c("`",Bug, "`"), collapse="")
    for ( pheno in pheno_list){
      message(paste("Testing outcome:", Bug2, "| Phenotype:", pheno, "| Timepoint:", tp))
      pheno2 = paste(c("`",pheno, "`"), collapse="")
      df[is.na(df[colnames(df) == pheno]) == F, ID] -> To_keep
      df_pheno = filter(df, !!sym(ID) %in% To_keep )
      
      Model2 = as.formula(paste( c(Bug2,  " ~ clean_reads_FQ_1 + dna_conc + BATCH_NUMBER+",pheno2), collapse="" ))
      lm(Model2, df_pheno) -> resultmodel2
      
      as.data.frame(summary(resultmodel2)$coefficients)[grep(pheno, row.names(as.data.frame(summary(resultmodel2)$coefficients))),] -> Summ_simple
      Summ_simple %>% rownames_to_column("Feature") %>% as_tibble() %>% mutate(P = Summ_simple$`Pr(>|t|)`, Outcome =Bug, Timepoint=tp, Correction="techincal") -> temp_output
      rbind(Overall_result_phenos, temp_output) -> Overall_result_phenos
    }
  }
  
  p=as.data.frame(Overall_result_phenos)
  
  return(p)
}

# ------------------- Load main dataset ------------------- #
setwd("/Users/trishlasinha/Desktop/Bile_acids_NEXT")
BA<-read.delim("All_BA_2026.txt")
BA$NEXT_ID <- sub("^ID:([^-]+).*", "\\1", BA$TUBE_barcode)
BA$NEXT_ID <- paste0("LLNEXT",BA$NEXT_ID)
row.names(BA)<-BA$TUBE_barcode

# ------------------- Bile acids proportion ------------------- #
bile_acids_only <- grep("^(?!.*birthcard).*CA.*$", names(BA), ignore.case = TRUE, perl = TRUE, value = TRUE)
bile_acids_selection <- BA[, bile_acids_only]
bile_acids_selection <- na.omit(bile_acids_selection)
bile_acids_prop <- bile_acids_selection / rowSums(bile_acids_selection) * 100
bile_acids_prop <- bile_acids_prop %>%
  mutate(across(everything(), ~ replace(., is.nan(.), NA)))
bile_acids_prop<-run_invrank_dataFrame(bile_acids_prop)
# ------------------- Non-BA data ------------------- #
all_basic <- BA[, !names(BA) %in% bile_acids_only]
all_basic <- all_basic %>%
  select(-starts_with("Total"))

# ------------------- Merge ------------------- #
all_basic_bile_acids_prop <- merge(bile_acids_prop, all_basic, by="row.names")
row.names(all_basic_bile_acids_prop) <- all_basic_bile_acids_prop$Row.names
all_basic_bile_acids_prop$Row.names <- NULL
all_basic_bile_acids_prop$FK_ID=NULL
all_basic_bile_acids_prop$Family_ID=NULL
all_basic_bile_acids_prop$TUBE_barcode=NULL
all_basic_bile_acids_prop$Type=NULL

# ------------------- Subset by type ------------------- #
all_basic_bile_acids_prop_infants <- all_basic_bile_acids_prop[all_basic_bile_acids_prop$Timepoint=="B", ]
all_basic_bile_acids_prop_infants <-
  all_basic_bile_acids_prop_infants[
    , !colnames(all_basic_bile_acids_prop_infants) %in%
      c("TUDCA", "LCA", "TLCA", "GLCA", "LCA_3S",
        "GLCA_3S_quant", "TLCA_3S_qual")
  ]

# ------------------- Load and prepare taxa file ------------------- #
taxa <- read.delim("~/Desktop/LLNEXT/Analysis/taxa/LLNEXT_metaphlan_4_CLR_transformed_fil_SGB_infants_20_07_2023.txt")
metadata<-read.delim("~/Desktop/LLNEXT/Analysis/metadata/LLNEXT_metadata_15_04_2024.txt")
metadata<-metadata[, c("NG_ID", "SAMPLE_ID", "clean_reads_FQ_1", "dna_conc","BATCH_NUMBER")]
row.names(metadata)<-metadata$NG_ID
taxa_all<-merge(taxa, metadata, by="row.names")
taxa_all <- taxa_all[!duplicated(taxa_all$SAMPLE_ID), ]
row.names(taxa_all)<-taxa_all$SAMPLE_ID
taxa_all$Row.names=NULL
taxa_all$NG_ID=NULL
taxa_all$SAMPLE_ID=NULL


# ------------------- Infant processing M1 ------------------- #

all_basic_bile_acids_prop_infants$SAMPLE_ID<-paste0(all_basic_bile_acids_prop_infants$NEXT_ID, "_", "M1") 
row.names(all_basic_bile_acids_prop_infants)<-all_basic_bile_acids_prop_infants$SAMPLE_ID
phenos_test <- setdiff(colnames(all_basic_bile_acids_prop_infants),
                       c("NEXT_ID", "SAMPLE_ID" , "Timepoint" ))
all_basic_bile_acids_prop_infants$Timepoint="M1"
taxa_infant <- taxa_all[rownames(taxa_all) %in% rownames(all_basic_bile_acids_prop_infants), ]
all_basic_bile_acids_prop_infants<-all_basic_bile_acids_prop_infants[row.names(all_basic_bile_acids_prop_infants) %in% rownames(taxa_infant),]

# ------------------- infant association ------------------- #
association_infant <- linear_model_taxa_simple_per_timepoint(
  metadata = all_basic_bile_acids_prop_infants,
  ID = "SAMPLE_ID",
  CLR_transformed_data = taxa_infant ,
  pheno_list = phenos_test,
  tp = "M1"
)

association_infant$FDR<-p.adjust(association_infant$P, method = "BH")
all_infant<-merge(all_basic_bile_acids_prop_infants,taxa_infant, by="row.names" )
sig_association_infant<-association_infant[association_infant$FDR<0.05,]
write.table(association_infant, "Submission_2026/association_microbiome_BA_infant_M1.txt", row.names = F, sep = "\t")


library(ggplot2)

plot_list <- list()

for(i in 1:nrow(sig_association_infant)) {
  
  taxon <- sig_association_infant$Outcome[i]
  ba    <- sig_association_infant$Feature[i]
  
  df <- data.frame(
    BA = all_basic_bile_acids_prop_infants[[ba]],
    Taxon = taxa_infant[[taxon]]
  )
  
  p <- ggplot(df, aes(x = BA, y = Taxon)) +
    geom_point(alpha = 0.7) +
    geom_smooth(method = "lm", color = "purple", se = TRUE) +
    theme_bw() +
    labs(
      title = paste0("Infant M1: ", taxon, " vs ", ba),
      x = ba,
      y = taxon
    )
  
  plot_list[[i]] <- p
}





plot_list   





# ------------------- Infant processing W2 ------------------- #

all_basic_bile_acids_prop_infants$SAMPLE_ID<-paste0(all_basic_bile_acids_prop_infants$NEXT_ID, "_", "W2") 
row.names(all_basic_bile_acids_prop_infants)<-all_basic_bile_acids_prop_infants$SAMPLE_ID
phenos_test <- setdiff(colnames(all_basic_bile_acids_prop_infants),
                       c("NEXT_ID", "SAMPLE_ID" , "Timepoint" ))
all_basic_bile_acids_prop_infants$Timepoint="W2"
taxa_infant <- taxa_all[rownames(taxa_all) %in% rownames(all_basic_bile_acids_prop_infants), ]
all_basic_bile_acids_prop_infants<-all_basic_bile_acids_prop_infants[row.names(all_basic_bile_acids_prop_infants) %in% rownames(taxa_infant),]

# ------------------- infant association ------------------- #
association_infant_W2 <- linear_model_taxa_simple_per_timepoint(
  metadata = all_basic_bile_acids_prop_infants,
  ID = "SAMPLE_ID",
  CLR_transformed_data = taxa_infant ,
  pheno_list = phenos_test,
  tp = "W2"
)

association_infant_W2$FDR<-p.adjust(association_infant_W2$P, method = "BH")
all_infant<-merge(all_basic_bile_acids_prop_infants,taxa_infant, by="row.names" )
sig_association_infant_W2<-association_infant_W2[association_infant_W2$FDR<0.05,]

write.table(association_infant_W2, "Submission_2026/association_microbiome_BA_infant_M1.txt", row.names = F, sep = "\t")

library(ggplot2)

plot_list <- list()

for(i in 1:nrow(sig_association_infant_W2)) {
  
  taxon <- sig_association_infant_W2$Outcome[i]
  ba    <- sig_association_infant_W2$Feature[i]
  
  df <- data.frame(
    BA = all_basic_bile_acids_prop_infants[[ba]],
    Taxon = taxa_infant[[taxon]]
  )
  
  p <- ggplot(df, aes(x = BA, y = Taxon)) +
    geom_point(alpha = 0.7) +
    geom_smooth(method = "lm", color = "purple", se = TRUE) +
    theme_bw() +
    labs(
      title = paste0("Infant W2: ", taxon, " vs ", ba),
      x = ba,
      y = taxon
    )
  
  plot_list[[i]] <- p
}





plot_list   


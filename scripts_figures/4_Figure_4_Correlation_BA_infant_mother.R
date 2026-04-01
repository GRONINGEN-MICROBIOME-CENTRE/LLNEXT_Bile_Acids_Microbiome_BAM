# Figures for chapter bile acids
# Author: Trishla Sinha 
# Date: 26th of May, 2023
# Last update: 11th january 2025 

setwd("~/Desktop/Bile_acids_NEXT/")
library(vegan)
library(ggplot2)
library(reshape2)
library(tidyverse)
library(ggpubr)
library(plyr)
library(wesanderson)

####### Load functions #######
spearman <- function(x,y) {
  matchID <- intersect(rownames(x), rownames(y))
  x1 <- x[matchID,]
  y1 <- y[matchID,]
  result_cor <- matrix(nrow=ncol(x1), ncol=ncol(y1))
  rownames(result_cor) <- colnames(x1)
  colnames(result_cor) <- colnames(y1)
  result_pvalue <- matrix(nrow=ncol(x1), ncol=ncol(y1))
  for (i in 1:ncol(y1)) {
    for (j in 1:ncol(x1)) {
      cor1<-try(cor.test(x1[,j], y1[,i], method = "spearman"))
      if(class(cor1)[1] != "try-error") {
        result_cor[j,i]= cor1$estimate
        result_pvalue[j,i]= cor1$p.value
      } else {
        result_cor[j,i]= NA
        result_pvalue[j,i]= NA
      }
    }}
  result = list()
  result$p.val= result_pvalue
  result$cor= result_cor
  return(result)
}
result <- function(x,y){
  correlation <- spearman(x,y)
  a<- melt(correlation$cor)
  a<- cbind(a, melt(correlation$p.val)[,"value"])
  result= a[order(a[,4]),]
  colnames(result)=c("factor1", "factor2", "CorCoefficient","pvalue")
  return(result)
}

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

linear_model_taxa <- function(metadata, ID, CLR_transformed_data, pheno_list) {
  
  df <- metadata
  df<-merge(df, CLR_transformed_data, by='row.names')
  row.names(df) <- df$Row.names
  df$Row.names <- NULL
  
  Prevalent= c(colnames(CLR_transformed_data))
  #pheno_list= phenotypes
  
  Overall_result_phenos =tibble() 
  
  for (Bug in Prevalent){
    if (! Bug %in% colnames(df)){ next }
    Bug2 = paste(c("`",Bug, "`"), collapse="")
    for ( pheno in pheno_list){
      pheno2 = paste(c("`",pheno, "`"), collapse="")
      df[is.na(df[colnames(df) == pheno]) == F, ID] -> To_keep
      df_pheno = filter(df, !!sym(ID) %in% To_keep )
      
      Model2 = as.formula(paste( c(Bug2,  " ~ ",pheno2), collapse="" ))
      lm(Model2, df_pheno) -> resultmodel2
      
      as.data.frame(summary(resultmodel2)$coefficients)[grep(pheno, row.names(as.data.frame(summary(resultmodel2)$coefficients))),] -> Summ_simple
      Summ_simple %>% rownames_to_column("Feature") %>% as_tibble() %>% mutate(P = Summ_simple$`Pr(>|t|)`, Bug =Bug, Pheno=pheno, Model="simple") -> temp_output
      rbind(Overall_result_phenos, temp_output) -> Overall_result_phenos
    }
  }
  
  p=as.data.frame(Overall_result_phenos)
  p$FDR<-p.adjust(p$P, method = "BH")
  
  return(p)
}

############ Figure 4 Birth correlations ##########

#Figure 4
lijst_2<-read.delim("lijst_2_mother_v_baby_phenotypes_TS_NEW_NOV.txt")
row.names(lijst_2)<-lijst_2$SAMPLE_ID
bile_acids_only <- grep("^(?!.*birthcard).*CA.*$", names(lijst_2), ignore.case = TRUE, perl = TRUE, value = TRUE)
bile_acids_selection<-lijst_2[,bile_acids_only]
bile_acids_selection<-na.omit(bile_acids_selection)
bile_acids_conc<-bile_acids_selection
bile_acids_prop <- bile_acids_conc / rowSums(bile_acids_conc) * 100

lijst_2$Mother_or_Baby <- gsub("K", "Infant", lijst_2$Mother_or_Baby)
lijst_2$Mother_or_Baby <- gsub("M", "Mother", lijst_2$Mother_or_Baby)
lijst_2$ID<-paste0(lijst_2$Family_ID,"_", lijst_2$Mother_or_Baby)
lijst_2$ID <- make.unique(lijst_2$ID, sep = "_")
row.names(lijst_2)<-lijst_2$ID

numeric_cols <- names (bile_acids_prop)

# Applying a filter on bile acids 
# Calculate the number of non-zero values in each column
count_non_zero <- sapply(numeric_cols, function(x) sum(x != 0, na.rm = TRUE))

# Calculate the total number of samples (rows) in the dataframe
total_samples <- nrow(numeric_cols)

# Identify columns where the count of non-zero values is greater than 10% of the total samples
cols_to_keep <- names(count_non_zero[count_non_zero > 0.1 * total_samples])
# Removed "UDCA"  "TUDCA" "LCA"   "TLCA"  "GLCA"  "LCA3S" due to a prevalence of non-zero values <10%
# Subset the dataframe to keep only these columns
numeric_cols <- numeric_cols[, cols_to_keep]
numeric_cols$frac_12aOH=NULL # Double 
numeric_cols$frac_sulfated=NULL # Double 

infants <- numeric_cols[!grepl("_Mother$", rownames(numeric_cols)), ]
rownames(infants) <- sub("_Infant$", "", rownames(infants))

mothers <- numeric_cols[grepl("_Mother$", rownames(numeric_cols)), ]
rownames(mothers) <- sub("_Mother$", "", rownames(mothers))

infants <- infants[row.names(infants) %in% row.names(mothers), ]
mothers <- mothers[row.names(mothers) %in% row.names(infants), ]
infants <- infants[order(row.names(infants)), ]
mothers <- mothers[order(row.names(mothers)), ]

# Now you can run spearman correlations on these two datasets 
corrlation_spearman <- result(infants, mothers)
names (corrlation_spearman)[1]<-"Infant"
names (corrlation_spearman)[2]<-"Mother"

pheno_spearman_FDR <- corrlation_spearman
pheno_spearman_FDR$FDR<-p.adjust(pheno_spearman_FDR$pvalue, method = "BH")
setwd("/Users/trishlasinha/Desktop/Bile_acids_NEXT/Submission_2026")
write.table(pheno_spearman_FDR, "spearman_correlation_infant_mother.txt", sep = "\t")

pheno_spearman_FDR$stars <- cut(pheno_spearman_FDR$FDR, breaks=c(-Inf, 0.0001, 0.001, 0.01, Inf), label=c("***", "**", "*", ""))  # Create column of significance labels

pdf('chapter_hilde/Figures/Figure_4.pdf', width=10, height=8)
ggplot(pheno_spearman_FDR, aes(Infant, Mother, fill = CorCoefficient))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Spearman\nCorrelation") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 10, hjust = 1))+
  theme(axis.text.y = element_text(size = 10))+
  geom_text(aes(label=stars), color="black", size=5) + 
  coord_fixed()
dev.off()

#write.table(pheno_spearman_FDR, "spearman_correlation_sig_mother_markers_infant_markers.txt", quote = F, row.names = F, sep = "\t")

# Run linear regressions on both mom and baby separately 

# Filter the dataframe to rows ending in "_Mother" to first run with mother birth
phenos <- lijst_2[grep("_Mother$", rownames(lijst_2)), ]
phenos<-phenos[,c(43:46, 52:54, 57:61 )]
rownames(phenos) <- sub("_Mother$", "", rownames(phenos))
BA<-mothers[,c(3,22,23,24,25,27)]

overlap<-intersect(row.names(phenos), row.names(BA))
BA<-BA[overlap,]
phenos$ID<-row.names(phenos)
phenos[sapply(phenos, is.character)] <- lapply(phenos[sapply(phenos, is.character)], 
                                       as.factor)
# phenos should have a column with ID's 

crossectional_mothers <- linear_model_taxa(phenos, "ID", BA,c("mother_lifelinesbirthcards_prepreg_bmi_kg_m2", "mother_birthcard_age_at_delivery"))


#write.table(crossectional_mothers, "mother_B_association_phenotypes_BA.txt", sep = "\t", row.names = F, quote = F)


# Filter the dataframe to rows ending in "_Mother" to first run with mother birth
phenos <- lijst_2[grep("_Infant$", rownames(lijst_2)), ]
phenos<-phenos[,c(43:46, 52:54, 57:61 )]
rownames(phenos) <- sub("_Infant$", "", rownames(phenos))
BA<-infants[,c(3,22,23,24,25,27)]

overlap<-intersect(row.names(phenos), row.names(BA))
BA<-BA[overlap,]
phenos$ID<-row.names(phenos)
phenos[sapply(phenos, is.character)] <- lapply(phenos[sapply(phenos, is.character)], 
                                               as.factor)
# phenos should have a column with ID's 


crossectional_infants <- linear_model_taxa(phenos, "ID", BA,c("mother_lifelinesbirthcards_prepreg_bmi_kg_m2","mother_birthcard_gestational_age_weeks", "birth_deliverybirthcard_mode_binary",
                                             "mother_birthcard_age_at_delivery", "infant_misc_sex","infant_growth_birth_weight_kg", "infant_growth_rate_weight"  , 
                                            "infant_birthcard_feeding_mode_after_birth"  ,"infant_birthcard_apgar_score_1min", "infant_birthcard_apgar_score_5min"))

#write.table(crossectional_infants, "infant_B_association_phenotypes_BA.txt", sep = "\t", row.names = F, quote = F)

all_infants<-merge(infants, phenos, by="row.names")
ggplot(all_infants, aes(birth_deliverybirthcard_mode_binary, y=Ratio_sulfated_nonsulfated, fill=birth_deliverybirthcard_mode_binary, color=birth_deliverybirthcard_mode_binary)) +
  scale_fill_manual(values = wes_palette("Royal1", n = 2))+
  scale_color_manual(values = wes_palette("Royal1", n = 2))+
  #scale_fill_manual(values=c("#e60049", "#0bb4ff", "#50e991", "#f46a9b", "#9b19f5", "#ffa300", "#dc0ab4", "#b3d4ff", "#00bfa0", "#b30000", "#7c1158", "#4421af"))+
  #scale_color_manual(values=c("#e60049", "#0bb4ff", "#50e991", "#f46a9b", "#9b19f5", "#ffa300", "#dc0ab4", "#b3d4ff", "#00bfa0", "#b30000", "#7c1158", "#4421af"))+
  geom_boxplot(alpha=0.4, outlier.colour = NA)+
  #geom_point(alpha=0.6,
             #position = position_jitterdodge(jitter.width = 0.3, jitter.height = 0, outlier.colo))+
  # geom_jitter() +
  ggtitle("")+
  #geom_path(aes(group =NEXT_ID, col=infant_place_delivery), size = .8) +
  theme_bw()+labs(x="", y = "Ratio_sulfated_nonsulfated")+
  theme(
    plot.title = element_text(color="black", size=18, face="bold"),
    axis.title.x = element_text(color="black", size=18, face="bold"),
    axis.title.y = element_text(color="black", size=18, face="bold"),
    axis.text.y = element_text(face="bold", size=10),
    axis.text.x = element_text(size=22, angle = 60, hjust = 1),
    legend.position = "none")
#strip.text.x = element_text(size = 10))
library(wesanderson)



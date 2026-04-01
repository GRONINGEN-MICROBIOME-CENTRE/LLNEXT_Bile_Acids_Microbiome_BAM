############ BA Microbiome associations mothers ##########
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

# ------------------- Subset by timepoint ------------------- #
all_basic_bile_acids_prop_P12 <- all_basic_bile_acids_prop[all_basic_bile_acids_prop$Timepoint=="P12", ]
all_basic_bile_acids_prop_P40 <- all_basic_bile_acids_prop[all_basic_bile_acids_prop$Timepoint=="P40", ]
all_basic_bile_acids_prop_infants <- all_basic_bile_acids_prop[all_basic_bile_acids_prop$Timepoint=="B", ]

# ------------------- Load and prepare taxa file ------------------- #
taxa <- read.delim("~/Desktop/LLNEXT/Analysis/taxa/NEXT_metaphlan_4_CLR_transformed_fil_30_percent_SGB_mothers_03_08_2023.txt")
metadata<-read.delim("~/Desktop/LLNEXT/Analysis/metadata/LLNEXT_metadata_15_04_2024.txt")
metadata<-metadata[, c("NG_ID", "SAMPLE_ID", "clean_reads_FQ_1", "dna_conc","BATCH_NUMBER")]
row.names(metadata)<-metadata$NG_ID
taxa_all<-merge(taxa, metadata, by="row.names")
row.names(taxa_all)<-taxa_all$SAMPLE_ID
taxa_all$Row.names=NULL
taxa_all$NG_ID=NULL
taxa_all$SAMPLE_ID=NULL

# ------------------- P12 processing ------------------- #

all_basic_bile_acids_prop_P12$SAMPLE_ID<-paste0(all_basic_bile_acids_prop_P12$NEXT_ID, "_", "P12") 
row.names(all_basic_bile_acids_prop_P12)<-all_basic_bile_acids_prop_P12$SAMPLE_ID
phenos_test <- setdiff(colnames(all_basic_bile_acids_prop_P12),
                       c("NEXT_ID", "SAMPLE_ID" , "Timepoint" ))
all_basic_bile_acids_prop_P12$Timepoint="P12"
taxa_P12 <- taxa_all[rownames(taxa_all) %in% rownames(all_basic_bile_acids_prop_P12), ]
all_basic_bile_acids_prop_P12<-all_basic_bile_acids_prop_P12[row.names(all_basic_bile_acids_prop_P12) %in% rownames(taxa_P12),]

# ------------------- P12 association ------------------- #
association_P12 <- linear_model_taxa_simple_per_timepoint(
  metadata = all_basic_bile_acids_prop_P12,
  ID = "SAMPLE_ID",
  CLR_transformed_data = taxa_P12,
  pheno_list = phenos_test,
  tp = "P12"
)

association_P12$FDR<-p.adjust(association_P12$P, method = "BH")
all_P12<-merge(all_basic_bile_acids_prop_P12,taxa_P12, by="row.names" )
sig_association_P12<-association_P12[association_P12$FDR<0.05,]

library(ggplot2)

plot_list <- list()

for(i in 1:nrow(sig_association_P12)) {
  
  taxon <- sig_association_P12$Outcome[i]
  ba    <- sig_association_P12$Feature[i]
  
  df <- data.frame(
    BA = all_basic_bile_acids_prop_P12[[ba]],
    Taxon = taxa_P12[[taxon]]
  )
  
  p <- ggplot(df, aes(x = BA, y = Taxon)) +
    geom_point(alpha = 0.7) +
    geom_smooth(method = "lm", color = "purple", se = TRUE) +
    theme_bw() +
    labs(
      title = paste0("P12: ", taxon, " vs ", ba),
      x = ba,
      y = taxon
    )
  
  plot_list[[i]] <- p
}

plot_list   

wide_data <- dcast(sig_association_P12, Outcome ~ Feature, value.var = "t value" )  
row.names(wide_data)<-wide_data$Outcome
wide_data$Outcome=NULL
wide_data[is.na(wide_data)] <- 0

heatmap_data <- as.matrix(t(wide_data))

# Plot clustered heatmap
pdf('Figures/Figure_2_a.pdf', width=10, height=8)
BA_SGB_P12 <- pheatmap(
  heatmap_data,
  cluster_rows = TRUE,
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  color = viridis(100, option = "D"),
  display_numbers = FALSE,
  main = "",
  angle_col = 90,
  fontsize_col = 10
)

BA_SGB_P12
dev.off()


write.table(association_P12, "Submission_2026/association_microbiome_BA_P12.txt", row.names = F, sep = "\t")

# ------------------- P40 processing ------------------- #

all_basic_bile_acids_prop_P40$SAMPLE_ID <-
  paste0(all_basic_bile_acids_prop_P40$NEXT_ID, "_", "B")

row.names(all_basic_bile_acids_prop_P40) <-
  all_basic_bile_acids_prop_P40$SAMPLE_ID

phenos_test <- setdiff(
  colnames(all_basic_bile_acids_prop_P40),
  c("NEXT_ID", "SAMPLE_ID", "Timepoint")
)

all_basic_bile_acids_prop_P40$Timepoint <- "P40"

taxa_P40 <- taxa_all[
  rownames(taxa_all) %in% rownames(all_basic_bile_acids_prop_P40),
]

all_basic_bile_acids_prop_P40 <-
  all_basic_bile_acids_prop_P40[
    row.names(all_basic_bile_acids_prop_P40) %in% rownames(taxa_P40),
  ]

# ------------------- P40 association ------------------- #

association_P40 <- linear_model_taxa_simple_per_timepoint(
  metadata = all_basic_bile_acids_prop_P40,
  ID = "SAMPLE_ID",
  CLR_transformed_data = taxa_P40,
  pheno_list = phenos_test,
  tp = "P40"
)

association_P40$FDR <- p.adjust(association_P40$P, method = "BH")

all_P40 <- merge(
  all_basic_bile_acids_prop_P40,
  taxa_P40,
  by = "row.names"
)

sig_association_P40 <- association_P40[association_P40$FDR < 0.05, ]

write.table(association_P40, "Submission_2026/association_microbiome_BA_P40.txt", row.names = F, sep = "\t")



########### Bile Acid Dehydrogenase Genes #####################

setwd("/Users/trishlasinha/Desktop/Bile_acids_NEXT")

BADHGENES<-read.delim("/Users/trishlasinha/Desktop/Bile_acids_NEXT/NEXT_Bile_Acid_shortbred.csv")
BADHGENES$NG_ID<-substr(BADHGENES$Sample, 0, 13) 
BADHGENES$NG_ID<-gsub("_", "", BADHGENES$NG_ID)
BADHGENES<-BADHGENES[!duplicated(BADHGENES$NG_ID),]
row.names(BADHGENES)<-BADHGENES$NG_ID
BADHGENES$Sample=NULL
BADHGENES$NG_ID=NULL

metadata<-read.delim("LLNEXT_metadata_15_04_2024.txt")
row.names(metadata)<-metadata$NG_ID
all<-merge(BADHGENES, metadata, by="row.names")
row.names(all)<-all$Row.names
all$Row.names=NULL

lijst_1<-read.delim("lijst_1_MB_v_P12_phenotypes_TS_NEW_NOV.txt")
row.names(lijst_1)<-lijst_1$SAMPLE_ID
lijst_1$SAMPLE_ID<-row.names(lijst_1)
lijst_1$Timepoint[lijst_1$Timepoint == "B"] <- "P40"
lijst_1Timepoint <- factor(lijst_1$Timepoint, levels = c("P12", "P40"))
lijst_1 <- lijst_1[!is.na(lijst_1$CA), ]

filtered_all <- all[all$SAMPLE_ID %in% lijst_1$SAMPLE_ID, ]
filtered_all$sum_choloylglycine_hydrolase <- rowSums(filtered_all[, grepl("^tr\\.", names(filtered_all))])
filtered_all$Timepoint <-filtered_all$Timepoint_categorical
filtered_all$Timepoint[filtered_all$Timepoint == "B"] <- "P40"
filtered_all$Timepoint <- factor(filtered_all$Timepoint, levels = c("P12", "P40"))


# Violin plot for secondary bile acids using Timepoint
ggplot(lijst_1, 
       aes(x = Timepoint, y = Fraction_secondary, fill = Timepoint, color = Timepoint)) + 
  geom_violin(aes(fill = Timepoint, alpha = 0.2, width = 0.3, trim = FALSE)) + 
  geom_boxplot(alpha = 0.9, outlier.colour = NA, width = 0.03, fill = "white") + 
  geom_point(alpha = 0.6, 
             position = position_jitterdodge(jitter.width = 0.3, jitter.height = 0), size = 0.2) + 
  labs(x = "", y = "Fraction Secondary Bile Acids") + 
  scale_fill_manual(values = c("#FF0000", "#00007F")) + 
  scale_color_manual(values = c("#FF0000", "#00007F")) + 
  theme_classic() + 
  theme(plot.title = element_text(size = 18, hjust = 0.5, face = "bold"), 
        axis.text = element_text(size = 14), 
        axis.title = element_text(size = 16), 
        axis.title.y = element_text(margin = margin(r = 10)), 
        axis.title.x = element_text(margin = margin(t = 10))) +
  guides(fill = FALSE, color = FALSE, alpha = FALSE)


# Violin plot for  choloylglycine_hydrolase using Timepoint
ggplot(filtered_all, 
       aes(x = Timepoint, y = sum_choloylglycine_hydrolase, fill = Timepoint, color = Timepoint)) + 
  geom_violin(aes(fill = Timepoint, alpha = 0.2, width = 0.3, trim = FALSE)) + 
  geom_boxplot(alpha = 0.9, outlier.colour = NA, width = 0.03, fill = "white") + 
  geom_point(alpha = 0.6, 
             position = position_jitterdodge(jitter.width = 0.3, jitter.height = 0), size = 0.2) + 
  labs(x = "", y = "Total Choloylglycine Hydrolase") + 
  scale_fill_manual(values = c("#FF0000", "#00007F")) + 
  scale_color_manual(values = c("#FF0000", "#00007F")) + 
  theme_classic() + 
  theme(plot.title = element_text(size = 18, hjust = 0.5, face = "bold"), 
        axis.text = element_text(size = 14), 
        axis.title = element_text(size = 16), 
        axis.title.y = element_text(margin = margin(r = 10)), 
        axis.title.x = element_text(margin = margin(t = 10))) +
  guides(fill = FALSE, color = FALSE, alpha = FALSE)
names (filtered_all)


# Load necessary libraries
library(dplyr)
library(tidyr)
library(stats)


bile_acid_columns <- grep("^tr\\.", colnames(filtered_all), value = TRUE)

# Violin plot for  tr.A0A414INL9.A0A414INL9_BACUN using Timepoint
ggplot(filtered_all, 
       aes(x = Timepoint, y = tr.A0A414INL9.A0A414INL9_BACUN, fill = Timepoint, color = Timepoint)) + 
  geom_violin(aes(fill = Timepoint, alpha = 0.2, width = 0.3, trim = FALSE)) + 
  geom_boxplot(alpha = 0.9, outlier.colour = NA, width = 0.03, fill = "white") + 
  geom_point(alpha = 0.6, 
             position = position_jitterdodge(jitter.width = 0.3, jitter.height = 0), size = 0.2) + 
  labs(x = "", y = "tr.A0A414INL9.A0A414INL9_BACUN") + 
  scale_fill_manual(values = c("#FF0000", "#00007F")) + 
  scale_color_manual(values = c("#FF0000", "#00007F")) + 
  theme_classic() + 
  theme(plot.title = element_text(size = 18, hjust = 0.5, face = "bold"), 
        axis.text = element_text(size = 14), 
        axis.title = element_text(size = 16), 
        axis.title.y = element_text(margin = margin(r = 10)), 
        axis.title.x = element_text(margin = margin(t = 10))) +
  guides(fill = FALSE, color = FALSE, alpha = FALSE)
names (filtered_all)



names (lijst_1)


ggplot(lijst_1, 
       aes(x = Timepoint, y = Fraction_conjugated, fill = Timepoint, color = Timepoint)) + 
  geom_violin(aes(fill = Timepoint, alpha = 0.2, width = 0.3, trim = FALSE)) + 
  geom_boxplot(alpha = 0.9, outlier.colour = NA, width = 0.03, fill = "white") + 
  geom_point(alpha = 0.6, 
             position = position_jitterdodge(jitter.width = 0.3, jitter.height = 0), size = 0.2) + 
  labs(x = "", y = "Total choloylglycine hydrolase") + 
  scale_fill_manual(values = c("#FF0000", "#00007F")) + 
  scale_color_manual(values = c("#FF0000", "#00007F")) + 
  theme_classic() + 
  theme(plot.title = element_text(size = 18, hjust = 0.5, face = "bold"), 
        axis.text = element_text(size = 14), 
        axis.title = element_text(size = 16), 
        axis.title.y = element_text(margin = margin(r = 10)), 
        axis.title.x = element_text(margin = margin(t = 10))) +
  guides(fill = FALSE, color = FALSE, alpha = FALSE)

#conjugated_to_sum <- colnames(lijst_1)[grepl("T", colnames(lijst_1)) | grepl("G", colnames(lijst_1))]
#lijst_1$check_congugated<-lijst_1$TCA + lijst_1$GCA + lijst_1$TCDCA + lijst_1$GCDCA + lijst_1$TDCA + lijst_1$GDCA + lijst_1$TLCA + lijst_1$GLCA +lijst_1$TUDCA +lijst_1$GUDCA+lijst_1$TLCA3S +lijst_1$GLCA3S
#lijst_1$check_uncongugated<-lijst_1$CA + lijst_1$CDCA + lijst_1$DCA + lijst_1$LCA +lijst_1$UDCA
#lijst_1$total_sum <-lijst_1$TCA + lijst_1$GCA + lijst_1$TCDCA + lijst_1$GCDCA + lijst_1$TDCA + lijst_1$GDCA + lijst_1$TLCA + lijst_1$GLCA +lijst_1$TUDCA +lijst_1$GUDCA+lijst_1$CA + lijst_1$CDCA + lijst_1$DCA + lijst_1$LCA +lijst_1$UDCA +lijst_1$LCA3S +lijst_1$GLCA3S+ lijst_1$TLCA3S
#lijst_1$check_fraction_conjugated<-lijst_1$check_congugated/lijst_1$total_sum

#lijst_1$ratio_unconjugated_conjugated_concentrations = (lijst_1$CA + lijst_1$CDCA + lijst_1$DCA + lijst_1$LCA+lijst_1$UDCA)/(lijst_1$TCA + lijst_1$GCA + lijst_1$TCDCA + lijst_1$GCDCA + lijst_1$TDCA + lijst_1$GDCA + lijst_1$TLCA + lijst_1$GLCA + lijst_1$TUDCA+lijst_1$GUDCA)

row.names(filtered_all)<-filtered_all$SAMPLE_ID
all<-merge(filtered_all, lijst_1, by="row.names")

ggplot(all, aes(x =Fraction_conjugated , y = sum_choloylglycine_hydrolase)) +
  geom_point(color = "blue", size = 2) +  # Scatter plot points
  geom_smooth(method = "lm", color = "red", se = FALSE) +  # Add a linear regression line
  theme_minimal()


# Extract subsets of lijst_1 based on the Timepoint column
data_P12_full <- lijst_1[lijst_1$Timepoint == "P12", ]
data_P40_full <- lijst_1[lijst_1$Timepoint == "P40", ]

# Order the subsets by NEXT_ID
data_P12_ordered <- data_P12_full[order(data_P12_full$NEXT_ID), ]
data_P40_ordered <- data_P40_full[order(data_P40_full$NEXT_ID), ]

# Identify numeric columns for testing (excluding NEXT_ID and Timepoint)
numeric_columns <- names(data_P12_ordered)[sapply(data_P12_ordered, is.numeric)]

# Initialize a list to store the test results
test_results <- list()

# Initialize a data frame to store the test results
test_results <- data.frame(
  Column = character(),
  W_Statistic = numeric(),
  P_Value = numeric(),
  stringsAsFactors = FALSE
)

# Loop over each numeric column and perform the paired Wilcoxon test
for (col in numeric_columns) {
  # Perform paired Wilcoxon test for each column
  test <- wilcox.test(data_P12_ordered[[col]], data_P40_ordered[[col]], paired = TRUE)
  
  # Append the results to the test_results data frame
  test_results <- rbind(test_results, data.frame(
    Column = col,
    W_Statistic = test$statistic,
    P_Value = test$p.value
  ))
}

test_results$FDR <- p.adjust(test_results$P_Value, method = "BH")
write.table(test_results, "wilcoxen_paired.txt", row.names = F, sep = "\t")

ggplot(lijst_1, 
      aes(x = Timepoint, y = Total_BAs, fill = Timepoint, color = Timepoint)) + 
  geom_violin(aes(fill = Timepoint, alpha = 0.2, width = 0.3, trim = FALSE)) + 
  geom_boxplot(alpha = 0.9, outlier.colour = NA, width = 0.03, fill = "white") + 
  geom_point(alpha = 0.6, 
             position = position_jitterdodge(jitter.width = 0.3, jitter.height = 0), size = 0.2) + 
  labs(x = "", y = "Total choloylglycine hydrolase") + 
  scale_fill_manual(values = c("#FF0000", "#00007F")) + 
  scale_color_manual(values = c("#FF0000", "#00007F")) + 
  theme_classic() + 
  theme(plot.title = element_text(size = 18, hjust = 0.5, face = "bold"), 
        axis.text = element_text(size = 14), 
        axis.title = element_text(size = 16), 
        axis.title.y = element_text(margin = margin(r = 10)), 
        axis.title.x = element_text(margin = margin(t = 10))) +
  guides(fill = FALSE, color = FALSE, alpha = FALSE)


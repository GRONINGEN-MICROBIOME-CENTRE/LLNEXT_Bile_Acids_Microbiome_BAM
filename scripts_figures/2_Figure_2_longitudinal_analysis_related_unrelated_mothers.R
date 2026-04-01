### Longitudinal analysis Maternal Bile Acids ###

library(plyr)

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

linear_model_taxa_simple <- function(metadata, ID, CLR_transformed_data, pheno_list) {
  
  df <- metadata
  row.names(df) <- df[,ID]
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


# Analysis 
lijst_1<-read.delim("lijst_1_MB_v_P12_phenotypes_TS_NEW_NOV.txt")
row.names(lijst_1)<-lijst_1$SAMPLE_ID
bile_acids_only <- grep("^(?!.*birthcard).*CA.*$", names(lijst_1), ignore.case = TRUE, perl = TRUE, value = TRUE)
bile_acids_selection<-lijst_1[,bile_acids_only]
bile_acids_selection<-na.omit(bile_acids_selection)
bile_acids_conc<-bile_acids_selection
bile_acids_prop <- bile_acids_conc / rowSums(bile_acids_conc) * 100

all_basic<- lijst_1[, !names(lijst_1) %in% bile_acids_only]

all_basic_bile_acids_prop <- merge(bile_acids_prop, all_basic, by="row.names")
row.names(all_basic_bile_acids_prop)<-all_basic_bile_acids_prop$Row.names
all_basic_bile_acids_prop$Row.names=NULL

all_basic_bile_acids_prop<- all_basic_bile_acids_prop[match(rownames(bile_acids_prop), rownames(all_basic_bile_acids_prop)), ]

#my_pseudocount_normal=min(bile_acids_prop [bile_acids_prop !=0])/2
#all_ba_p_dis<-vegdist(bile_acids_prop, method = "aitchison", na.rm = T, pseudocount=my_pseudocount_normal)
all_ba_p_dis<-vegdist(bile_acids_prop, method = 'canberra', na.rm = T)
dist_matrix <-as.matrix(all_ba_p_dis)
dist_matrix [upper.tri(dist_matrix )] <- NA 
distances_long<-melt(dist_matrix)
head (distances_long)

# Remove duplicate distance 
distances_long=distances_long %>%
  drop_na()

# Add related or unrelated on basis of the FAMID
distances_long$relationship <- ifelse(substr(distances_long$Var1, 1, 12) == substr(distances_long$Var2, 1, 12), "related", "unrelated")


# Removing same sample pairs 
distances_long <- distances_long[distances_long$Var1 != distances_long$Var2, ]


ggplot(distances_long, aes(x = relationship, y = value, color=relationship)) +
  geom_boxplot()+
  geom_violin(alpha = 0.5) +
  ylab("Eucledian distances") +
  geom_jitter(position = position_jitter(seed = 1, width = 0.2), alpha = 0.2) +
  xlab("Relationship") +
  theme_bw()+
  ggtitle("Sharing of Bile Acid Profile between related mother-mother pairs")


ggplot(distances_long, aes(x=value,fill=relationship)) + 
  geom_density(alpha=0.4)+
  theme_classic()+labs(x="Eucledian distances", y = "Density")+
  ggtitle("Sharing of Bile Acid Profile between mother-mother pairs")
theme(legend.title=element_blank())

#Sharing_related_unrelated_mother_infant_pairs

p<-ggplot(distances_long, aes(x=value, fill=relationship, color=relationship)) +
  geom_density(alpha=0.4)+scale_color_manual(values=c("#56B4E9", "#E69F00" ))+
  scale_fill_manual(values=c("#56B4E9", "#E69F00" ))+theme_bw()
p

mu <- ddply(distances_long, "relationship", summarise, grp.median=median(value))
# Add mean lines
p+geom_vline(data=mu, aes(xintercept=grp.median, color=relationship),
             linetype="dashed")+
  theme_bw()+labs(x="Canberra distances", y = "Density")+
  ggtitle("Sharing of Bile Acid Profile between mother-infant pairs")+
  theme(legend.title=element_blank())


# Preallocate a vector for efficiency
p_values <- numeric(1000)

for (i in 1:1000) {
  p_values[i] <- wilcox.test(distances_long$value ~ sample(distances_long$relationship))$p.value
}

# Convert the vector to a dataframe
results_perm_wilcoxen <- data.frame(p_values = p_values)

# Plot the histogram
ggplot(results_perm_wilcoxen, aes(x = p_values)) + 
  geom_histogram(color = "black", fill = "white")

# Actual Wilcoxon test
actual_test <- wilcox.test(distances_long$value ~ distances_long$relationship)

# Compare the actual p-value to the distribution
actual_p_value <- actual_test$p.value
# actual_p_value=0.004717323

# Plot the histogram with a line indicating the actual p-value
ggplot(results_perm_wilcoxen, aes(x = p_values)) + 
  geom_histogram(color = "black", fill = "white", bins = 30) +  # Adjust number of bins as needed
  geom_vline(xintercept = actual_p_value, color = "red", linetype = "dashed") +
  ggtitle("Distribution of Permutation P-values with Actual P-value") +
  xlab("P-values") +
  ylab("Frequency") +
  theme_minimal()

# Calculate the proportion of permuted p-values less than or equal to the actual p-value
proportion_less_than_actual <- mean(results_perm_wilcoxen$p_values <= actual_p_value)

# Print the proportion
print(proportion_less_than_actual)


############## Now look at the correlations between specific bile acids #####################

numeric_cols <- all_basic_bile_acids_prop[, sapply(all_basic_bile_acids_prop, is.numeric)]
numeric_cols <- numeric_cols[, !grepl("^mother_|^infant_", names(numeric_cols))]
P12 <- numeric_cols[!grepl("_B$", rownames(numeric_cols)), ]
rownames(P12) <- sub("_P12", "", rownames(P12))

MB <- numeric_cols[grepl("_B$", rownames(numeric_cols)), ]
rownames(MB) <- sub("_B$", "", rownames(MB))

P12 <- P12[row.names(P12) %in% row.names(MB), ]
MB<- MB[row.names(MB) %in% row.names(P12), ]

P12 <- P12[order(row.names(P12)), ]
MB <- MB[order(row.names(MB)), ]

# Now you can run spearman correlations on these two datasets 
corrlation_spearman <- result(P12, MB)
names (corrlation_spearman)[1]<-"P12"
names (corrlation_spearman)[2]<-"MB"

pheno_spearman_FDR <- corrlation_spearman
pheno_spearman_FDR$FDR<-p.adjust(pheno_spearman_FDR$pvalue, method = "BH")
#pheno_spearman_FDR <- pheno_spearman_FDR %>% 
 # filter(FDR<0.01)
pheno_spearman_FDR$stars <- cut(pheno_spearman_FDR$FDR, breaks=c(-Inf, 0.0001, 0.001, 0.01, Inf), label=c("***", "**", "*", ""))  # Create column of significance labels

ggplot(pheno_spearman_FDR, aes(P12, MB, fill = CorCoefficient))+
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

merged<-merge(P12, MB, by="row.names")

q <- ggplot(merged, aes(x = TLCA3S.x, y = TLCA3S.y)) +
  geom_point() + 
  geom_smooth(method = "lm", se = T, aes()) + # Added regression line
  labs(title = "",
       x = "Mother P12", 
       y = "Mother B") +
  theme_minimal() 

q

############## Checking the effect of timepoint on the bile acid composition ################




# Figures for chapter bile acids
# Author: Trishla Sinha 
# Date: 25th December, 2023

setwd("~/Desktop/Bile_acids_NEXT/")
library(vegan)
library(ggplot2)
library(reshape2)
library(tidyverse)
library(ggpubr)
library(dplyr)

############ Mother and Infant birth ##########

#Figure 3
lijst_2<-read.delim("lijst_2_mother_v_baby_phenotypes_TS_NEW_NOV.txt")
row.names(lijst_2)<-lijst_2$SAMPLE_ID
bile_acids_only <- grep("^(?!.*birthcard).*CA.*$", names(lijst_2), ignore.case = TRUE, perl = TRUE, value = TRUE)
bile_acids_selection<-lijst_2[,bile_acids_only]
bile_acids_selection<-na.omit(bile_acids_selection)
bile_acids_conc<-bile_acids_selection
bile_acids_prop <- bile_acids_conc / rowSums(bile_acids_conc) * 100

#all_basic<-lijst_2[,c(1,3,6,31:41,43:61)]
all_basic<-lijst_2

all_basic_bile_acids_prop <- all_basic[match(rownames(bile_acids_prop), rownames(all_basic)),]

all_ba_p_dis<-vegdist(bile_acids_prop, method = 'canberra', na.rm = T)
all_ba_p_dis_mds<-cmdscale(all_ba_p_dis, k=5, eig = T)
all_ba_p_pcoa <- data.frame(all_ba_p_dis_mds$points)

pdf('chapter_hilde/Figures/Figure_3_c.pdf', width=6, height=5)
p_ba_p_pcoa_timepoint<-ggplot(all_ba_p_pcoa ,aes(X1,X2, color = all_basic_bile_acids_prop$Mother_or_Baby))+
  geom_point(size = 2,alpha = 0.5)+
  stat_ellipse(aes(group = all_basic_bile_acids_prop$Mother_or_Baby, fill = all_basic_bile_acids_prop$Mother_or_Baby, color = all_basic_bile_acids_prop$Mother_or_Baby) ,type = "norm",linetype = 2,geom = "polygon",alpha = 0.05,show.legend = F)+
  xlab(paste("PCo1=",round(all_ba_p_dis_mds$eig[1],digits = 2),"%",sep = ""))+
  ylab(paste("PCo2=",round(all_ba_p_dis_mds$eig[2],digits = 2),"%",sep = ""))+
  scale_color_manual(name=NULL, 
                     breaks = c("M",  "K"),
                     labels = c("Mother P40         ", "Cord blood"),
                     values = c("#00007F", "#008000"))+
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.line.x =  element_line(),
        axis.line.y = element_line(),
        legend.position = 'bottom',
        legend.title = element_blank(),
        legend.key = element_rect(fill = NA), 
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA),
        panel.background = element_rect(fill = NA))
p_ba_p_pcoa_timepoint<-ggExtra::ggMarginal(p_ba_p_pcoa_timepoint, type = "histogram", groupColour = F, groupFill = TRUE,
                                           xparams = list(bins = 60, alpha = 0.5,position = 'identity', color = 'white'),
                                           yparams = list(bins = 60, alpha = 0.5,position = 'identity', color = 'white'))
p_ba_p_pcoa_timepoint
dev.off()
adonis2(bile_acids_prop~all_basic_bile_acids_prop$Mother_or_Baby)
wilcox.test(all_ba_p_pcoa$X1 ~all_basic_bile_acids_prop$Mother_or_Baby)
#p-value < 2.2e-16
wilcox.test(all_ba_p_pcoa$X2 ~all_basic_bile_acids_prop$Mother_or_Baby)
#p-value = 0.05967

adonis2(bile_acids_prop~all_basic_bile_acids_prop$Mother_or_Baby)
#adonis2(formula = bile_acids_prop ~ all_basic_bile_acids_prop$Mother_or_Baby)
#Df SumOfSqs      R2      F Pr(>F)    
#all_basic_bile_acids_prop$Mother_or_Baby   1   17.165 0.36947 280.09  0.001 ***

# Replace "K" with "Cord Blood" in the Timepoint column
all_basic$Mother_or_Baby[all_basic$Mother_or_Baby == "K"] <- "Cord Blood"
all_basic$Mother_or_Baby[all_basic$Mother_or_Baby == "M"] <- "Mother P40"
all_basic$Mother_or_Baby <- factor(all_basic$Mother_or_Baby, levels = c("Mother P40", "Cord Blood"))

#figure_3a

pdf('chapter_hilde/Figures/Figure_3_a.pdf', width=6, height=5)
figure_3_a<-ggplot(all_basic, 
                   aes(x=Mother_or_Baby, y=C4_nM, fill = Mother_or_Baby, color=Mother_or_Baby)) +
  geom_violin(aes(fill=Mother_or_Baby, alpha=0.2, width=0.3, trim = F)) +
  geom_boxplot(alpha=0.9, outlier.colour = NA, width=0.03, fill="white") +
  geom_point(alpha=0.6,
             position = position_jitterdodge(jitter.width = 0.3, jitter.height = 0), size=0.2) +
  labs(x = "", y = "C4 (nM)") +
  scale_fill_manual(values = c("#00007F", "#008000")) +
  scale_color_manual(values = c("#00007F", "#008000")) +
  theme_classic() +
  theme(plot.title = element_text(size=18, hjust = 0.5, face="bold"), 
        axis.text = element_text(size=14),
        axis.title = element_text(size=16), 
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10))) +
  guides(fill = FALSE, color = FALSE, alpha=FALSE) 
figure_3_a
dev.off()
#figure_3b

pdf('chapter_hilde/Figures/Figure_3_b.pdf', width=6, height=5)
figure_3_b<-ggplot(all_basic, 
                   aes(x=Mother_or_Baby, y=Total_BAs, fill = Mother_or_Baby, color=Mother_or_Baby)) +
  geom_violin(aes(fill=Mother_or_Baby, alpha=0.2, width=0.3, trim = F)) +
  geom_boxplot(alpha=0.9, outlier.colour = NA, width=0.03, fill="white") +
  geom_point(alpha=0.6,
             position = position_jitterdodge(jitter.width = 0.3, jitter.height = 0), size=0.2) +
  labs(x = "", y = "TotalBAs (µM)") +
  scale_fill_manual(values = c("#00007F", "#008000")) +
  scale_color_manual(values = c("#00007F", "#008000")) +
  theme_classic() +
  theme(plot.title = element_text(size=18, hjust = 0.5, face="bold"), 
        axis.text = element_text(size=14),
        axis.title = element_text(size=16), 
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10))) +
  guides(fill = FALSE, color = FALSE, alpha=FALSE) 
figure_3_b
dev.off()


# Just this one left to do

#figure_3_D
type_fam<-all_basic[,c("Mother_or_Baby", "Family_ID")]
bile_acids_prop_all<-merge(bile_acids_prop,type_fam, by="row.names" )
bile_acids_prop_all$ID<-paste0(bile_acids_prop_all$Family_ID,"_", bile_acids_prop_all$Mother_or_Baby)
row.names(bile_acids_prop_all)<-bile_acids_prop_all$ID
bile_acids_prop_all$Row.names=NULL
bile_acids_prop_all$Mother_or_Baby=NULL
bile_acids_prop_all$ID=NULL
bile_acids_prop_all$Family_ID=NULL

mother_BA_prop<-bile_acids_prop_all[grep("_Mother P40$", rownames(bile_acids_prop_all)), ]
rownames(mother_BA_prop) <- gsub("_Mother P40", "", rownames(mother_BA_prop))
mother_BA_selection_order <- mother_BA_prop %>%
  arrange(desc(CA+GCA+TCA+DCA+GDCA+TDCA)) %>%
  rownames()
mother_BA_prop$SAMPLE_ID<-row.names(mother_BA_prop)
mother_BA_prop_long<-melt(mother_BA_prop)

mother_BA_prop_long <- mother_BA_prop_long %>%
  dplyr::rename(bile_acid = variable)


bile_acid_levels <- c("CA", "GCA", "TCA", "DCA", "GDCA", "TDCA", 
                     "CDCA", "GCDCA", "TCDCA", "LCA", "GLCA", "TLCA", 
                     "UDCA", "GUDCA", "TUDCA", "LCA3S", "GLCA3S", "TLCA3S")

mother_BA_prop_long$bile_acid <- factor(mother_BA_prop_long$bile_acid, 
                                            levels = bile_acid_levels)

# Define the new colors for each variable
new_colors <- c(
  "CA" = "#130661",
  "GCA" = "#231381",
  "TCA" = "#474184",
  "DCA" = "#6b6c9f",
  "GDCA" = "#8898b8",
  "TDCA" = "#bfd0e0",
  "CDCA" = "#f2ff37",
  "GCDCA" = "#e4cf37",
  "TCDCA" = "#d9a031",
  "LCA" = "#744700",
  "GLCA" = "#5e3f0f",
  "TLCA" = "#715d3e",
  "UDCA" = "#e3b1c2",
  "GUDCA" = "#a25ebf",
  "TUDCA" = "#6c0ecf",
  "LCA3S" = "#c8742a",
  "GLCA3S" = "#c62b1b",
  "TLCA3S" = "#c11d1e"
)


stalked_bar_RA_mother <- ggplot(mother_BA_prop_long, aes(x = SAMPLE_ID, y = value,  group=bile_acid, fill = bile_acid))+
  geom_bar(width = 1, position = "fill", stat = "identity")+
  scale_x_discrete(limits = mother_BA_selection_order)+
  scale_y_continuous(expand = c(0,0)) +
  ylab("Mother P40")+
  theme_bw()+
  theme(legend.position = "right",
        legend.justification = c(0, 1),
        legend.text = element_text(size = 18),
        legend.title = element_blank(),
        plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.line = element_line(colour = 'white'), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = NA),
        panel.border = element_blank(),
        axis.line.x = element_line(colour = 'white'),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        #axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size =4 ),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size= 20))+
  scale_fill_manual(values = new_colors)

stalked_bar_RA_mother


infant_BA_prop<-bile_acids_prop_all[grep("Cord Blood$", rownames(bile_acids_prop_all)), ]
rownames(infant_BA_prop) <- gsub("_Cord Blood", "", rownames(infant_BA_prop))
infant_BA_prop$SAMPLE_ID<-row.names(infant_BA_prop)
infant_BA_prop_long<-melt(infant_BA_prop)
infant_BA_prop_long <- infant_BA_prop_long %>%
  dplyr::rename(bile_acid = variable)

stalked_bar_RA_infant <- ggplot(infant_BA_prop_long, aes(x = SAMPLE_ID, y = value,  group=bile_acid, fill = bile_acid))+
  geom_bar(width = 1, position = "fill", stat = "identity")+
  scale_x_discrete(limits = mother_BA_selection_order)+
  scale_y_continuous(expand = c(0,0)) +
  ylab("Cord Blood")+
  theme_bw()+
  theme(legend.position = "right",
        legend.justification = c(0, 1),
        legend.text = element_text(size = 18),
        legend.title = element_blank(),
        plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.line = element_line(colour = 'white'), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = NA),
        panel.border = element_blank(),
        axis.line.x = element_line(colour = 'white'),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        #axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size =4 ),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=20))+
  scale_fill_manual(values = new_colors)

stalked_bar_RA_infant

pdf('chapter_hilde/Figures/Figure_3_d.pdf', width=10, height=5)
RA_BA_birth<-ggarrange(stalked_bar_RA_mother,stalked_bar_RA_infant, nrow = 2, ncol = 1, common.legend = TRUE)
annotate_figure(RA_BA_birth, top = text_grob("", 
                                                 color = "black", face = "bold", size = 14))

RA_BA_birth
dev.off()

# Figure 3e
pdf('chapter_hilde/Figures/Figure_3_e.pdf', width=6, height=5)
figure_3_e<-ggplot(all_basic, 
                   aes(x=Mother_or_Baby, y=Fraction_secondary, fill = Mother_or_Baby, color=Mother_or_Baby)) +
  geom_violin(aes(fill=Mother_or_Baby, alpha=0.2, width=0.3, trim = F)) +
  geom_boxplot(alpha=0.9, outlier.colour = NA, width=0.03, fill="white") +
  geom_point(alpha=0.6,
             position = position_jitterdodge(jitter.width = 0.3, jitter.height = 0), size=0.2) +
  labs(x = "", y = "Fraction secondary BAs of \n total BAs") +
  scale_fill_manual(values = c("#00007F", "#008000")) +
  scale_color_manual(values = c("#00007F", "#008000")) +
  theme_classic() +
  theme(plot.title = element_text(size=18, hjust = 0.5, face="bold"), 
        axis.text = element_text(size=14),
        axis.title = element_text(size=16), 
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10))) +
  guides(fill = FALSE, color = FALSE, alpha=FALSE) 
figure_3_e
dev.off()


# Figure 3f
pdf('chapter_hilde/Figures/Figure_3_f.pdf', width=6, height=5)
figure_3_f<-ggplot(all_basic, 
                   aes(x=Mother_or_Baby, y=Fraction_conjugated, fill = Mother_or_Baby, color=Mother_or_Baby)) +
  geom_violin(aes(fill=Mother_or_Baby, alpha=0.2, width=0.3, trim = F)) +
  geom_boxplot(alpha=0.9, outlier.colour = NA, width=0.03, fill="white") +
  geom_point(alpha=0.6,
             position = position_jitterdodge(jitter.width = 0.3, jitter.height = 0), size=0.2) +
  labs(x = "", y = "Fraction conjugated BAs of total BAs") +
  scale_fill_manual(values = c("#00007F", "#008000")) +
  scale_color_manual(values = c("#00007F", "#008000")) +
  theme_classic() +
  theme(plot.title = element_text(size=18, hjust = 0.5, face="bold"), 
        axis.text = element_text(size=14),
        axis.title = element_text(size=16), 
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10))) +
  guides(fill = FALSE, color = FALSE, alpha=FALSE) 
figure_3_f
dev.off()

# Figure 3g
pdf('chapter_hilde/Figures/Figure_3_g.pdf', width=6, height=5)
figure_3_g<-ggplot(all_basic, 
                   aes(x=Mother_or_Baby, y=frac_Taurine, fill = Mother_or_Baby, color=Mother_or_Baby)) +
  geom_violin(aes(fill=Mother_or_Baby, alpha=0.2, width=0.3, trim = F)) +
  geom_boxplot(alpha=0.9, outlier.colour = NA, width=0.03, fill="white") +
  geom_point(alpha=0.6,
             position = position_jitterdodge(jitter.width = 0.3, jitter.height = 0), size=0.2) +
  labs(x = "", y = "Fraction taurine conjugated BAs of \n conjugated BAs") +
  scale_fill_manual(values = c("#00007F", "#008000")) +
  scale_color_manual(values = c("#00007F", "#008000")) +
  theme_classic() +
  theme(plot.title = element_text(size=18, hjust = 0.5, face="bold"), 
        axis.text = element_text(size=14),
        axis.title = element_text(size=16), 
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10))) +
  guides(fill = FALSE, color = FALSE, alpha=FALSE) 
figure_3_g
dev.off()



# Figure 3h
pdf('chapter_hilde/Figures/Figure_3_h.pdf', width=6, height=5)
figure_3_h<-ggplot(all_basic, 
                   aes(x=Mother_or_Baby, y=Ratio_12OH_non12OH, fill = Mother_or_Baby, color=Mother_or_Baby)) +
  geom_violin(aes(fill=Mother_or_Baby, alpha=0.2, width=0.3, trim = F)) +
  geom_boxplot(alpha=0.9, outlier.colour = NA, width=0.03, fill="white") +
  geom_point(alpha=0.6,
             position = position_jitterdodge(jitter.width = 0.3, jitter.height = 0), size=0.2) +
  labs(x = "", y = "Ratio of 12alphaOH/ \n non 12alphaOH BAs") +
  scale_fill_manual(values = c("#00007F", "#008000")) +
  scale_color_manual(values = c("#00007F", "#008000")) +
  theme_classic() +
  theme(plot.title = element_text(size=18, hjust = 0.5, face="bold"), 
        axis.text = element_text(size=14),
        axis.title = element_text(size=16), 
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10))) +
  guides(fill = FALSE, color = FALSE, alpha=FALSE) 
figure_3_h
dev.off()


# Figure 3i
pdf('chapter_hilde/Figures/Figure_3_i.pdf', width=6, height=5)
figure_3_i<-ggplot(all_basic, 
                   aes(x=Mother_or_Baby, y=Ratio_sulfated_nonsulfated, fill = Mother_or_Baby, color=Mother_or_Baby)) +
  geom_violin(aes(fill=Mother_or_Baby, alpha=0.2, width=0.3, trim = F)) +
  geom_boxplot(alpha=0.9, outlier.colour = NA, width=0.03, fill="white") +
  geom_point(alpha=0.6,
             position = position_jitterdodge(jitter.width = 0.3, jitter.height = 0), size=0.2) +
  labs(x = "", y = "Ratio sulfated BAs/non sulfated BAs") +
  scale_fill_manual(values = c("#00007F", "#008000")) +
  scale_color_manual(values = c("#00007F", "#008000")) +
  theme_classic() +
  theme(plot.title = element_text(size=18, hjust = 0.5, face="bold"), 
        axis.text = element_text(size=14),
        axis.title = element_text(size=16), 
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10))) +
  guides(fill = FALSE, color = FALSE, alpha=FALSE) 
figure_3_i
dev.off()





# Extra related versus unrelated analysis 

dist_matrix <- vegdist(bile_acids_prop_all, method = "euclidean")
dist_matrix <-as.matrix(dist_matrix)
dist_matrix [upper.tri(dist_matrix )] <- NA 
distances_long<-melt(dist_matrix)
head (distances_long)

# Remove duplicate distance 
distances_long=distances_long %>%
  drop_na()

# Add related or unrelated on basis of the FAMID
distances_long$relationship <- ifelse(substr(distances_long$Var1, 1, 7) == substr(distances_long$Var2, 1, 7), "related", "unrelated")


# Removing same sample pairs 
distances_long <- distances_long[distances_long$Var1 != distances_long$Var2, ]


# Add column mother or baby for each pair of samples 
distances_long$Var1=as.character(distances_long$Var1)
distances_long$Var2=as.character(distances_long$Var2)
distances_long$type1=sapply(strsplit(distances_long$Var1,"_"),function(x)x[2])
distances_long$type2=sapply(strsplit(distances_long$Var2,"_"),function(x)x[2])

# Pick only mother-baby pairs 
distances_long=distances_long %>%
  filter(type1!=type2)

#values = c("#00007F", "#008000"))+

ggplot(distances_long, aes(x = relationship, y = value, color=relationship)) +
  geom_boxplot()+
  geom_violin(alpha = 0.5) +
  ylab("Eucledian distances") +
  geom_jitter(position = position_jitter(seed = 1, width = 0.2, alpha=0.2)) +
  xlab("Relationship") +
  theme_bw()+
  ggtitle("Sharing of Bile Acid Profile between related mother-infant pairs")


ggplot(distances_long, aes(x=value,fill=relationship)) + 
  geom_density(alpha=0.4)+
  theme_classic()+labs(x="Eucledian distances", y = "Density")+
  ggtitle("Sharing of Bile Acid Profile between mother-infant pairs")
  theme(legend.title=element_blank())

#Sharing_related_unrelated_mother_infant_pairs

  p<-ggplot(distances_long, aes(x=value, fill=relationship, color=relationship)) +
    geom_density(alpha=0.4)+scale_color_manual(values=c("#56B4E9", "#E69F00" ))+scale_fill_manual(values=c("#56B4E9", "#E69F00" ))
  p
  # Add mean lines
  p+geom_vline(data=mu, aes(xintercept=grp.median, color=relationship),
               linetype="dashed")+
    theme_classic()+labs(x="Euclidean distances", y = "Density")+
    ggtitle("Sharing of Bile Acid Profile between mother-infant pairs")+
  theme(legend.title=element_blank())
  

  mu <- ddply(distances_long, "relationship", summarise, grp.median=median(value))

  results_perm_wilcoxen = data.frame()
  
  for(i in 1:1000){
    res=wilcox.test(distances_long$value ~ sample(distances_long$relationship))$p.value
    results_perm_wilcoxen  = rbind(results_perm_wilcoxen  , res)
  }
  
  names (results_perm_wilcoxen) [1]<-"p_values"
  
  ggplot(results_perm_wilcoxen, aes(x=p_values)) + 
    geom_histogram(color="black", fill="white")
  p
  
  res$p.value <- results_perm_wilcoxen$p_values
  sum(res$p.value <- results_perm_wilcoxen$p_values)
  5/1000
  P<1/1000
  
  wilcox.test(distances_long$value ~ distances_long$relationship)
  #p-value = 0.2471

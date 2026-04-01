# Analysis Bile Acids composition structure 
# Author: Trishla Sinha 
# Date: December 21st, 2024
# Last update: 29th November, 2025

setwd("/Users/trishlasinha/Desktop/Bile_acids_NEXT")
library(vegan)
library(ggplot2)
library(reshape2)
library(tidyverse)
library(ggpubr)
library(lmerTest)
library(dplyr)
library(broom.mixed)  # For tidying model outputs


############ Figure 1: Mother bile acid composition according to timepoint ##########

lijst_1<-read.delim("lijst_1_MB_v_P12_phenotypes_TS_NEW_NOV.txt")
lijst_1$Plasma_Timepoint[lijst_1$Plasma_Timepoint == "B"] <- "P40"
lijst_1$Timepoint[lijst_1$Timepoint == "B"] <- "P40"

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

all_ba_p_dis<-vegdist(bile_acids_prop, method = 'canberra', na.rm = T)
all_ba_p_dis_mds<-cmdscale(all_ba_p_dis, k=5, eig = T)
all_ba_p_pcoa <- data.frame(all_ba_p_dis_mds$points)


#figure_1_a

pdf('Figures/Figure_1_a.pdf', width=6, height=5)
figure_1_a<-ggplot(all_basic, 
       aes(x=Timepoint, y=C4_nM, fill = Timepoint, color=Timepoint)) +
  geom_violin(aes(fill=Timepoint), alpha=0.2, width=0.3, trim = F) +
  geom_boxplot(alpha=0.9, outlier.colour = NA, width=0.05, fill="white") +
  geom_point(alpha=0.6,
             position = position_jitterdodge(jitter.width = 0.3, jitter.height = 0), size=0.2) +
  labs(x = "", y = "C4 (nM)", fill = "Timepoint") +
  scale_fill_manual(values = c("#FF0000", "#00007F")) +
  scale_color_manual(values = c("#FF0000", "#00007F")) +
  theme_classic() +
  theme(plot.title = element_text(size=18, hjust = 0.5, face="bold"), 
        axis.text = element_text(size=14),
        axis.title = element_text(size=16), 
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10))) +
  guides(fill = FALSE, color = FALSE) 
figure_1_a
dev.off()

pdf('Figures/Figure_1_b.pdf', width=6, height=5)
p_ba_p_pcoa_timepoint<-ggplot(all_ba_p_pcoa ,aes(X1,X2, color = all_basic_bile_acids_prop$Timepoint))+
  geom_point(size = 2,alpha = 0.5)+
  stat_ellipse(aes(group = all_basic_bile_acids_prop$Timepoint, fill = all_basic_bile_acids_prop$Timepoint, color = all_basic_bile_acids_prop$Timepoint) ,type = "norm",linetype = 2,geom = "polygon",alpha = 0.05,show.legend = F)+
  xlab(paste("PCo1=",round(all_ba_p_dis_mds$eig[1],digits = 2),"%",sep = ""))+
  ylab(paste("PCo2=",round(all_ba_p_dis_mds$eig[2],digits = 2),"%",sep = ""))+
  scale_color_manual(name=NULL, 
                     breaks = c("P12",  "P40"),
                     labels = c("P12          ", "P40"),
                     values = c("#FF0000", "#00007F"))+
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

adonis2(bile_acids_prop~all_basic_bile_acids_prop$Timepoint)
# p=0.001 ***
adonis2(bile_acids_prop~all_basic_bile_acids_prop$infant_misc_sex)


wilcox.test(all_ba_p_pcoa$X1 ~all_basic_bile_acids_prop$Timepoint)
#p-value < 2.2e-16
wilcox.test(all_ba_p_pcoa$X2 ~all_basic_bile_acids_prop$Timepoint)
#p-value = 0.8716

figure_1_b<-p_ba_p_pcoa_timepoint
figure_1_b
dev.off()


# Figure 1c
# First load and process data from LLD (for conjugated )

metadata_LLD<-read.delim("/Users/trishlasinha/Desktop/Bile_acids_NEXT/LLD_300OB/20200801_LLD_300OB_basic_1437samples.tsv")
LLD<-read.delim("/Users/trishlasinha/Desktop/Bile_acids_NEXT/LLD_300OB/20200801_LLD_39BA_1135samples.tsv")
all<-merge(metadata_LLD, LLD, by="row.names")
women_reproductive_age <-all[all$Gender==0 & all$Age<45,]
unconjugated  <- c('CA','CDCA','DCA','LCA')
conjugated    <- c('TCA','GCA','TCDCA','GCDCA','TDCA', 'GDCA', 'TLCA','GLCA') # 'TLCA_3S','GLCA_3S'

# Compute total conjugated and unconjugated bile acids per sample
women_reproductive_age$conjugated_sum <- rowSums(women_reproductive_age[, conjugated], na.rm = TRUE)
women_reproductive_age$unconjugated_sum <- rowSums(women_reproductive_age[, unconjugated], na.rm = TRUE)

# Compute total bile acids
women_reproductive_age$total_BA <- women_reproductive_age$conjugated_sum + women_reproductive_age$unconjugated_sum

# Compute *fraction* of conjugated bile acids to total bile acids
women_reproductive_age$Fraction_conjugated <- women_reproductive_age$conjugated_sum / women_reproductive_age$total_BA
summary (women_reproductive_age$Fraction_conjugated)
women_reproductive_age$Timepoint<-"Non_pregnant_women"

all_women<- rbind(
  women_reproductive_age[, c("Timepoint", "Fraction_conjugated")],
  all_basic[, c("Timepoint", "Fraction_conjugated")]
)

pdf('Figures/Figure_1_c.pdf', width=6, height=5)
figure_1_c<-ggplot(all_women, 
                   aes(x=Timepoint, y=Fraction_conjugated, fill = Timepoint, color=Timepoint)) +
  geom_violin(aes(fill=Timepoint), alpha=0.2, width=0.3, trim = F) +
  geom_boxplot(alpha=0.9, outlier.colour = NA, width=0.1, fill="white") +
  geom_point(alpha=0.6,
             position = position_jitterdodge(jitter.width = 0.3, jitter.height = 0), size=0.2) +
  labs(x = "", y = "Fraction conjugated BAs of total BAs", fill = "Timepoint") +
  scale_fill_manual(values = c("darkgreen", "#FF0000", "#00007F")) +
  scale_color_manual(values = c("darkgreen","#FF0000", "#00007F")) +
  theme_classic() +
  theme(plot.title = element_text(size=18, hjust = 0.5, face="bold"), 
        axis.text = element_text(size=14),
        axis.title = element_text(size=16), 
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10))) +
  guides(fill = FALSE, color = FALSE)
figure_1_c
dev.off()



#Figure 1d
mother_P12_BA_prop<-bile_acids_prop[grep("_P12$", rownames(bile_acids_prop)), ]
rownames(mother_P12_BA_prop) <- gsub("_P12", "", rownames(mother_P12_BA_prop))
rows_to_remove <- c("LLNEXT000866", "LLNEXT000109", "LLNEXT010126") # Removing samples that do not have a corresponding pair at P28
mother_P12_BA_prop<- mother_P12_BA_prop[!rownames(mother_P12_BA_prop) %in% rows_to_remove, ]
mother_BA_P12_selection_order <- mother_P12_BA_prop %>%
  arrange(desc(CA+GCA+TCA+DCA+GDCA+TDCA)) %>%
  rownames()
mother_P12_BA_prop$SAMPLE_ID<-row.names(mother_P12_BA_prop)
mother_P12_BA_prop_long<-melt(mother_P12_BA_prop)
mother_P12_BA_prop_long <- mother_P12_BA_prop_long %>%
  dplyr::rename(bile_acid = variable)


bile_acid_levels <- c("CA", "GCA", "TCA", "DCA", "GDCA", "TDCA", 
                      "CDCA", "GCDCA", "TCDCA", "LCA", "GLCA", "TLCA", 
                      "UDCA", "GUDCA", "TUDCA", "LCA3S", "GLCA3S", "TLCA3S")

mother_P12_BA_prop_long$bile_acid <- factor(mother_P12_BA_prop_long$bile_acid, 
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


# Remove ID's not present in the next one 
stalked_bar_RA_P12 <- ggplot(mother_P12_BA_prop_long, aes(x = SAMPLE_ID, y = value,  group=bile_acid, fill = bile_acid))+
  geom_bar(width = 1, position = "fill", stat = "identity")+
  scale_x_discrete(limits = mother_BA_P12_selection_order)+
  scale_y_continuous(expand = c(0,0)) +
  ylab("P12")+
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


stalked_bar_RA_P12



mother_B_BA_prop<-bile_acids_prop[grep("_B$", rownames(bile_acids_prop)), ]
rownames(mother_B_BA_prop) <- gsub("_B", "", rownames(mother_B_BA_prop))
mother_B_BA_prop$SAMPLE_ID<-row.names(mother_B_BA_prop)
mother_B_BA_prop_long <- melt(mother_B_BA_prop)
mother_B_BA_prop_long <- mother_B_BA_prop_long %>%
  rename(bile_acid = variable)


mother_B_BA_prop_long$bile_acid <- factor(mother_B_BA_prop_long$bile_acid, 
                                            levels = bile_acid_levels)

stalked_bar_RA_B <- ggplot(mother_B_BA_prop_long, aes(x = SAMPLE_ID, y = value,  group=bile_acid, fill = bile_acid))+
  geom_bar(width = 1, position = "fill", stat = "identity")+
  #scale_x_discrete(limits = mother_BA_P12_selection_order)+ #Note: here using the order of P12
  scale_x_discrete(limits = mother_BA_P12_selection_order)+
  scale_y_continuous(expand = c(0,0)) +
  ylab("P40")+
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
 
stalked_bar_RA_B

pdf('Figures/Figure_1_d.pdf', width=10, height=4)
RA_BA_pregnancy<-ggarrange(stalked_bar_RA_P12,stalked_bar_RA_B, nrow = 2, ncol = 1, common.legend = TRUE)
annotate_figure(RA_BA_pregnancy, top = text_grob("", 
                                 color = "black", face = "bold", size = 20))
figure_1_d<-RA_BA_pregnancy
figure_1_d
dev.off()



# Figure 1e
pdf('Figures/Figure_1_e.pdf', width=6, height=5)
figure_1_e<-ggplot(all_basic, 
       aes(x=Timepoint, y=Taurine_frac, fill = Timepoint, color=Timepoint)) +
  geom_violin(aes(fill=Timepoint), alpha=0.2, width=0.3, trim = F) +
  geom_boxplot(alpha=0.9, outlier.colour = NA, width=0.1, fill="white") +
  geom_point(alpha=0.6,
             position = position_jitterdodge(jitter.width = 0.3, jitter.height = 0), size=0.2) +
  labs(x = "", y = "Fraction taurine conjugated BAs of \n conjugated BAs", fill = "Timepoint") +
  
  scale_fill_manual(values = c("#FF0000", "#00007F")) +
  scale_color_manual(values = c("#FF0000", "#00007F")) +
  theme_classic() +
  theme(plot.title = element_text(size=18, hjust = 0.5, face="bold"), 
        axis.text = element_text(size=14),
        axis.title = element_text(size=16), 
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10))) +
  guides(fill = FALSE, color = FALSE)
figure_1_e
dev.off()

# Figure 1f
pdf('Figures/Figure_1_f.pdf', width=6, height=5)
figure_1_f<-ggplot(all_basic, 
       aes(x=Timepoint, y=Ratio_12OH_non12OH, fill = Timepoint, color=Timepoint)) +
  geom_violin(aes(fill=Timepoint), alpha=0.2, width=0.3, trim = F) +
  geom_boxplot(alpha=0.9, outlier.colour = NA, width=0.1, fill="white") +
  geom_point(alpha=0.6,
             position = position_jitterdodge(jitter.width = 0.3, jitter.height = 0), size=0.2) +
  labs(x = "", y = "Ratio of 12alphaOH/ \n non 12alphaOH BAs", fill = "Timepoint") +
  scale_fill_manual(values = c("#FF0000", "#00007F")) +
  scale_color_manual(values = c("#FF0000", "#00007F")) +
  theme_classic() +
  theme(plot.title = element_text(size=18, hjust = 0.5, face="bold"), 
        axis.text = element_text(size=14),
        axis.title = element_text(size=16), 
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10))) +
  guides(fill = FALSE, color = FALSE)
figure_1_f
dev.off()

# Figure 1g: Ratio_sulfated_nonsulfated
pdf('Figures/Figure_1_g.pdf', width=6, height=5)
figure_1_g<-ggplot(all_basic, 
       aes(x=Timepoint, y=Ratio_sulfated_nonsulfated, fill = Timepoint, color=Timepoint)) +
  geom_violin(aes(fill=Timepoint), alpha=0.2, width=0.3, trim = F) +
  geom_boxplot(alpha=0.9, outlier.colour = NA, width=0.05, fill="white") +
  geom_point(alpha=0.6,
             position = position_jitterdodge(jitter.width = 0.3, jitter.height = 0), size=0.2) +
  labs(x = "", y = "Ratio sulfated BAs/non sulfated BAs", fill = "Timepoint") +
  scale_fill_manual(values = c("#FF0000", "#00007F")) +
  scale_color_manual(values = c("#FF0000", "#00007F")) +
  theme_classic() +
  theme(plot.title = element_text(size=18, hjust = 0.5, face="bold"), 
        axis.text = element_text(size=14),
        axis.title = element_text(size=16), 
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10))) +
  guides(fill = FALSE, color = FALSE)
figure_1_g
dev.off()



pdf('Figures/Figure_1_e_f_g.pdf', width=14, height=5)
figure_1_e_f_g<-ggarrange(figure_1_e, figure_1_f,figure_1_g,
                    labels = c( "e", "f","g" ), 
                    ncol = 3, nrow = 1)

figure_1_e_f_g
dev.off()




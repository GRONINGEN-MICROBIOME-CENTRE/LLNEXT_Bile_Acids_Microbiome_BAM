# Phenotype associations 

# ------------------- Load libraries ------------------- #
library(tidyverse)
library(reshape2)
library(dplyr)
library(pheatmap)

# ------------------- Drop phenotypes function ------------------- #
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
phenos_test<-c(
  "mother_lifelinesbirthcards_prepreg_bmi_kg_m2",  
  "mother_education_p18"  ,
  "mother_income_net_p18"  ,
  "mother_health_smoked_one_whole_year_p18" ,   
  "mother_food_avoid_gluten_preg" ,
  "mother_food_avoid_dairy_preg"   ,
  "birth_deliverybirthcard_mode_binary",
  "mother_birthcardself_gestational_age_weeks",
  "mother_birthcard_parity",
  "mother_birthcardhealth_gravidity",
  "mother_deliverybirthcard_med_del_epidural",
  "birth_deliverybirthcard_place_delivery_simple",
  "mother_deliveryhealth_preg_complaint_flu_rti",
  "mother_deliveryhealth_preg_complaint_gastroenteritis",
  "mother_deliveryhealth_preg_complaint_uti",
  "mother_deliveryhealth_preg_complaint_vaginal_fungal_inf",
  "mother_deliverybirthcard_preg_comp_hypertension",
  "mother_deliverybirthcard_preg_comp_GDM",
  "mother_deliverybirthcardhealth_preg_complaint_hyperemesis_gravidarum",
  "infant_deliverybirthcard_meconium_amniotic_fluid",
  "mother_birthcard_age_at_delivery",
  "infant_misc_sex",
  "family_pets_any",
  "infant_ffq_ever_never_breastfed",
  "infant_birthcard_apgar_score_1min" ,
  "infant_growth_birth_weight_kg",
  "infant_growth_standardized_weight_slope_kg",
  "mother_ffq_aMED_score"   ,
  "mother_ffq_energy_kcal_m2"     ,
  "mother_ffq_pattern_vegetarian_foods"     ,
  "infant_health_food_allergy"        
)

# ------------------- Linear model function ------------------- #
linear_model_taxa_simple_per_timepoint <- function(metadata, ID, CLR_transformed_data, pheno_list, tp = "W2") {
  metadata %>% filter(Timepoint == tp) -> metadata
  df <- merge(metadata, CLR_transformed_data, by='row.names')
  row.names(df) <- df$Row.names
  df$Row.names <- NULL
  Prevalent <- colnames(CLR_transformed_data)
  
  Overall_result_phenos <- tibble()
  
  for (Bug in Prevalent) {
    if (!Bug %in% colnames(df)) next
    Bug2 <- paste0("`", Bug, "`")
    
    for (pheno in pheno_list) {
      message(paste("Testing outcome:", Bug, "| Phenotype:", pheno, "| Timepoint:", tp))
      pheno2 <- paste0("`", pheno, "`")
      To_keep <- df[!is.na(df[[pheno]]), ID]
      df_pheno <- filter(df, !!sym(ID) %in% To_keep)
      Model2 <- as.formula(paste(Bug2, "~", pheno2))
      resultmodel2 <- lm(Model2, df_pheno)
      
      Summ_simple <- as.data.frame(summary(resultmodel2)$coefficients)[
        grep(pheno, row.names(as.data.frame(summary(resultmodel2)$coefficients))), ]
      
      temp_output <- Summ_simple %>% 
        rownames_to_column("Feature") %>% 
        as_tibble() %>% 
        mutate(P = Summ_simple$`Pr(>|t|)`, Outcome = Bug, Timepoint = tp, Correction = "technical")
      
      Overall_result_phenos <- bind_rows(Overall_result_phenos, temp_output)
    }
  }
  
  return(as.data.frame(Overall_result_phenos))
}

# ------------------- Set working directory ------------------- #
setwd("/Users/trishlasinha/Desktop/Bile_acids_NEXT")

# ------------------- Load main dataset ------------------- #
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
all_basic_bile_acids_prop <- na.omit(all_basic_bile_acids_prop)

# ------------------- Subset by timepoint ------------------- #
all_basic_bile_acids_prop_P12 <- all_basic_bile_acids_prop[all_basic_bile_acids_prop$Timepoint=="P12", ]
all_basic_bile_acids_prop_P40 <- all_basic_bile_acids_prop[all_basic_bile_acids_prop$Timepoint=="P40", ]
all_basic_bile_acids_prop_infants <- all_basic_bile_acids_prop[all_basic_bile_acids_prop$Timepoint=="B", ]

# ------------------- Load phenotype file ------------------- #
cross_phenotypes <- read.delim("~/Desktop/LLNEXT/Analysis/phenotypes/masterfile_cross_sectional_2023_11_15.txt")

# ------------------- P12 processing ------------------- #
names(all_basic_bile_acids_prop_P12)[24] <- "next_id_mother"
cross_phenotypes_P12 <- cross_phenotypes %>%
  filter(next_id_mother %in% all_basic_bile_acids_prop_P12$next_id_mother)
cross_phenotypes_P12 <- cross_phenotypes_P12[!duplicated(cross_phenotypes_P12$next_id_mother), ]
row.names(cross_phenotypes_P12) <- cross_phenotypes_P12$next_id_mother
cross_phenotypes_P12$Timepoint="P12"
cross_phenotypes_P12$next_id_infant=NULL
row.names(all_basic_bile_acids_prop_P12)<-all_basic_bile_acids_prop_P12$next_id_mother
all_basic_bile_acids_prop_P12$next_id_mother=NULL
all_basic_bile_acids_prop_P12$Timepoint=NULL

# ------------------- P12 association ------------------- #
association_P12 <- linear_model_taxa_simple_per_timepoint(
  metadata = cross_phenotypes_P12,
  ID = "next_id_mother",
  CLR_transformed_data = all_basic_bile_acids_prop_P12,
  pheno_list = phenos_test,
  tp = "P12"
)

association_P12$FDR<-p.adjust(association_P12$P, method = "BH")
write.table(association_P12, "Submission_2026/association_phenotypes_BA_P12.txt", row.names = F, sep = "\t")
#all_P12<-merge(cross_phenotypes_P12, all_basic_bile_acids_prop_P12, by="row.names")

# ------------------- P40 processing ------------------- #
names(all_basic_bile_acids_prop_P40)[24] <- "next_id_mother"

cross_phenotypes_P40 <- cross_phenotypes %>%
  filter(next_id_mother %in% all_basic_bile_acids_prop_P40$next_id_mother)
cross_phenotypes_P40 <- cross_phenotypes_P40[!duplicated(cross_phenotypes_P40$next_id_mother), ]
row.names(cross_phenotypes_P40) <- cross_phenotypes_P40$next_id_mother
cross_phenotypes_P40$Timepoint <- "P40"
cross_phenotypes_P40$next_id_infant = NULL
cross_phenotypes_P40$next_id_partner=NULL

row.names(all_basic_bile_acids_prop_P40) <- all_basic_bile_acids_prop_P40$next_id_mother
all_basic_bile_acids_prop_P40$next_id_mother <- NULL
all_basic_bile_acids_prop_P40$Timepoint <- NULL

# ------------------- P40 association & figure ------------------- #
association_P40 <- linear_model_taxa_simple_per_timepoint(
  metadata = cross_phenotypes_P40,
  ID = "next_id_mother",
  CLR_transformed_data = all_basic_bile_acids_prop_P40,
  pheno_list = phenos_test,
  tp = "P40"
)

association_P40$FDR <- p.adjust(association_P40$P, method = "BH")
write.table(association_P40, "Submission_2026/association_phenotypes_BA_P40.txt", row.names = F, sep = "\t")
all_P40<-merge(cross_phenotypes_P40, all_basic_bile_acids_prop_P40, by="row.names")

TUDCA_gestation_P40 <- ggplot(
  all_P40,
  aes(
    x = mother_birthcardself_gestational_age_weeks,
    y = TUDCA
  )
) +
  geom_point(
    alpha = 0.7,
    size = 2,
    color = "#1B9E77"
  ) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    color = "black",
    linewidth = 1
  ) +
  labs(
    title = "",
    x = "Gestational Age (weeks)",
    y = "TUDCA"
  ) +
  theme_minimal(base_size = 14)

TUDCA_gestation_P40

# ------------------- Infant cord blood processing ------------------- #
names(all_basic_bile_acids_prop_infants)[24] <- "next_id_infant"
cross_phenotypes_infant <- cross_phenotypes %>%
  filter(next_id_infant %in% all_basic_bile_acids_prop_infants$next_id_infant)

row.names(cross_phenotypes_infant) <- cross_phenotypes_infant$next_id_infant
cross_phenotypes_infant$Timepoint="B"
cross_phenotypes_infant$next_id_mother=NULL

row.names(all_basic_bile_acids_prop_infants)<-all_basic_bile_acids_prop_infants$next_id_infant
all_basic_bile_acids_prop_infants$next_id_infant=NULL
all_basic_bile_acids_prop_infants$Timepoint=NULL

all_basic_bile_acids_prop_infants <-
  all_basic_bile_acids_prop_infants[
    , !colnames(all_basic_bile_acids_prop_infants) %in%
      c("TUDCA", "LCA", "TLCA", "GLCA", "LCA_3S",
        "GLCA_3S_quant", "TLCA_3S_qual")
  ]
# Removed "UDCA"  "TUDCA" "LCA"   "TLCA"  "GLCA"  "LCA3S" due to a prevalence of non-zero values <10%
# ------------------- Infant cord blood association ------------------- #
association_infant_cord_blood <- linear_model_taxa_simple_per_timepoint(
  metadata = cross_phenotypes_infant,
  ID = "next_id_infant",
  CLR_transformed_data = all_basic_bile_acids_prop_infants,
  pheno_list = phenos_test,
  tp = "B"
)

association_infant_cord_blood$FDR<-p.adjust(association_infant_cord_blood$P, method = "BH")
write.table(association_infant_cord_blood, "Submission_2026/association_infant_cord_blood.txt", row.names = F, sep = "\t")

all_infants<-merge(cross_phenotypes_infant, all_basic_bile_acids_prop_infants, by="row.names")

TCDCA_gestation_infants <- ggplot(
  all_infants,
  aes(
    x = mother_birthcardself_gestational_age_weeks,
    y = TCDCA
  )
) +
  geom_point(
    alpha = 0.7,
    size = 2,
    color = "#1B9E77"
  ) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    color = "black",
    linewidth = 1
  ) +
  labs(
    title = "",
    x = "Gestational Age (weeks)",
    y = "TCDCA"
  ) +
  theme_minimal(base_size = 14)

TCDCA_gestation_infants


ggplot(
  all_infants[ 
    !is.na(all_infants$mother_deliverybirthcard_preg_comp_GDM) &
      !is.na(all_infants$TUDCA),
  ],
  aes(
    x = mother_deliverybirthcard_preg_comp_GDM,
    y = TUDCA,
    fill  = mother_deliverybirthcard_preg_comp_GDM,
    color = mother_deliverybirthcard_preg_comp_GDM
  )
) +
  geom_boxplot(alpha = 0.9, outlier.colour = NA, width = 0.03, fill = "white") +
  geom_point(
    alpha = 0.6,
    position = position_jitter(width = 0.25, height = 0),
    size = 0.4
  ) +
  labs(x = "", y = "LCA") +
  scale_fill_manual(values = c("#FF0000", "#00007F")) +
  scale_color_manual(values = c("#FF0000", "#00007F")) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.title.x = element_text(margin = margin(t = 10))
  ) +
  guides(fill = FALSE, color = FALSE, alpha = FALSE)


# Makign summary stats table of the data for table S1 and S2 
all_basic_summary<-all_basic_bile_acids_prop[!all_basic_bile_acids_prop$Timepoint=="P28",]
cross_phenotypes_summary <- cross_phenotypes[
  cross_phenotypes$next_id_mother %in% all_basic_summary$NEXT_ID,
]

cross_phenotypes_summary <-cross_phenotypes_summary[,c("next_id_mother",phenos_test)]

cross_phenotypes_summary <- cross_phenotypes_summary[
  !duplicated(cross_phenotypes_summary$next_id_mother),
]

cross_phenotypes_summary <- cross_phenotypes_summary %>%
  mutate(across(where(is.character), as.factor))

setwd("/Users/trishlasinha/Desktop/Bile_acids_NEXT/Submission_2026")
generate_summary_statistics <- function(df) {
  
  #subset numeric variables
  numeric_vars <- df %>%
    select_if(is.numeric)
  
  #generate summary statistics for numeric variables
  numeric_summary <- numeric_vars %>%
    pivot_longer(everything(), names_to = "variable", values_to = "value") %>%
    group_by(variable) %>%
    summarise(
      total_n = sum(!is.na(value)),
      percentage_of_total_n = sum(!is.na(value)) / nrow(df) * 100,
      total_number_missings = sum(is.na(value)),
      percentage_of_total_number_missings = sum(is.na(value)) / nrow(df) * 100,
      mean = mean(value, na.rm = TRUE),
      sd = sd(value, na.rm = TRUE),
      min = min(value, na.rm = TRUE),
      Q1 = quantile(value, 0.25, na.rm = TRUE),
      median = median(value, na.rm = TRUE),
      Q3 = quantile(value, 0.75, na.rm = TRUE),
      max = max(value, na.rm = TRUE)
      
    ) %>%
    ungroup()
  
  #subset categorical variables
  factor_vars <- df %>%
    select_if(is.factor)
  
  #generate summary statistics for categorical variables
  factor_summary <- factor_vars %>%
    pivot_longer(cols = everything(),
                 names_to = "variable",
                 values_to = "value") %>%
    group_by(variable) %>%
    summarise(
      n_factor_levels = n_distinct(value, na.rm = TRUE),
      total_n = sum(!is.na(value)),
      percentage_total_n = sum(!is.na(value)) / nrow(df) * 100,
      total_number_missings = sum(is.na(value)),
      percentage_of_total_number_missings = sum(is.na(value)) / nrow(df) * 100)
  factor_summary_levels <- #this script was updated here for factor_summary_levels and how it was defined
    factor_vars %>%
    pivot_longer(cols = everything(),
                 names_to = "variable",
                 values_to = "value") %>%
    group_by(variable, value) %>%
    summarise(
      n = n(),
      perc = n() / nrow(factor_vars) * 100
    ) %>%
    drop_na()
  factor_summary_levels <- factor_summary_levels %>%
    mutate(perc=n/sum(n)*100) %>%
    mutate(
      level = paste0("level_", row_number())
    ) %>%
    filter(value!="NA") %>%
    pivot_wider(names_from = level,
                values_from = c(value, n, perc),
                names_glue = "{level}_{.value}")
  
  #determine the maximum number of levels in the factor_summary_levels (for colname reordering)
  level_cols <- grep("^level_\\d+_value$", names(factor_summary_levels), value = TRUE)
  n <- max(as.integer(gsub("^level_(\\d+)_value$", "\\1", level_cols)))
  
  #define order of colnames for categorical variables
  colnames_order <- c("variable",
                      unlist(lapply(1:n, function(i) paste0("level_", i, "_", c("value", "n", "perc")))))
  
  #reorder the columns of the dataframe
  factor_summary_levels <- factor_summary_levels[, colnames_order]
  
  #combine summary statistics data for factor summary
  factor_summary <- left_join(factor_summary,factor_summary_levels, by = "variable") %>%
    ungroup()
  
  #export the new table
  write.table (numeric_summary, file = "./meta_data_numeric_summary_stats.txt" , quote = F, sep = "\t", row.names = F)
  write.table (factor_summary, file = "./meta_data_categorical_summary_stats.txt" , quote = F, sep = "\t", row.names = F)
}

generate_summary_statistics(cross_phenotypes_summary)


# create survival curves for Multi-modal clusters, RNA and methyl-derived subtypes

suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(survival)
  library(survminer)
  library(data.table)
  library(gtsummary)
})

# parse command line options
option_list <- list(
  make_option(c("--cluster_file"), type = "character", help = "File with multi-modal derived clusters (.tsv)"),
  make_option(c("--histology_file"), type = "character", help = "Histology file (.tsv)"),
  make_option(c("--plots_dir"), type = "character", help = "path to plots directory"),
  make_option(c("--output_dir"), type = "character", help = "path to output directory")
)
opt <- parse_args(OptionParser(option_list = option_list, add_help_option = TRUE))

# plots directory
plots_dir <- opt$plots_dir
dir.create(plots_dir, showWarnings = F, recursive = T)

# output directory
output_dir <- opt$output_dir
dir.create(output_dir, showWarnings = F, recursive = T)

# multi-modal derived clusters
mm_clusters <- read_tsv(file = opt$cluster_file)

# combine Multi-modal clusters with RNA-derived molecular subtypes
histology_file <- opt$histology_file
anno_file_rna <- read_tsv(file = histology_file) %>%
  dplyr::select(
    Kids_First_Participant_ID,
    Kids_First_Biospecimen_ID,
    OS_days,
    EFS_days,
    EFS_event_type,
    molecular_subtype,
    extent_of_tumor_resection,
    age_at_diagnosis_days
  ) %>%
  dplyr::rename("Kids_First_Biospecimen_ID_RNA" = "Kids_First_Biospecimen_ID") %>%
  unique() %>%
  dplyr::inner_join(mm_clusters, by = "Kids_First_Biospecimen_ID_RNA")

# code a EFS_status column with the following logic:
anno_file_rna <- anno_file_rna %>%
  dplyr::mutate(OS_days = as.numeric(OS_days), EFS_days = as.numeric(EFS_days))
anno_file_rna <- anno_file_rna %>%
  dplyr::mutate(
    EFS_status = case_when(
      EFS_event_type == 'Not Applicable' ~ 0,
      EFS_event_type != 'Not Applicable' ~ 1
    )
  )

# combine Multi-modal clusters with methylation-derived subclass
anno_file_methyl <- read_tsv(file = histology_file) %>%
  dplyr::select(
    Kids_First_Participant_ID,
    Kids_First_Biospecimen_ID,
    OS_days,
    EFS_days,
    EFS_event_type,
    dkfz_v11_methylation_subclass,
    dkfz_v12_methylation_subclass
  ) %>%
  dplyr::rename("Kids_First_Biospecimen_ID_Methyl" = "Kids_First_Biospecimen_ID") %>%
  unique() %>%
  dplyr::inner_join(mm_clusters, by = "Kids_First_Biospecimen_ID_Methyl")

# code a EFS_status column with the following logic:
anno_file_methyl <- anno_file_methyl %>%
  dplyr::mutate(OS_days = as.numeric(OS_days), EFS_days = as.numeric(EFS_days))
anno_file_methyl <- anno_file_methyl %>%
  dplyr::mutate(
    EFS_status = case_when(
      EFS_event_type == 'Not Applicable' ~ 0,
      EFS_event_type != 'Not Applicable' ~ 1
    )
  )

# combine both and create one standardized annotation file
anno_file <- anno_file_rna %>%
  dplyr::inner_join(anno_file_methyl)
anno_file$dkfz_v11_methylation_subclass <- gsub(", | ", "_", anno_file$dkfz_v11_methylation_subclass)

# format survival data
surv_data <- anno_file
surv_data <- surv_data %>%
  filter(EFS_days <= 4000)

# custom table theme
custom_table_theme <- theme_survminer(base_size = 12) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())

# survival curves stratified by RNA-derived molecular subtype
surv_data$molecular_subtype <- factor(surv_data$molecular_subtype, levels = sort(unique(surv_data$molecular_subtype)))
fit <- survival::survfit(formula = Surv(as.numeric(EFS_days), EFS_status) ~ molecular_subtype,
                         data = surv_data)
pdf(
  file = file.path(plots_dir, "survival_molsubtype.pdf"),
  height = 8,
  width = 8,
  onefile = FALSE
)
p <- ggsurvplot(
  fit,
  title = "Survival stratified by RNA-derived molecular subtypes",
  data = surv_data,
  font.x = c(12),
  font.tickslab = c(12),
  font.y = c(12),
  pval = TRUE,
  ggtheme = ggpubr::theme_pubr(),
  tables.theme = custom_table_theme,
  legend.lab = gsub(".*=", "", summary(fit)$table %>% rownames()),
  ylab = "Event-Free survival probability",
  legend = "none",
  risk.table = TRUE,
  break.x.by = 500
) %++%
  guides(colour = guide_legend(ncol = 1))
print(p)
dev.off()

# survival curves stratified by Methylation-derived molecular subtype (v11)
surv_data$dkfz_v11_methylation_subclass <- factor(surv_data$dkfz_v11_methylation_subclass, levels = sort(unique(
  surv_data$dkfz_v11_methylation_subclass
)))

surv_data_temp <- surv_data %>%
  dplyr::group_by(dkfz_v11_methylation_subclass) %>%
  dplyr::mutate(n = n()) %>%
  dplyr::filter(n >= 5) %>%
  dplyr::mutate(dkfz_v11_methylation_subclass = factor(dkfz_v11_methylation_subclass
  ))

fit <- survival::survfit(Surv(as.numeric(EFS_days), EFS_status) ~ dkfz_v11_methylation_subclass,
                         data = surv_data_temp)
pdf(
  file = file.path(plots_dir, "survival_methyl_subtype_v11.pdf"),
  height = 8,
  width = 8,
  onefile = FALSE
)
p <- ggsurvplot(
  fit,
  title = "Survival stratified by Methylation-derived molecular subtypes",
  data = surv_data,
  font.x = c(12),
  font.tickslab = c(12),
  font.y = c(12),
  pval = TRUE,
  pval.coord = c(2000, 0.25),
  ggtheme = ggpubr::theme_pubr(),
  tables.theme = custom_table_theme,
  legend.lab = gsub(".*=", "", summary(fit)$table %>% rownames()),
  ylab = "Event-Free survival probability",
  legend = "none",
  risk.table = TRUE,
  break.x.by = 500
) %++%
  guides(colour = guide_legend(ncol = 1))
print(p)
dev.off()

# survival curves stratified by Methylation-derived molecular subtype (v12)
surv_data$dkfz_v12_methylation_subclass <- factor(surv_data$dkfz_v12_methylation_subclass, levels = sort(unique(
  surv_data$dkfz_v12_methylation_subclass
)))
surv_data_temp <- surv_data %>%
  dplyr::group_by(dkfz_v12_methylation_subclass) %>%
  dplyr::mutate(n = n()) %>%
  dplyr::filter(n >= 5) %>%
  dplyr::mutate(dkfz_v12_methylation_subclass = factor(dkfz_v12_methylation_subclass
  ))
fit <- survival::survfit(Surv(as.numeric(EFS_days), EFS_status) ~ dkfz_v12_methylation_subclass,
                         data = surv_data_temp)
pdf(
  file = file.path(plots_dir, "survival_methyl_subtype_v12.pdf"),
  height = 9,
  width = 7,
  onefile = FALSE
)
p <- ggsurvplot(
  fit,
  title = "Survival stratified by Methylation-derived molecular subtypes",
  data = surv_data,
  font.x = c(12),
  font.tickslab = c(12),
  font.y = c(12),
  pval = TRUE,
  pval.coord = c(2000, 0.25),
  ggtheme = ggpubr::theme_pubr(),
  tables.theme = custom_table_theme,
  legend.lab = gsub(".*=", "", summary(fit)$table %>% rownames()),
  ylab = "Event-Free survival probability",
  legend = "none",
  risk.table = TRUE,
  risk.table.height = 0.4,
  break.x.by = 500
) %++%
  guides(colour = guide_legend(ncol = 1))
print(p)
dev.off()

# survival curves stratified by Multi-modal subtypes
surv_data$mm_cluster <- factor(surv_data$mm_cluster, levels = sort(as.numeric(unique(
  surv_data$mm_cluster
))))
fit <- survival::survfit(Surv(as.numeric(EFS_days), EFS_status) ~ mm_cluster, data = surv_data)
pdf(
  file = file.path(plots_dir, "survival_mm_clusters.pdf"),
  height = 9,
  width = 7,
  onefile = FALSE
)
p <- ggsurvplot(
  fit,
  palette = c(
    "#d227da",
    "#beaed4",
    "#ff1493",
    "#abcd21",
    "#8a2be2",
    "#ffe135",
    "#228b22",
    "#967117",
    "#21abcd",
    "#3a56ca",
    "#162e95",
    "#ff8c00",
    "#ff4040",
    "#8b0000"
  ),
  title = "Survival stratified by Multi-modal derived clusters",
  data = surv_data,
  font.x = c(12),
  font.tickslab = c(12),
  font.y = c(12),
  pval = TRUE,
  ggtheme = ggpubr::theme_pubr(),
  tables.theme = custom_table_theme,
  legend.lab = gsub(".*=", "", summary(fit)$table %>% rownames()),
  ylab = "Event-Free survival probability",
  legend = "none",
  risk.table = TRUE,
  risk.table.height = 0.4,
  break.x.by = 500
) %++%
  guides(colour = guide_legend(ncol = 1))
print(p)
dev.off()


# EFS
surv_data_tmp <- anno_file %>%
  dplyr::filter(!is.na(EFS_days)) %>%
  dplyr::group_by(mm_cluster) %>%
  dplyr::mutate(n = n()) %>%
  dplyr::filter(n >= 5) %>%
  dplyr::mutate(mm_cluster = as.factor(mm_cluster)) %>%
  dplyr::mutate(EFS_status = ifelse(EFS_event_type == 'Not Applicable', 0, 1))

res_cox <- survival::coxph(Surv(EFS_days, EFS_status) ~ mm_cluster + 
                             age_at_diagnosis_days, 
                           data = surv_data_tmp)

# write summary

sink(file = file.path(output_dir, 'coxph_summary_EFS.txt'))
  summary(res_cox)
sink()

# save risk scores
fwrite(
  data.table(
    'Kids_First_Participant_ID' = surv_data_tmp$Kids_First_Participant_ID,
    'risk_score' = predict(res_cox, type = 'risk')
  ),
  file.path(output_dir, 'coxph_risk_score_EFS.txt'),
  sep = '\t'
)


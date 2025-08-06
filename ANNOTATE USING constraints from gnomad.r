# ANNOTATE USING constraints from gnomad v4.1

# USE prepossecing_step_function_for_all_groups to process the data


constraint_data <- gnomad.v4.1.constraint_metrics %>%
  filter(constraint_flags == "[]")

# Step 2: Annotate all variant dataframes using transcript info
PTV<- processed_data$PTV %>%
  left_join(constraint_data, by = c("Feature" = "transcript"))

PTV_canonic<- processed_data$PTV_canonic %>%
  left_join(constraint_data, by = c("Feature" = "transcript"))

Missense<- processed_data$Missense %>%
  left_join(constraint_data, by = c("Feature" = "transcript"))

Missense_canonical<- processed_data$Missense_canonical %>%
  left_join(constraint_data, by = c("Feature" = "transcript"))

MPC_only<- processed_data$MPC_only %>%
  left_join(constraint_data, by = c("Feature" = "transcript"))


# For PTVs (filter by lof.oe_ci.upper < 0.6)
PTV_lof <- PTV %>% filter(lof.oe_ci.upper < 0.6)
PTV_canonic_lof <- PTV_canonic %>% filter(lof.oe_ci.upper < 0.6)

# For Missense variants (filter by both lof.oe_ci.upper and mis.oe_ci.upper)

# lof constraint (across all for completeness)
Missense_lof <- Missense %>% filter(lof.oe_ci.upper < 0.6)
Missense_canonical_lof <- Missense_canonical %>% filter(lof.oe_ci.upper < 0.6)
MPC_only_lof <- MPC_only %>% filter(lof.oe_ci.upper < 0.6)

# missense constraint (specific to missense relevance)
Missense_mis <- Missense %>% filter(mis.oe_ci.upper < 0.6)
Missense_canonical_mis <- Missense_canonical %>% filter(mis.oe_ci.upper < 0.6)
MPC_only_mis <- MPC_only %>% filter(mis.oe_ci.upper < 0.6)



Missense_canonical <-Missense_canonical_lof
Missense <-Missense_lof
PTV <- PTV_lof
PTV_canonic <- PTV_canonic_lof
MPC_only <- MPC_only_lof

PTV <- PTV_lof
PTV_canonic <- PTV_canonic_lof
Missense <- Missense_mis
Missense_canonical <- Missense_canonical_mis
MPC_only <- MPC_only_mis
setwd("/home/rachele/SNVs/results_pasteur_tsv_With_FEATURE/Constraint/group2/not_private")
setwd("/home/rachele/SNVs/results_pasteur_tsv_With_FEATURE/Constraint/group2/private")


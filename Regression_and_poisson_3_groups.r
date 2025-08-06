#REGRESSIUON AND POISSON GROUP3




library(readr)
library(data.table)
# library(tidyverse)
library(dplyr)
library(gprofiler2)
library(httr)
library(readxl)
library(tidyr)
library(stringr)
library(dplyr)
library(purrr)
library(broom)
library(gridExtra)
library(ggplot2)
library(ggsignif)  # For adding significance annotations
library(patchwork)
library(FSA)
library(dunn.test)
library(pscl)
library(ggrepel)

#  Load my datasets 
SCHEMA <- read_csv("/home/rachele/Documents/old/geneset_confrontation/SCHEMA.csv")
GWAS_120 <- read.delim("~/GWAS_120.csv")       
BipEx_Bipolar <- read_csv("/home/rachele/Documents/old/geneset_confrontation/Bipolar_Disorder.csv")
SFARI <- read_csv("/home/rachele/Documents/old/geneset_confrontation/SFARI.csv")

# !Clean GENESETS
    # Convert gene names 
    # BipEx_Bipolar
    convert_col <- gconvert(query = BipEx_Bipolar$Gene, organism = "hsapiens", target = "ENSG", mthreshold = Inf, filter_na = TRUE)
    BipEx_Bipolar <- merge(BipEx_Bipolar, convert_col[, c("target", "name")], by.x = "Gene", by.y = "target", all.x = TRUE)
    BipEx_Bipolar <- BipEx_Bipolar %>% select(Gene, name, everything())

    # SCHEMA
    convert_col <- gconvert(query = SCHEMA$Gene, organism = "hsapiens", target = "ENSG", mthreshold = Inf, filter_na = TRUE)
    SCHEMA <- merge(SCHEMA, convert_col[, c("target", "name")], by.x = "Gene", by.y = "target", all.x = TRUE)
    SCHEMA <- SCHEMA %>% select(Gene, name, everything())

    # General clean 
    BipEx_Bipolar <- BipEx_Bipolar %>% subset(select = -Gene) %>% rename('Gene' = name) %>% mutate(Gene = gsub(' ', '', Gene)) %>% setNames(make.names(colnames(.))) 
    BipEx_Bipolar_p_val_PTV <- BipEx_Bipolar[!is.na(BipEx_Bipolar$PTV.Fisher.p.val) & BipEx_Bipolar$PTV.Fisher.p.val <= 0.05,]
    BipEx_Bipolar_p_val_Missense <- BipEx_Bipolar[!is.na(BipEx_Bipolar$Damaging.Missense.Fisher.p.val) & BipEx_Bipolar$Damaging.Missense.Fisher.p.val <= 0.05,]

    SCHEMA <- SCHEMA %>% subset(select = -Gene) %>% rename('Gene' = name) %>% mutate(Gene = gsub(' ', '', Gene)) %>% setNames(make.names(colnames(.)))  
    SCHEMA_OR <- SCHEMA[SCHEMA$OR..Class.I. >= 1, ]
    SCHEMA_pVal <- SCHEMA[SCHEMA$P.meta <= 0.05, ]
    SFARI <- SFARI %>% subset(select = -c(status, `ensembl-id`)) %>% rename("Gene" = `gene-symbol`, "gene_name" = `gene-name`) %>% mutate(Gene = gsub(' ', '', Gene))
    SFARI_syndromic <- SFARI 
    SFARI_non_syndromic_lower_3 <- SFARI[SFARI$`gene-score` < 3, ]
    GWAS_120 <- GWAS_120 %>% rename('Gene' = GENE_name)

    # Unify all the dataframes 
    genes_schema_pval <- unique(SCHEMA_pVal$Gene)
    genes_bipolar <- unique(c(BipEx_Bipolar_p_val_PTV$Gene, BipEx_Bipolar_p_val_Missense$Gene))
    genes_sfari_non_syndromic <- unique(SFARI_non_syndromic_lower_3$Gene)
    genes_gwas <- unique(GWAS_120$Gene)
    genes_schema_or <- unique(SCHEMA_OR$Gene)
    genes_sfari_syndromic <- unique(SFARI_syndromic$Gene)

# LOAD MINE     
# LOAD PERSONAL DATA 
    # Specify the folder containing your TSV files
    folder_path <- "/home/rachele/SNVs/results_pasteur"

    # List all the TSV files in the folder
    tsv_files <- list.files(path = folder_path, pattern = "\\.tsv$", full.names = TRUE)

    # Loop through each file and assign it to an object
    for (file in tsv_files) {
    file_name <- tools::file_path_sans_ext(basename(file))
    assign(file_name, read_tsv(file))
    }

    # Get the names of all dataframes in the environment
    dataframes <- ls(pattern = ".*")

    #! General clean up  
    # Loop through each dataframe and modify the SAMPLES column
    for (df_name in dataframes) {
    df <- get(df_name)
    if ("SAMPLES" %in% names(df)) {
        df$SAMPLES <- gsub("_pool", "", df$SAMPLES)
        assign(df_name, df)
    }
    }

#modify manifest 
manifest_correct$Status <- gsub("FEP-SCZ", "FEP", manifest_correct$Status)
manifest_correct$Status <- gsub("FEP-BD", "FEP", manifest_correct$Status)

# ADD LABELS 
# Function to get the labels (case or Non_Converter) for the outlier values
get_outlier_labels <- function(outlier_value, manifest_df) {
  outlier_values <- strsplit(outlier_value, ",")[[1]]
  labels <- sapply(outlier_values, function(val) {
    if (val %in% manifest_df$Sequencing_number) {
      return(manifest_df$Status[manifest_df$Sequencing_number == val])
    } else {
      return(NA)
    }
  })
  return(paste(labels, collapse = ", "))
}

# Function to check the label for the outlier values
check_outlier_label <- function(outlier_value, manifest_df) {
    # Split the outlier values by ',' and trim whitespace
    outlier_values <- strsplit(outlier_value, ",")[[1]]
    
    # Retrieve corresponding labels from the manifest dataframe
    labels <- sapply(outlier_values, function(val) {
        # Find the label in the manifest dataframe (either "case" or "Non_Converter")
        if (val %in% manifest_df$Sequencing_number) {
            return(manifest_df$Status[manifest_df$Sequencing_number == val])
        } else {
            return(NA) # Return NA if no match is found
        }
    })

    # If there's more than one label, return the unique ones
    unique_labels <- unique(labels)
    
    # If there's both "case" and "Non_Converter", return "mixed"; otherwise return the label
    if (length(unique_labels) > 1) {
        return("mixed")
    } else if (length(unique_labels) == 1) {
        return(unique_labels)
    } else {
        return(NA) # If no label found
    }
}

# DO the UHR NA Non_Converter

    clean_uh_na <- function(df) {
    # Step 1: Remove rows where sample_label contains "UHR_NA"
    df <- df %>%
        filter(!grepl("UHR_NA", sample_label))%>%
        filter(!grepl("UHR-NA", sample_label))
    
    # Step 2: Remove "UHR_NA" and clean up sample_label_2
    df <- df %>%
        mutate(
        # Remove "UHR_NA" and any surrounding commas/spaces
        sample_label_2 = gsub("\\s*,?\\s*UHR_NA\\s*,?\\s*", ",", sample_label_2),
        # Remove leading/trailing commas and spaces
        sample_label_2 = gsub("^,|,$", "", sample_label_2),
        # Replace multiple commas with a single comma
        sample_label_2 = gsub("\\s*,+\\s*", ",", sample_label_2),
        # Trim leading/trailing spaces
        sample_label_2 = trimws(sample_label_2)
        )

     df <- df %>%
    mutate(
      # First remove all UHR_NA instances
      sample_label_2 = gsub("UHR_NA", "", sample_label_2),
      sample_label_2 = gsub("UHR-NA", "", sample_label_2),
      # Then clean up resulting comma patterns
      sample_label_2 = gsub(",+", ",", sample_label_2),  # Replace multiple commas with single
      sample_label_2 = gsub("^,|,$", "", sample_label_2),  # Remove leading/trailing commas
      sample_label_2 = trimws(sample_label_2),  # Trim whitespace
      # Handle cases where we're left with just a comma
      sample_label_2 = ifelse(sample_label_2 == ",", NA, sample_label_2)
    )
      # Additional check for empty strings
    df$sample_label_2[df$sample_label_2 == ""] <- NA    
    
    return(df)
    }


    derive_group_label <- function(sample_label_2) {
    # Split the sample_label_2 by commas and trim whitespace
    labels <- trimws(strsplit(sample_label_2, ",")[[1]])
    allowed <- c("FEP", "Converter", "Non_Converter")
    
    # Check if all labels are valid
    if (!all(labels %in% allowed)) {
        return(NA)
    }
    
    # Get unique labels
    unique_labels <- unique(labels)
    
    # Determine the group based on unique labels
    if (length(unique_labels) == 1) {
        return(unique_labels)  # All labels are the same
    } else {
        return("mixed")  # Multiple different labels
    }
}

clean_na_samples <- function(df) {
  df <- df %>%
    # Remove rows where sample_label_2 is NA or the string "NA"
    filter(!is.na(sample_label_2) & sample_label_2 != "NA") %>%
    
    # Clean sample_label_2 by removing any lingering "NA" strings and tidy up commas/spaces
    mutate(
      sample_label_2 = gsub("\\bNA\\b", "", sample_label_2),           # Remove standalone 'NA'
      sample_label_2 = gsub(",+", ",", sample_label_2),               # Replace multiple commas with one
      sample_label_2 = gsub("^,|,$", "", sample_label_2),             # Remove leading/trailing commas
      sample_label_2 = trimws(sample_label_2),                        # Trim whitespace
      sample_label_2 = ifelse(sample_label_2 == "", NA, sample_label_2) # Convert empty strings to NA
    ) %>%
    
    # After cleaning, remove rows where sample_label_2 became NA or empty again
    filter(!is.na(sample_label_2))
  
  return(df)
}


MPC_only <- `Rubiu-GRCh38.deepvariant.splitted.norm.vep.merged.filtered_gnomad_mpc2.filtered_rarity_stringent_AF.missense_mpc.vcf.gz` 
Missense <- `Rubiu-GRCh38.deepvariant.splitted.norm.vep.merged.filtered_gnomad_mpc2.filtered_rarity_stringent_AF.missense_mpc_AM.vcf.gz` 
Missense_canonical <- `Rubiu-GRCh38.deepvariant.splitted.norm.vep.merged.filtered_gnomad_mpc2.filtered_rarity_stringent_AF.missense_mpc_AM_CANONIC.vcf.gz`
PTV <- `Rubiu-GRCh38.deepvariant.splitted.norm.vep.merged.filtered_gnomad_mpc2.filtered_rarity_stringent_AF.PTV_HC.vcf.gz`
PTV_canonic <- `Rubiu-GRCh38.deepvariant.splitted.norm.vep.merged.filtered_gnomad_mpc2.filtered_rarity_stringent_AF.PTV_HC.vcf.gz`
rm('Rubiu-GRCh38.deepvariant.splitted.norm.vep.merged.filtered_gnomad_mpc2.filtered_rarity_stringent_AF.missense_mpc_AM.vcf.gz')
rm('Rubiu-GRCh38.deepvariant.splitted.norm.vep.merged.filtered_gnomad_mpc2.filtered_rarity_stringent_AF.missense_mpc_AM_CANONIC.vcf.gz')
rm('Rubiu-GRCh38.deepvariant.splitted.norm.vep.merged.filtered_gnomad_mpc2.filtered_rarity_stringent_AF.missense_mpc.vcf.gz')
rm(list = c("SCHEMA", "SCHEMA_OR", "SCHEMA_pVal", "SFARI", "SFARI_non_syndromic_lower_3", "SFARI_syndromic", "GWAS_120", "BipEx_Bipolar", "BipEx_Bipolar_p_val_Missense", "BipEx_Bipolar_p_val_PTV", "convert_col"))
PTV_canonic <- PTV_canonic %>% filter(CANONICAL == "YES")

# Apply the get_outlier_labels function to the dataframes

# Apply the get_outlier_labels function to the dataframes
    Missense$sample_label_2 <- sapply(Missense$SAMPLES, get_outlier_labels, manifest_correct)
    Missense_canonical$sample_label_2 <- sapply(Missense_canonical$SAMPLES, get_outlier_labels, manifest_correct)
    MPC_only$sample_label_2 <- sapply(MPC_only$SAMPLES, get_outlier_labels, manifest_correct)
    PTV_canonic$sample_label_2 <- sapply(PTV_canonic$SAMPLES, get_outlier_labels, manifest_correct)
    PTV$sample_label_2 <- sapply(PTV$SAMPLES, get_outlier_labels,manifest_correct)

    #other label
    Missense$sample_label <- sapply(Missense$SAMPLES,check_outlier_label,manifest_correct)
    Missense_canonical$sample_label <- sapply(Missense_canonical$SAMPLES,check_outlier_label,manifest_correct)
    MPC_only$sample_label <- sapply(MPC_only$SAMPLES,check_outlier_label,manifest_correct)
    PTV_canonic$sample_label <- sapply(PTV_canonic$SAMPLES,check_outlier_label,manifest_correct)
    PTV$sample_label <- sapply(PTV$SAMPLES,check_outlier_label,manifest_correct)

    # appply clean sample fucntion 
    Missense <- clean_uh_na(Missense)
    Missense_canonical <- clean_uh_na(Missense_canonical)
    MPC_only <- clean_uh_na(MPC_only)
    PTV_canonic <- clean_uh_na(PTV_canonic)
    PTV <- clean_uh_na(PTV)   

    Missense <- clean_na_samples(Missense)
    Missense_canonical <- clean_na_samples(Missense_canonical)
    MPC_only <- clean_na_samples(MPC_only)
    PTV_canonic <- clean_na_samples(PTV_canonic)
    PTV <- clean_na_samples(PTV)




    Missense$sample_label_3 <- sapply(Missense$sample_label_2,derive_group_label)
    Missense_canonical$sample_label_3 <- sapply(Missense_canonical$sample_label_2,derive_group_label)
    MPC_only$sample_label_3 <- sapply(MPC_only$sample_label_2,derive_group_label)
    PTV_canonic$sample_label_3 <- sapply(PTV_canonic$sample_label_2,derive_group_label)
    PTV$sample_label_3 <- sapply(PTV$sample_label_2,derive_group_label)

    Missense <-Missense %>% mutate(count = sapply(strsplit(sample_label_2, ","), length))
    Missense_canonical <-Missense_canonical %>% mutate(count = sapply(strsplit(sample_label_2, ","), length))
    MPC_only <-MPC_only%>% mutate(count = sapply(strsplit(sample_label_2, ","), length))
    PTV_canonic <-PTV_canonic%>% mutate(count = sapply(strsplit(sample_label_2, ","), length))
    PTV <-PTV%>% mutate(count = sapply(strsplit(sample_label_2, ","), length))


# manifest for UHR_NA
    UHR_NA_SAMPLE_IDS <- read.csv("~/UHR_NA_SAMPLE_IDS.csv", sep="")
    #example S36827
 Low_Quality_SAMPLES_AND_UHR_NA <- read.delim("~/Low_Quality_SAMPLES_AND_UHR_NA.tsv")
colnames(manifest_correct) <- c("Status", "sample_id", "Case_control")
#modify manifest 
manifest_correct$Status <- gsub("FEP-SCZ", "FEP", manifest_correct$Status)
manifest_correct$Status <- gsub("FEP-BD", "FEP", manifest_correct$Status)




# function for df preparation 
process_df <- function(df, name, brain_gene_consensus_filtered_consensus_no_pitular, brain_gene_consensus_ntm_consensus_no_pitular, genes_schema_pval, genes_bipolar, genes_sfari_non_syndromic, genes_schema_or,private=FALSE) {

    if(private){
        df <- df %>% filter(count == 1)  # Only filter if private=TRUE
    } 
    # Assume expanded_df is already created as per your code
   expanded_df <- df %>%
    separate_rows(SAMPLES, sep = ",") %>%
    mutate(sample_id = trimws(SAMPLES)) %>%
    left_join(
        manifest_correct %>% distinct(sample_id, Status), 
        by = "sample_id"
    ) %>%
    filter(!sample_id %in% UHR_NA_SAMPLE_IDS$V1) %>%
    filter(!sample_id %in% Low_Quality_SAMPLES_AND_UHR_NA$Sequencing_number) %>%
    distinct(SYMBOL, CHROM, REF, ALT, POS, sample_id, .keep_all = TRUE)


    expanded_df_pure <- df %>% 
        separate_rows(SAMPLES, sep = ",") %>%  # Split SAMPLES only
        mutate(
            sample_id = trimws(SAMPLES),
            Status = sample_label_3  # Use pre-computed labels
        ) %>% filter(Status %in% c("FEP", "Non_Converter", "Converter")) %>% 
        filter(!sample_id %in% UHR_NA_SAMPLE_IDS$V1) %>%  # Remove UHR_NA
        filter(!sample_id %in% Low_Quality_SAMPLES_AND_UHR_NA$Sequencing_number) %>%
        distinct(SYMBOL, CHROM,REF,ALT, POS, sample_id, .keep_all = TRUE)  # Unique variant/sample
        
   
    # Filter based on the genes in the brain, brain filter_ntpm, and additional gene sets
    expanded_df <- expanded_df %>%
        mutate(brain = ifelse(SYMBOL %in% brain_gene_consensus_filtered_consensus_no_pitular$Gene.name, 1, 0),
               brain_filter_ntpm = ifelse(SYMBOL %in% brain_gene_consensus_ntm_consensus_no_pitular$Gene.name, 1, 0),
               schema_pval = ifelse(SYMBOL %in% genes_schema_pval, 1, 0),
               bipolar = ifelse(SYMBOL %in% genes_bipolar, 1, 0),
               sfari_non_syndromic = ifelse(SYMBOL %in% genes_sfari_non_syndromic, 1, 0),
               schema_or = ifelse(SYMBOL %in% genes_schema_or, 1, 0))  

    expanded_df_pure <- expanded_df_pure %>%
        mutate(brain = ifelse(SYMBOL %in% brain_gene_consensus_filtered_consensus_no_pitular$Gene.name, 1, 0),
               brain_filter_ntpm = ifelse(SYMBOL %in% brain_gene_consensus_ntm_consensus_no_pitular$Gene.name, 1, 0),
               schema_pval = ifelse(SYMBOL %in% genes_schema_pval, 1, 0),
               bipolar = ifelse(SYMBOL %in% genes_bipolar, 1, 0),
               sfari_non_syndromic = ifelse(SYMBOL %in% genes_sfari_non_syndromic, 1, 0),
               schema_or = ifelse(SYMBOL %in% genes_schema_or, 1, 0))     

    # Count the number of variants per individual for each category
    result_df <- expanded_df %>%
        group_by(sample_id, Status) %>%
        summarise(
            Number_of_Variants = n(),
            Number_of_Brain_Variants = sum(brain),
            Number_of_Brain_Filter_NTPM_Variants = sum(brain_filter_ntpm),
            Number_of_Schema_Pval_Variants = sum(schema_pval),
            Number_of_Bipolar_Variants = sum(bipolar),
            Number_of_Sfari_Non_Syndromic_Variants = sum(sfari_non_syndromic),
            Number_of_Schema_Or_Variants = sum(schema_or),
            .groups = "drop"
        )

    result_df_pure <- expanded_df_pure %>%  
    group_by(sample_id, Status) %>%
        summarise(
            Number_of_Variants_pure = n(),
            Number_of_Brain_Variants_pure = sum(brain),
            Number_of_Brain_Filter_NTPM_Variants_pure = sum(brain_filter_ntpm),
            Number_of_Schema_Pval_Variants_pure = sum(schema_pval),
            Number_of_Bipolar_Variants_pure = sum(bipolar),
            Number_of_Sfari_Non_Syndromic_Variants_pure = sum(sfari_non_syndromic),
            Number_of_Schema_Or_Variants_pure = sum(schema_or),
            .groups = "drop"
        )   
    merged_results <- result_df %>%
    full_join(result_df_pure, by = "sample_id", suffix = c("_status", "_group")) %>%
    mutate(
        # Create a single Status column, only if both are equal or one is NA
        Status = case_when(
            !is.na(Status_status) & !is.na(Status_group) & Status_status == Status_group ~ Status_status,
            is.na(Status_status) ~ Status_group,
            is.na(Status_group) ~ Status_status,
            TRUE ~ NA_character_  # mismatch case
        )
    ) %>%
    select(-Status_status, -Status_group) %>%
    distinct()


    return(merged_results)
}



# Create separate data frames for each call to process_df
     # Process each dataframe and create separate result tables
    result_df_Missense <- process_df(Missense, "Missense", brain_gene_consensus_filtered_consensus_no_pitular, brain_gene_consensus_ntm_consensus_no_pitular, genes_schema_pval, genes_bipolar, genes_sfari_non_syndromic, genes_schema_or,private=FALSE)
    result_df_Missense_canonical <- process_df(Missense_canonical, "Missense_canonical", brain_gene_consensus_filtered_consensus_no_pitular, brain_gene_consensus_ntm_consensus_no_pitular, genes_schema_pval, genes_bipolar, genes_sfari_non_syndromic, genes_schema_or,private=FALSE)
    result_df_MPC_only <- process_df(MPC_only, "MPC_only", brain_gene_consensus_filtered_consensus_no_pitular, brain_gene_consensus_ntm_consensus_no_pitular, genes_schema_pval, genes_bipolar, genes_sfari_non_syndromic, genes_schema_or,private=FALSE)
    result_df_PTV_canonic <- process_df(PTV_canonic, "PTV_canonic", brain_gene_consensus_filtered_consensus_no_pitular, brain_gene_consensus_ntm_consensus_no_pitular, genes_schema_pval, genes_bipolar, genes_sfari_non_syndromic, genes_schema_or,private=FALSE)
    result_df_PTV <- process_df(PTV, "PTV", brain_gene_consensus_filtered_consensus_no_pitular, brain_gene_consensus_ntm_consensus_no_pitular, genes_schema_pval, genes_bipolar, genes_sfari_non_syndromic, genes_schema_or,private=FALSE)
  

    result_df_Missense <- process_df(Missense, "Missense", brain_gene_consensus_filtered_consensus_no_pitular, brain_gene_consensus_ntm_consensus_no_pitular, genes_schema_pval, genes_bipolar, genes_sfari_non_syndromic, genes_schema_or,private=TRUE)
    result_df_Missense_canonical <- process_df(Missense_canonical, "Missense_canonical", brain_gene_consensus_filtered_consensus_no_pitular, brain_gene_consensus_ntm_consensus_no_pitular, genes_schema_pval, genes_bipolar, genes_sfari_non_syndromic, genes_schema_or,private=TRUE)
    result_df_MPC_only <- process_df(MPC_only, "MPC_only", brain_gene_consensus_filtered_consensus_no_pitular, brain_gene_consensus_ntm_consensus_no_pitular, genes_schema_pval, genes_bipolar, genes_sfari_non_syndromic, genes_schema_or,private=TRUE)
    result_df_PTV_canonic <- process_df(PTV_canonic, "PTV_canonic", brain_gene_consensus_filtered_consensus_no_pitular, brain_gene_consensus_ntm_consensus_no_pitular, genes_schema_pval, genes_bipolar, genes_sfari_non_syndromic, genes_schema_or,private=TRUE)
    result_df_PTV <- process_df(PTV, "PTV", brain_gene_consensus_filtered_consensus_no_pitular, brain_gene_consensus_ntm_consensus_no_pitular, genes_schema_pval, genes_bipolar, genes_sfari_non_syndromic, genes_schema_or,private=TRUE)




manifest_correct <- manifest_correct %>% filter(!grepl("UHR-NA", Status))


# ADD ALL SAMPLES INDEPENDENTLY OF PRESENCE/ABSENCE OF VARIANTS 
# Ensure consistency in column names and perform left join for all dataframes
complete_df_Missense <- manifest_correct %>%
    left_join(result_df_Missense, by = c("sample_id", "Status")) %>%
    mutate(across(starts_with("Number_of"), ~ ifelse(is.na(.), 0, .)))

complete_df_Missense_canonical <- manifest_correct %>%
    left_join(result_df_Missense_canonical, by = c("sample_id", "Status")) %>%
    mutate(across(starts_with("Number_of"), ~ ifelse(is.na(.), 0, .)))    

complete_df_MPC_only <- manifest_correct %>%
    left_join(result_df_MPC_only, by = c("sample_id", "Status")) %>%
    mutate(across(starts_with("Number_of"), ~ ifelse(is.na(.), 0, .)))

complete_df_PTV_canonic <- manifest_correct %>%
    left_join(result_df_PTV_canonic, by = c("sample_id", "Status")) %>%
    mutate(across(starts_with("Number_of"), ~ ifelse(is.na(.), 0, .)))

complete_df_PTV <- manifest_correct %>%
    left_join(result_df_PTV, by = c("sample_id", "Status")) %>%
    mutate(across(starts_with("Number_of"), ~ ifelse(is.na(.), 0, .)))


complete_df_Missense$Status <- relevel(factor(complete_df_Missense$Status), ref = "Non_Converter")
complete_df_Missense_canonical$Status <- relevel(factor(complete_df_Missense_canonical$Status), ref = "Non_Converter")
complete_df_MPC_only$Status <- relevel(factor(complete_df_MPC_only$Status), ref = "Non_Converter")
complete_df_PTV_canonic$Status <- relevel(factor(complete_df_PTV_canonic$Status), ref = "Non_Converter")
complete_df_PTV$Status <- relevel(factor(complete_df_PTV$Status), ref = "Non_Converter")


multinom_model <- function(data) {
  results <- purrr::map_dfr(
    .x = names(data)[grepl("Number_of_", names(data))],
    .f = function(col_name) {
      formula <- as.formula(paste("Status ~", col_name))
      model <- nnet::multinom(formula, data = data, trace = FALSE)

      # Tidy output (returns one row per coefficient per outcome level vs. baseline)
      tidy_model <- broom::tidy(model, conf.int = TRUE, exponentiate = TRUE) %>%
        mutate(predictor = col_name)
      
      return(tidy_model)
    }
  )

  results_clean <- results %>%
    rename(
      odds_ratio = estimate,
      lower_OR = conf.low,
      upper_OR = conf.high,
      p_val = p.value
    ) %>%
    mutate(
      log_odds = log(odds_ratio),
      adj_p_val = ifelse(term == "(Intercept)", NA, p.adjust(p_val, method = "fdr"))
    ) %>%
    filter(term != "(Intercept)") %>%
    select(!predictor)

  return(results_clean)
}

Missense_multinom_results <- multinom_model(complete_df_Missense)
Missense_canonical_multinom_results <- multinom_model(complete_df_Missense_canonical)
MPC_only_multinom_results <- multinom_model(complete_df_MPC_only)       
PTV_canonic_multinom_results <- multinom_model(complete_df_PTV_canonic)
PTV_multinom_results <- multinom_model(complete_df_PTV)

# Combine all multinom results into one dataframe
all_multinom_results <- bind_rows(
  Missense_multinom_results %>% mutate(variant_type = "Missense"),
  Missense_canonical_multinom_results %>% mutate(variant_type = "Missense_canonical"),
  MPC_only_multinom_results %>% mutate(variant_type = "MPC_only"),
  PTV_canonic_multinom_results %>% mutate(variant_type = "PTV_canonic"),
  PTV_multinom_results %>% mutate(variant_type = "PTV")
)


#zeroinflation poisson regression

zero_poisson <- function(data) {
  results <- purrr::map_dfr(
    .x = names(data)[grepl("Number_of_", names(data))],
    .f = function(col_name) {
      formula <- as.formula(paste(col_name, "~ Status | 1"))
      model <- zeroinfl(formula, data = data, dist = "poisson")
      
      # Extract count model coefficients
      coefs <- summary(model)$coefficients$count
      terms <- rownames(coefs)
      
      df <- as.data.frame(coefs)
      df$term <- terms
      df$predictor <- col_name
      
      df
    }
  )
  
  # Clean and calculate odds ratio + CIs
  results_clean <- results %>%
    filter(term != "(Intercept)") %>%
    mutate(
      log_odds = Estimate,
      SE = `Std. Error`,
      p_val = `Pr(>|z|)`,
      odds_ratio = exp(log_odds),
      lower_OR = exp(log_odds - 1.96 * SE),
      upper_OR = exp(log_odds + 1.96 * SE),
      adj_p_val = p.adjust(p_val, method = "fdr")
    ) %>%
    select(predictor, term, log_odds, SE, odds_ratio, lower_OR, upper_OR, p_val, adj_p_val)
  
  return(results_clean)
}
Missense_poisson_results <- zero_poisson(complete_df_Missense)
Missense_canonical_poisson_results <- zero_poisson(complete_df_Missense_canonical)
MPC_only_poisson_results <- zero_poisson(complete_df_MPC_only)
PTV_canonic_poisson_results <- zero_poisson(complete_df_PTV_canonic)
PTV_poisson_results <- zero_poisson(complete_df_PTV)

#combine all results into a single dataframe of poisson results
all_poisson_results <- bind_rows(
  Missense_poisson_results %>% mutate(variant_type = "Missense"),
  Missense_canonical_poisson_results %>% mutate(variant_type = "Missense_canonical"),
  MPC_only_poisson_results %>% mutate(variant_type = "MPC_only"),
  PTV_canonic_poisson_results %>% mutate(variant_type = "PTV_canonic"),
  PTV_poisson_results %>% mutate(variant_type = "PTV")
)



volcano_poisson <- all_poisson_results %>%
  mutate(
    log2_OR = log2(odds_ratio),
    neg_log10_p = -log10(adj_p_val),
    Label = ifelse(adj_p_val < 0.05, term, NA)
  ) %>%
  ggplot(aes(x = log2_OR, y = neg_log10_p)) +
  geom_point(aes(shape=variant_type ,color = predictor), alpha = 0.7, size = 3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_label_repel(aes(label = Label), size = 2.5, box.padding = 0.5, max.overlaps = 20) +
  labs(
    title = "Volcano Plot: Zero-Inflated Poisson Results",
    x = "log2(Odds Ratio)",
    y = "-log10(FDR-adjusted p-value)",
    color = "Variant Type"
  ) +
  theme_bw()

volcano_multinom <- all_multinom_results %>%
    mutate(
        neg_log10_p = -log10(adj_p_val),
        Label = ifelse(adj_p_val < 0.05, paste(term, y.level, sep = " - "), NA)
    ) %>%
    ggplot(aes(x = odds_ratio, y = neg_log10_p)) +
    geom_point(aes(shape=y.level ,color = variant_type), alpha = 0.7, size = 3) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_label_repel(aes(label = Label), size = 2.5, box.padding = 0.5, max.overlaps = 20) +
    labs(
        title = "Volcano Plot: Multinomial Regression Results",
        x = "Log-Odds",
        y = "-log10(FDR-adjusted p-value)",
        color = "Variant Type"
    ) +
    theme_bw()


setwd("/home/rachele/SNVs/results_pasteur/group3/private")
# Save the plots
ggsave("volcano_poisson_plot.png", plot = volcano_poisson, width = 10, height = 6)
ggsave("volcano_multinom_plot.png", plot = volcano_multinom, width = 10, height = 6)

# Save the results to CSV files
write.csv(all_poisson_results, "all_poisson_results.csv", row.names = FALSE)
write.csv(all_multinom_results, "all_multinom_results.csv", row.names = FALSE)
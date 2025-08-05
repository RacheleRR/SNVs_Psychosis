

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
library(ggplot2)



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

# ADD LABELS 
# Function to get the labels (case or control) for the outlier values
get_outlier_labels <- function(outlier_value, manifest_df) {
  outlier_values <- strsplit(outlier_value, ",")[[1]]
  labels <- sapply(outlier_values, function(val) {
    if (val %in% manifest_df$Sequencing_number) {
      return(manifest_df$new_column[manifest_df$Sequencing_number == val])
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
        # Find the label in the manifest dataframe (either "case" or "control")
        if (val %in% manifest_df$Sequencing_number) {
            return(manifest_df$new_column[manifest_df$Sequencing_number == val])
        } else {
            return(NA) # Return NA if no match is found
        }
    })

    # If there's more than one label, return the unique ones
    unique_labels <- unique(labels)
    
    # If there's both "case" and "control", return "mixed"; otherwise return the label
    if (length(unique_labels) > 1) {
        return("mixed")
    } else if (length(unique_labels) == 1) {
        return(unique_labels)
    } else {
        return(NA) # If no label found
    }
}

# DO the UHR NA CONTROL

    # Function to clean UHR_NA from a dataframe
    clean_uh_na <- function(df) {
    # Step 1: Remove rows where sample_label contains "UHR_NA"
    df <- df %>%
        filter(!grepl("UHR_NA", sample_label))
    
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
    
    # Determine the group based on the labels
    if (all(labels == "Case")) {
        return("Case")
    } else if (all(labels == "Control")) {
        return("Control")
    } else if (any(labels == "Case") && any(labels == "Control")) {
        return("mixed")
    } else {
        return(NA)  # If no valid labels are found
    }
    }

# Function to clean NA samples
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

rm('Rubiu-GRCh38.deepvariant.splitted.norm.vep.merged.filtered_gnomad_mpc2.filtered_rarity_stringent_AF.missense_mpc_AM.vcf.gz')
rm('Rubiu-GRCh38.deepvariant.splitted.norm.vep.merged.filtered_gnomad_mpc2.filtered_rarity_stringent_AF.missense_mpc_AM_CANONIC.vcf.gz')
rm('Rubiu-GRCh38.deepvariant.splitted.norm.vep.merged.filtered_gnomad_mpc2.filtered_rarity_stringent_AF.missense_mpc.vcf.gz')
rm(list = c("SCHEMA", "SCHEMA_OR", "SCHEMA_pVal", "SFARI", "SFARI_non_syndromic_lower_3", "SFARI_syndromic", "GWAS_120", "BipEx_Bipolar", "BipEx_Bipolar_p_val_Missense", "BipEx_Bipolar_p_val_PTV", "convert_col"))


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
    Missense <- clean_na_samples(Missense)
    Missense_canonical <- clean_na_samples(Missense_canonical)
    MPC_only <- clean_na_samples(MPC_only)
    PTV_canonic <- clean_na_samples(PTV_canonic)
    PTV <- clean_na_samples(PTV)

    #apply clean uhna and derive group label 
    Missense <- clean_uh_na(Missense)
    Missense_canonical <- clean_uh_na(Missense_canonical)
    MPC_only <- clean_uh_na(MPC_only)
    PTV_canonic <- clean_uh_na(PTV_canonic)
    PTV <- clean_uh_na(PTV)   

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

    UHR_NA_SAMPLE_IDS <- read.csv("~/UHR_NA_SAMPLE_IDS.csv", sep="")
    #example S36827

colnames(manifest_correct) <- c("Status", "sample_id", "group")


# stuff
process_df <- function(df, name, brain_gene_consensus_filtered_consensus_no_pitular, brain_gene_consensus_ntm_consensus_no_pitular, genes_schema_pval, genes_bipolar, genes_sfari_non_syndromic, genes_schema_or,private=FALSE) {

    if(private){
        df <- df %>% filter(count == 1)  # Only filter if private=TRUE
    }
    
    # Assume expanded_df is already created as per your code
   expanded_df <- df %>%
    separate_rows(SAMPLES, sep = ",") %>%
    mutate(sample_id = trimws(SAMPLES)) %>%
    left_join(
        manifest_correct %>% distinct(sample_id, group), 
        by = "sample_id"
    ) %>%
    filter(!sample_id %in% UHR_NA_SAMPLE_IDS$V1) %>%  # Remove UHR_NA
    filter(!sample_id %in% Low_Quality_SAMPLES_AND_UHR_NA$Sequencing_number) %>%
    distinct(SYMBOL, CHROM, REF, ALT, POS, sample_id, .keep_all = TRUE)


    expanded_df_pure <- df %>% 
        separate_rows(SAMPLES, sep = ",") %>%  # Split SAMPLES only
        mutate(
            sample_id = trimws(SAMPLES),
            group = sample_label_3  # Use pre-computed labels
        ) %>% filter(group %in% c("Case", "Control")) %>% 
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
        group_by(sample_id, group) %>%
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
    group_by(sample_id, group) %>%
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
    full_join( result_df_pure,
    by = "sample_id"  ) %>% mutate(
    group = coalesce(group.x, group.y)  # Combine group.x and group.y into "group"
    ) %>%select(-group.x, -group.y) %>%       # Remove the original group columns
    distinct()                            # Remove duplicates

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




manifest_correct <- manifest_correct %>% filter(!grepl("UHR_NA", group))


# ADD ALL SAMPLES INDEPENDENTLY OF PRESENCE/ABSENCE OF VARIANTS 
# Ensure consistency in column names and perform left join for all dataframes
complete_df_Missense <- manifest_correct %>%
    left_join(result_df_Missense, by = c("sample_id", "group")) %>%
    mutate(across(starts_with("Number_of"), ~ ifelse(is.na(.), 0, .)))

complete_df_Missense_canonical <- manifest_correct %>%
    left_join(result_df_Missense_canonical, by = c("sample_id", "group")) %>%
    mutate(across(starts_with("Number_of"), ~ ifelse(is.na(.), 0, .)))    

complete_df_MPC_only <- manifest_correct %>%
    left_join(result_df_MPC_only, by = c("sample_id", "group")) %>%
    mutate(across(starts_with("Number_of"), ~ ifelse(is.na(.), 0, .)))

complete_df_PTV_canonic <- manifest_correct %>%
    left_join(result_df_PTV_canonic, by = c("sample_id", "group")) %>%
    mutate(across(starts_with("Number_of"), ~ ifelse(is.na(.), 0, .)))

complete_df_PTV <- manifest_correct %>%
    left_join(result_df_PTV, by = c("sample_id", "group")) %>%
    mutate(across(starts_with("Number_of"), ~ ifelse(is.na(.), 0, .)))




# statistical test 
# without normalize

    analyze_metrics <- function(df) {
    # 1. Test each variant count metric
    columns_to_test <- grep("^Number_of", names(df), value = TRUE)
    
    results <- map_dfr(columns_to_test, function(col) {
        # Initialize result row
        result_row <- tibble(
        metric = col,
        case_zero = FALSE,
        control_zero = FALSE,
        shapiro_case_p = NA_real_,
        shapiro_control_p = NA_real_,
        primary_test = NA_character_,
        primary_p = NA_real_,
        secondary_test = NA_character_,
        secondary_p = NA_real_,
        note = NA_character_
        )
        
        # Split data
        case_data <- df[[col]][df$group == "Case"]
        control_data <- df[[col]][df$group == "Control"]
        
        # Check for all-zero cases
        result_row$case_zero <- all(case_data == 0)
        result_row$control_zero <- all(control_data == 0)
        
        if (result_row$case_zero || result_row$control_zero) {
        result_row$note <- case_when(
            result_row$case_zero && result_row$control_zero ~ "All zeros in both groups",
            result_row$case_zero ~ "All zeros in Case group",
            TRUE ~ "All zeros in Control group"
        )
        return(result_row)
        }
        
        # Handle constant values for Shapiro
        safe_shapiro <- function(x) {
        if (length(unique(x)) < 3) return(list(p.value = NA))
        tryCatch(shapiro.test(x), error = function(e) list(p.value = NA))
        }
        
        # Normality checks
        shapiro_case <- safe_shapiro(case_data)
        shapiro_control <- safe_shapiro(control_data)
        result_row$shapiro_case_p <- shapiro_case$p.value
        result_row$shapiro_control_p <- shapiro_control$p.value
        
        # Determine test strategy
        if (!any(is.na(c(shapiro_case$p.value, shapiro_control$p.value))) &&
        shapiro_case$p.value > 0.05 && 
        shapiro_control$p.value > 0.05) {
        primary <- "t-test"
        secondary <- "wilcox"
        } else {
        primary <- "wilcox"
        secondary <- "t-test"
        }
        
        # Perform tests
        tryCatch({
        t_res <- t.test(df[[col]] ~ df$group) %>% tidy()
        w_res <- wilcox.test(df[[col]] ~ df$group) %>% tidy()
        
        result_row$primary_test <- primary
        result_row$primary_p <- if(primary == "t-test") t_res$p.value else w_res$p.value
        result_row$secondary_test <- secondary
        result_row$secondary_p <- if(secondary == "t-test") t_res$p.value else w_res$p.value
        }, error = function(e) {
        result_row$note <- paste("Test error:", e$message)
        })



        result_row
        })
           # Add FDR-adjusted p-values
         results %>%
         mutate(
            primary_p_adj = p.adjust(primary_p, method = "fdr"),
            secondary_p_adj = p.adjust(secondary_p, method = "fdr")
           ) %>%
          select(
          metric, 
           primary_test, primary_p, primary_p_adj,
             secondary_test, secondary_p, secondary_p_adj,
             everything()
    )
    }
        
    # apply 
    results_pred_Missense <- analyze_metrics(complete_df_Missense) %>%  select(-case_zero, -control_zero)
    results_pred_Missense_canonical <- analyze_metrics(complete_df_Missense_canonical) %>%  select(-case_zero, -control_zero)
    results_pred_MPC_only <- analyze_metrics(complete_df_MPC_only) %>%  select(-case_zero, -control_zero)
    results_pred_PTV_canonic <- analyze_metrics(complete_df_PTV_canonic) %>%  select(-case_zero, -control_zero)
    results_pred_PTV <- analyze_metrics(complete_df_PTV) %>%  select(-case_zero, -control_zero)

#
# VISUALIZATION
setwd("/home/rachele/SNVs/results_pasteur/group2/private")
# Save all the plots indiscriminate from significance  
        plot <- function(df, results_df, title_suffix) {
                # 2. Prepare annotation data from test results
                annotation_data <- results_df %>%
                    filter(str_detect(metric, "^Number_of")) %>%
                    mutate(
                        metric_clean = gsub("Number_of_", "", metric),
                        label = paste0(primary_test, "\nFDR p = ", format.pval(primary_p_adj, digits = 2)),
                        group = NA  # Add dummy group variable
                    )
                
                # 3. Reshape for plotting
                plot_data <- df %>%
                    select(sample_id, group, starts_with("Number_of")) %>%
                    pivot_longer(
                        cols = -c(sample_id, group),
                        names_to = "metric",
                        values_to = "value"
                    ) %>%
                    mutate(
                        metric_clean = gsub("Number_of_", "", metric)
                    )
                
                # 4. Create the plot
                ggplot(plot_data, aes(x = group, y = value, fill = group)) +
                    geom_boxplot(outlier.shape = NA) +
                    geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
                    facet_wrap(~ metric_clean, scales = "free_y", ncol = 3) +
                    geom_text(
                        data = annotation_data,
                        aes(x = 1.5, y = Inf, label = label),
                        inherit.aes = FALSE,  # Don't inherit aesthetics from main plot
                        vjust = 1.5, size = 3, color = "black"
                    ) +
                    labs(
                        title = paste("Case vs. Control:", title_suffix),
                        x = "Group",
                        y = "Variants per individual (raw)",
                        fill = "Group"
                    ) +
                    theme_bw() +
                    theme(
                        strip.background = element_blank(),
                        strip.text = element_text(face = "bold")
                    )
            }

            plot_Missense <- plot(
            complete_df_Missense, 
            results_pred_Missense,
            "Missense"
            )

            plot_Missense_canonical <- plot(
            complete_df_Missense_canonical, 
            results_pred_Missense_canonical,
            "Missense canonical"
            )

            plot_MPC_only <- plot(
            complete_df_MPC_only, 
            results_pred_MPC_only,
            "MPC only"
            )

            plot_PTV_canonic <- plot(
            complete_df_PTV_canonic, 
            results_pred_PTV_canonic,
            "PTV canonical"
            )

            plot_PTV <- plot(
            complete_df_PTV, 
            results_pred_PTV,
            "PTV"
            )

            # Display one plot as example
            plot_Missense_canonical

            # To save all plots:
            plots <- list(
            plot_Missense, plot_Missense_canonical,
            plot_MPC_only, 
            plot_PTV_canonic, plot_PTV
            )


            walk2(plots, c("Missense", "Missense_canonical", "MPC_only", "PTV_canonic", "PTV"), 
                ~ ggsave(paste0("wilcoxon_variant_comparison_", .y, ".png"), .x, width = 10, height = 8))


        #save tabels 
        combined_results_non_norm <- bind_rows(
        "Missense" = results_pred_Missense,
        "Missense_canonical" = results_pred_Missense_canonical,
        "MPC_only" = results_pred_MPC_only,
        "PTV_canonic" = results_pred_PTV_canonic,
        "PTV" = results_pred_PTV,
        .id = "analysis_type"  # Adds a column to identify the source
        )

# ONLY  show and save significant results
    plot_significant <- function(df, results_df, title_suffix) {
    # Filter for significant metrics (Wilcoxon FDR p < 0.05)
    sig_metrics <- results_df %>%
        filter(
        str_detect(metric, "^Number_of"),
        primary_p_adj < 0.05, 
        primary_test == "wilcox"  # Focus on Wilcoxon results
        ) %>%
        pull(metric)
    
    # Return NULL if no significant results
    if (length(sig_metrics) == 0) {
        message("No significant results for: ", title_suffix)
        return(NULL)
    }
    
    # Prepare annotation data (only for significant metrics)
    annotation_data <- results_df %>%
        filter(metric %in% sig_metrics) %>%
        mutate(
        metric_clean = gsub("Number_of_", "", metric),
        label = paste0("Wilcoxon\nFDR p = ", format.pval(primary_p_adj, digits = 2)),
        group = NA  # Dummy group for positioning
        )
    
    # Reshape data for plotting
    plot_data <- df %>%
        select(sample_id, group, all_of(sig_metrics)) %>%
        pivot_longer(
        cols = -c(sample_id, group),
        names_to = "metric",
        values_to = "value"
        ) %>%
        mutate(
        metric_clean = gsub("Number_of_", "", metric)
        )
    
    # Create the plot
    ggplot(plot_data, aes(x = group, y = value, fill = group)) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
        facet_wrap(~ metric_clean, scales = "free_y") +
        geom_text(
        data = annotation_data,
        aes(x = 1.5, y = Inf, label = label),
        vjust = 1.5, size = 3, color = "black"
        ) +
        labs(
        title = paste("Significant results:", title_suffix),
        x = "Group",
        y = "Variants per individual",
        fill = "Group"
        ) +
        theme_bw() +
        theme(
        strip.background = element_blank(),
        strip.text = element_text(face = "bold")
        )
    }

    # Generate and save only significant plots
    plot_inputs <- list(
    "Missense" = list(df = complete_df_Missense, results = results_pred_Missense, title = "Missense"),
    "Missense_canonical" = list(df = complete_df_Missense_canonical, results = results_pred_Missense_canonical, title = "Missense (canonical)"),
    "MPC_only" = list(df = complete_df_MPC_only, results = results_pred_MPC_only, title = "MPC only"),
    "PTV_canonic" = list(df = complete_df_PTV_canonic, results = results_pred_PTV_canonic, title = "PTV (canonical)"),
    "PTV" = list(df = complete_df_PTV, results = results_pred_PTV, title = "PTV")
    )

    # Apply function and keep names
    sig_plots_named <- imap(plot_inputs, function(input, name) {
    p <- plot_significant(input$df, input$results, input$title)
    if (!is.null(p)) list(plot = p, name = name) else NULL
    }) %>% compact()

    # Save significant plots
    walk(sig_plots_named, function(entry) {
    ggsave(
        paste0("wilcoxon_SIGNIFICANT_variant_comparison_", entry$name, ".png"),
        entry$plot,
        width = 10,
        height = 6
    )
    })

# Save the combined results to a CSV file
write.csv(combined_results_non_norm, "wilcoxon_combined_results_non_norm.csv", row.names = FALSE)



#posssibilty to do more 


# Load the combined results (normalized or non-normalized)
combined_results <- combined_results_non_norm  # or combined_results_non_norm

# Filter for the metrics we care about (Number_of_Variants and Number_of_Variants_pure)
p_value_data <- combined_results %>%
    filter(metric %in% c("Number_of_Variants", "Number_of_Variants_pure")) %>%
    select(analysis_type, metric, primary_p_adj) %>%
    mutate(
        metric_clean = case_when(
            metric == "Number_of_Variants" ~ "All Samples",
            metric == "Number_of_Variants_pure" ~ "Pure Samples"
        )
    )


    # List of all complete dataframes
complete_dfs <- list(
  "Missense" = complete_df_Missense,
  "Missense_canonical" = complete_df_Missense_canonical,
  "MPC_only" = complete_df_MPC_only,
  "PTV_canonic" = complete_df_PTV_canonic,
  "PTV" = complete_df_PTV
)

# Extract and combine the relevant columns
combined_variants <- map_dfr(complete_dfs, ~ {
  .x %>%
    select(sample_id, group, Number_of_Variants, Number_of_Variants_pure) %>%
    pivot_longer(
      cols = c(Number_of_Variants, Number_of_Variants_pure),
      names_to = "variant_type",
      values_to = "count"
    ) %>%
    mutate(
      variant_type = case_when(
        variant_type == "Number_of_Variants" ~ "All Samples",
        variant_type == "Number_of_Variants_pure" ~ "Pure Samples"
      )
    )
}, .id = "analysis_type")





# Merge with the plotting data
plot_data <- combined_variants %>%
    left_join(p_value_data, by = c("analysis_type", "variant_type" = "metric_clean"))


# 3. Create the plot
variant_comparison_plot <- ggplot(plot_data, aes(x = group, y = count, fill = group)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
    facet_grid(variant_type ~ analysis_type, scales = "free_y") +
    
    # Add p-value annotations (using p_value_data)
    geom_text(
        data = plot_data %>%
            group_by(analysis_type, variant_type) %>%
            slice(1),  # Take one row per facet to annotate
        aes(
            x = 1.5, 
            y = Inf, 
            label = paste("FDR p =", format.pval(primary_p_adj, digits = 2))
        ),
        inherit.aes = FALSE,
        vjust = 1.5, 
        size = 3, 
        color = "black"
    )+

    
    labs(
        title = "Case vs. Control: Variant Counts (All vs. Pure Samples)",
        x = "Group",
        y = "Number of Variants (normalized per individual)",
        fill = "Group"
    ) +
    theme_bw() +
    theme(
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1)
    )

# Display the plot
print(variant_comparison_plot)

# Save the plot
ggsave("wilcoxon_variant_comparison_plot.png", variant_comparison_plot, width = 12, height = 8)


# ADDITIONAL REGRESSION AND POISSSON 
complete_df_Missense_canonical <- complete_df_Missense_canonical %>% mutate(num_group =  ifelse(group %in% "Case", 1, 0))
lt <- glm(num_group ~ Number_of_Variants , data = complete_df_Missense_canonical ,family = binomial)
summary(lt)
lt <- glm(num_group ~ Number_of_Variants_pure , data = complete_df_Missense_canonical ,family = binomial)
summary(lt)

lt <- glm(num_group ~ Number_of_Brain_Variants , data = complete_df_Missense_canonical ,family = binomial)
summary(lt)
lt <- glm(num_group ~ Number_of_Brain_Filter_NTPM_Variants , data = complete_df_Missense_canonical ,family = binomial)
summary(lt)
lt <- glm(num_group ~ Number_of_Schema_Pval_Variants , data = complete_df_Missense_canonical ,family = binomial)
summary(lt)
lt <- glm(num_group ~ Number_of_Bipolar_Variants , data = complete_df_Missense_canonical ,family = binomial)
summary(lt)

> summary(zeroinfl(Number_of_Variants ~ group | 1, 
+                  data = complete_df_Missense3, 
+                  dist = "poisson"))

> summary(zeroinfl(Number_of_Variants ~ num_group | 1, 
+                  data = complete_df_Missense_canonical3, 
+                  dist = "poisson"))





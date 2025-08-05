# wilcoxon_sNVs_PASTEUR_3GROUPS





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
rm('Rubiu-GRCh38.deepvariant.splitted.norm.vep.merged.filtered_gnomad_mpc2.filtered_rarity_stringent_AF.missense_mpc_AM.vcf.gz')
rm('Rubiu-GRCh38.deepvariant.splitted.norm.vep.merged.filtered_gnomad_mpc2.filtered_rarity_stringent_AF.missense_mpc_AM_CANONIC.vcf.gz')
rm('Rubiu-GRCh38.deepvariant.splitted.norm.vep.merged.filtered_gnomad_mpc2.filtered_rarity_stringent_AF.missense_mpc.vcf.gz')
rm(list = c("SCHEMA", "SCHEMA_OR", "SCHEMA_pVal", "SFARI", "SFARI_non_syndromic_lower_3", "SFARI_syndromic", "GWAS_120", "BipEx_Bipolar", "BipEx_Bipolar_p_val_Missense", "BipEx_Bipolar_p_val_PTV", "convert_col"))


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
colnames(manifest_correct) <- c("Status", "sample_id", "group")




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
            group = sample_label_3  # Use pre-computed labels
        ) %>% filter(group %in% c("FEP", "Non_Converter", "Converter")) %>% 
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
     full_join(
    result_df_pure,
    by = "sample_id",
    suffix = c("_status", "_group")
  ) %>%distinct()  # Remove duplicate combinations


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


#kruskal 
analyze_metrics_kruskal <- function(df, normalize = TRUE, df_name = NULL) {
  # Define group sizes
  total_individuals <- c(
    FEP = 255,
    Converter = 47,
    Non_Converter = 75
  )
  
  # 1. Filter data
  df_filtered <- df %>% 
    filter(Status %in% names(total_individuals))
  
  # 2. Conditionally normalize data
  if(normalize) {
    df_processed <- df_filtered %>%
      mutate(across(starts_with("Number_of"), 
                   ~ ./total_individuals[Status],
                   .names = "{.col}_per_indiv"))
    columns_to_test <- grep("_per_indiv$", names(df_processed), value = TRUE)
  } else {
    df_processed <- df_filtered
    columns_to_test <- grep("^Number_of", names(df_processed), value = TRUE)
  }
  
  # 3. Test each metric
  results <- map_dfr(columns_to_test, function(col) {
    # Initialize result row
    result_row <- tibble(
      dataset = ifelse(!is.null(df_name), df_name, deparse(substitute(df))),
      metric = col,
      normalization = normalize,
      kruskal_statistic = NA_real_,
      kruskal_p = NA_real_,
      kruskal_p_adj = NA_real_,
      comparison = NA_character_,
      dunn_z = NA_real_,
      dunn_p = NA_real_,
      dunn_p_adj = NA_real_,
      note = NA_character_
    )
    
    # Extract data
    test_data <- df_processed[c("Status", col)]
    colnames(test_data) <- c("group", "value")
    
    # Check for all-zero cases in any group
    zero_groups <- test_data %>%
      group_by(group) %>%
      summarise(all_zero = all(value == 0)) %>%
      filter(all_zero) %>%
      pull(group)
    
    if(length(zero_groups) > 0) {
      result_row$note <- paste("All zeros in groups:", paste(zero_groups, collapse = ", "))
      return(result_row)
    }
    
    # Perform Kruskal-Wallis test
    kruskal_res <- kruskal.test(value ~ group, data = test_data)
    result_row$kruskal_statistic <- kruskal_res$statistic
    result_row$kruskal_p <- kruskal_res$p.value
    result_row$kruskal_p_adj <- p.adjust(kruskal_res$p.value, method = "fdr")
    
    # If significant, perform Dunn's test and expand rows
    if(kruskal_res$p.value <= 0.05) {
      dunn_res <- dunn.test::dunn.test(
        x = test_data$value,
        g = test_data$group,
        method = "bh",
        list = TRUE
      )
      
      # Create rows for each comparison
      dunn_rows <- tibble(
        comparison = dunn_res$comparisons,
        dunn_z = dunn_res$Z,
        dunn_p = dunn_res$P,
        dunn_p_adj = p.adjust(dunn_res$P, method = "fdr")
      )
      
      # Combine with Kruskal results
      result_row <- result_row %>%
        slice(rep(1, nrow(dunn_rows))) %>%
        mutate(
          comparison = dunn_rows$comparison,
          dunn_z = dunn_rows$dunn_z,
          dunn_p = dunn_rows$dunn_p,
          dunn_p_adj = dunn_rows$dunn_p_adj
        )
    }
    
    result_row
  })
  
  return(results)
}

# Process all datasets
process_all_datasets <- function() {
  datasets <- list(
    complete_df_Missense,
    complete_df_Missense_canonical,
    complete_df_MPC_only,
    complete_df_PTV_canonic,
    complete_df_PTV
  )
  
  names <- c(
    "Missense",
    "Missense_canonical",
    "MPC_only",
    "PTV_canonic",
    "PTV"
  )
  
  # Run analysis for each dataset (both normalized and raw)
  all_results <- map2_dfr(datasets, names, function(df, name) {
    bind_rows( analyze_metrics_kruskal(df, normalize = FALSE, df_name = name)
    )
  })
  
  # Format final table
  final_table <- all_results %>%
    mutate(
      significance_kruskal = case_when(
        kruskal_p_adj < 0.001 ~ "***",
        kruskal_p_adj < 0.01 ~ "**",
        kruskal_p_adj < 0.05 ~ "*",
        TRUE ~ "ns"
      ),
      significance_dunn = case_when(
        is.na(dunn_p_adj) ~ NA_character_,
        dunn_p_adj < 0.001 ~ "***",
        dunn_p_adj < 0.01 ~ "**",
        dunn_p_adj < 0.05 ~ "*",
        TRUE ~ "ns"
      )
    ) %>%
    select(
      dataset, metric, normalization,
      kruskal_statistic, kruskal_p, kruskal_p_adj, significance_kruskal,
      comparison, dunn_z, dunn_p, dunn_p_adj, significance_dunn,
      note
    )
  
  return(final_table)
}

# Generate the complete results table
results_table <- process_all_datasets()

# View the results
View(results_table)

# Optional: Split into Kruskal and Dunn components if needed
kruskal_only <- results_table %>% 
  filter(is.na(comparison)) %>% 
  select(-starts_with("dunn"), -comparison)

dunn_only <- results_table %>% 
  filter(!is.na(comparison))

setwd("/home/rachele/SNVs/results_pasteur/group3/private")

# Save the results to a CSV file
write.csv(results_table, file = "results_table_kruskal_SCZ.csv", row.names = FALSE)
write.csv(kruskal_only, file = "kruskal_only_results_SCZ.csv", row.names = FALSE)
write.csv(dunn_only, file = "dunn_only_results_SCZ.csv", row.names = FALSE)

#plots 

# Function to prepare data for plotting
prepare_plot_data <- function(results_table, datasets) {
    # Create mapping between clean names and dataframe objects
    dataset_mapping <- list(
        "Missense" = complete_df_Missense,
        "Missense_canonical" = complete_df_Missense_canonical,
        "MPC_only" = complete_df_MPC_only,
        "PTV_canonic" = complete_df_PTV_canonic,  
        "PTV" = complete_df_PTV
    )
    
    # Prepare plot data
    plot_data <- map_dfr(names(dataset_mapping), function(ds_name) {
        df <- dataset_mapping[[ds_name]] %>%
            filter(Status %in% c("FEP", "Converter", "Non_Converter")) %>%
            select(sample_id, Status, starts_with("Number_of")) %>%
            pivot_longer(
                cols = -c(sample_id, Status),
                names_to = "metric",
                values_to = "count"
            ) %>%
            mutate(dataset = ds_name)
        
        # Add normalization information
        bind_rows(
            df %>% mutate(normalization = "FALSE"),
            df %>% 
                group_by(Status) %>%
                mutate(count = count/case_when(
                    Status == "FEP" ~ 255,
                    Status == "Converter" ~ 47,
                    Status == "Non_Converter" ~ 75
                )) %>%
                ungroup() %>%
                mutate(normalization = "TRUE")
        )
    })
    
    return(plot_data)
}

# Generate plot data
all_plot_data <- prepare_plot_data(results_table, list(
    "Missense", "Missense_canonical", "MPC_only", "PTV_canonic", "PTV"
))
# maybe not necessary 
# results_table$metric <- gsub("_per_indiv", "", results_table$metric)

# Function to create boxplot with Dunn's test results
create_boxplot_with_dunn <- function(plot_data, results_table, dataset_name, metric_name, is_normalized) {
  # Filter the specific data we want to plot
  current_data <- plot_data %>%
    filter(dataset == dataset_name,
           metric == metric_name,
           normalization == is_normalized)
  
  # Get the matching statistical results
  stats <- results_table %>%
    filter(dataset == dataset_name,
           metric == metric_name,
           normalization == is_normalized)
  
  # Create base plot
  p <- ggplot(current_data, aes(x = Status, y = count, fill = Status)) +
    geom_boxplot() +
    scale_fill_brewer(palette = "Set3") +
    labs(title = paste(dataset_name, "-", metric_name),
         subtitle = ifelse(is_normalized, "(Normalized)", "(Raw Counts)"),
         x = "Status",
         y = ifelse(is_normalized, "Count per individual", "Raw Count")) +
    theme_minimal() +
    theme(legend.position = "none")
  
  # Add Kruskal-Wallis annotation if available
  if (!is.na(stats$kruskal_p_adj[1])) {
    p <- p + annotate("text", x = 1, y = Inf, 
                     label = paste("K-W p =", formatC(stats$kruskal_p_adj[1], format = "e", digits = 2)), 
                     vjust = 2, hjust = -0.1)
    
    # Add Dunn's test results if significant and available
    if (stats$kruskal_p_adj[1] < 0.05 && 
        "comparison" %in% names(stats) && 
        !is.na(stats$comparison[1])) {
      
      # Prepare Dunn test annotations
      dunn_results <- stats %>%
        filter(!is.na(comparison)) %>%
        mutate(
          comparison = str_split(comparison, " - "),
          stars = case_when(
            dunn_p_adj < 0.001 ~ "***",
            dunn_p_adj < 0.01 ~ "**",
            dunn_p_adj < 0.05 ~ "*",
            TRUE ~ "ns"
          ),
          y_pos = max(current_data$count, na.rm = TRUE) * (1 + 0.1 * row_number())
        )
      
      # Add each comparison to the plot
      for (i in 1:nrow(dunn_results)) {
        p <- p + geom_signif(
          comparisons = list(dunn_results$comparison[[i]]),
          annotations = dunn_results$stars[i],
          y_position = dunn_results$y_pos[i],
          tip_length = 0.01
        )
      }
    }
  }
  
  return(p)
}

# Wrapper function to create all plots
create_all_plots <- function(plot_data, results_table) {
  # Get unique combinations of parameters
  combinations <- plot_data %>%
    distinct(dataset, metric, normalization)
  
  # Create a plot for each combination
  plots <- pmap(combinations, function(dataset, metric, normalization) {
    create_boxplot_with_dunn(
      plot_data = plot_data,
      results_table = results_table,
      dataset_name = dataset,
      metric_name = metric,
      is_normalized = normalization
    )
  })
  
  return(plots)
}

# Generate all plots
all_plots <- create_all_plots(all_plot_data, results_table)


# Create just one specific plot
single_plot <- create_boxplot_with_dunn(
  plot_data = all_plot_data,
  results_table = results_table,
  dataset_name = "Missense",
  metric_name = "Number_of_Variants",
  is_normalized = FALSE
)
print(single_plot)


non_normalized_plots <- all_plot_data %>%
  filter(normalization == "FALSE") %>%
  group_by(dataset, metric) %>%
  group_map(~ {
    create_boxplot_with_dunn(
      plot_data = all_plot_data,
      results_table = results_table,
      dataset_name = .y$dataset,
      metric_name = .y$metric,
      is_normalized = FALSE
    )
  }, .keep = TRUE)


# Set parameters
metrics_per_page <- 4
output_width <- 14  # inches
output_height <- 10 # inches

# Create a list to store all grids before saving
all_grids_non_normalized <- list()

for(current_dataset in unique(all_plot_data$dataset)) {
  
  # Get all metrics for current dataset (non-normalized)
  all_metrics <- all_plot_data %>%
    filter(dataset == current_dataset, normalization == "FALSE") %>%
    distinct(metric) %>%
    pull()
  
  # Split into chunks of 4 metrics
  metric_chunks <- split(
    all_metrics,
    ceiling(seq_along(all_metrics)/metrics_per_page)
  )
  
  # Store grids for this dataset
  dataset_grids <- list()
  
  for(page_num in seq_along(metric_chunks)) {
    
    # Create the 4 plots for this page
    page_plots <- lapply(metric_chunks[[page_num]], function(m) {
      create_boxplot_with_dunn(
        plot_data = all_plot_data,
        results_table = results_table,
        dataset_name = current_dataset,
        metric_name = m,
        is_normalized = FALSE
      ) +
        ggtitle(paste("Metric:", m)) +
        theme(plot.title = element_text(size = 10))
    })
    
    # Create and store the grid
    current_grid <- grid.arrange(
      grobs = page_plots,
      ncol = 2,
      nrow = 2,
      top = paste(current_dataset, "- Page", page_num)
    )
    
    # Display in RStudio Viewer
    print(current_grid)
    Sys.sleep(2)  # Brief pause for viewing
    
    # Store grid with unique name
    grid_name <- paste(current_dataset, "Page", page_num)
    dataset_grids[[page_num]] <- current_grid
    names(dataset_grids)[page_num] <- grid_name
  }
  
  # Add all pages for this dataset to master list
  all_grids_non_normalized[[current_dataset]] <- dataset_grids
}

# Save all grids after creation (external from the loop)
save_grids <- function(grid_list, format = "png") {
  for(dataset_name in names(grid_list)) {
    for(page_name in names(grid_list[[dataset_name]])) {
      filename <- paste0(dataset_name,"_","RAW" ,"_", gsub(" ", "_", page_name), ".", format)
      ggsave(
        filename = filename,
        plot = grid_list[[dataset_name]][[page_name]],
        width = output_width,
        height = output_height,
        units = "in",
        dpi = 300
      )
      message("Saved: ", filename)
    }
  }
}

# Execute saving (choose format)
save_grids(all_grids_non_normalized, format = "png")  # Can also use "pdf" for multi-page



#save only significant 
# Modified function to create all plots (only significant ones)
create_significant_plots <- function(plot_data, results_table) {
  # Get unique combinations of parameters with significant p-values
  significant_combinations <- results_table %>%
    filter(kruskal_p_adj < 0.05) %>%
    distinct(dataset, metric, normalization)
  
  # Create a plot for each significant combination
  plots <- pmap(significant_combinations, function(dataset, metric, normalization) {
    create_boxplot_with_dunn(
      plot_data = plot_data,
      results_table = results_table,
      dataset_name = dataset,
      metric_name = metric,
      is_normalized = normalization
    )
  })
  
  return(plots)
}


# Set parameters
metrics_per_page <- 4
output_width <- 14  # inches
output_height <- 10 # inches

# Create a list to store all significant grids before saving
all_significant_grids <- list()

for(current_dataset in unique(all_plot_data$dataset)) {
  # Get significant metrics for current dataset (non-normalized)
  significant_metrics <- results_table %>%
    filter(dataset == current_dataset, 
           normalization == "FALSE",
           kruskal_p_adj < 0.05) %>%
    distinct(metric) %>%
    pull()
  
  # Skip if no significant metrics for this dataset
  if(length(significant_metrics) == 0) next
  
  # Split into chunks of 4 metrics
  metric_chunks <- split(
    significant_metrics,
    ceiling(seq_along(significant_metrics)/metrics_per_page)
  )
  
  # Store grids for this dataset
  dataset_grids <- list()
  
  for(page_num in seq_along(metric_chunks)) {
    # Create the 4 plots for this page
    page_plots <- lapply(metric_chunks[[page_num]], function(m) {
      create_boxplot_with_dunn(
        plot_data = all_plot_data,
        results_table = results_table,
        dataset_name = current_dataset,
        metric_name = m,
        is_normalized = FALSE
      ) +
        ggtitle(paste("Metric:", m)) +
        theme(plot.title = element_text(size = 10))
    })
    
    # Create and store the grid
    current_grid <- grid.arrange(
      grobs = page_plots,
      ncol = 2,
      nrow = 2,
      top = paste(current_dataset, "- Page", page_num, "(Significant results only)")
    )
    
    # Display in RStudio Viewer
    print(current_grid)
    Sys.sleep(2)  # Brief pause for viewing
    
    # Store grid with unique name
    grid_name <- paste("Significant_",current_dataset, "Page", page_num)
    dataset_grids[[page_num]] <- current_grid
    names(dataset_grids)[page_num] <- grid_name
  }
  
  # Add all pages for this dataset to master list
  all_significant_grids[[current_dataset]] <- dataset_grids
}

# Save only significant grids
save_grids(all_significant_grids, format = "png")



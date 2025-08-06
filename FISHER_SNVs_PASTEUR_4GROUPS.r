# fisher_sNVs_PASTEUR_4GROUPS





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
manifest_correct$Status <- gsub("FEP-SCZ", "SCZ", manifest_correct$Status)
manifest_correct$Status <- gsub("FEP-BD", "BD", manifest_correct$Status)




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
    allowed <- c("SCZ", "BD", "Converter", "Non_Converter")
    
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







# CREATE DATAFRAME TO BE USED FOR STATISTICS 
    # Define a function to process each dataframe and return the result
    # This function processes the data to separate rows, clean up sample IDs, and filter based on specific gene sets
    process_df <- function(df, name, brain_gene_consensus_filtered_consensus_no_pitular, brain_gene_consensus_ntm_consensus_no_pitular, genes_schema_pval, genes_bipolar, genes_sfari_non_syndromic, genes_schema_or,private= FALSE) {
            
            
    if(private){
        df <- df %>% filter(count == 1)  # Only filter if private=TRUE
    }
    
    #dataframe to count individuals 
    expanded_df <- df %>%
        separate_rows(SAMPLES, sep = ",") %>%  # Split SAMPLES only
        mutate(
        sample_id = trimws(SAMPLES),
        group = sample_label_3  # Use pre-computed labels
        ) %>%
        filter(group %in% c("SCZ", "Non_Converter","BD","Converter")) %>%  # Remove mixed/NA
        filter(!sample_id %in% UHR_NA_SAMPLE_IDS$V1) %>%  # Remove UHR_NA
        filter(!sample_id %in% Low_Quality_SAMPLES_AND_UHR_NA$Sequencing_number) %>%
        distinct(SYMBOL, CHROM,REF,ALT, POS, sample_id, .keep_all = TRUE)  # Unique variant/sample

    #dataframe for count variants 
    variant_df <- df %>% 
    mutate(group = sample_label_3) %>%  # Use pre-computed labels
    filter(group %in% c("SCZ", "Non_Converter","BD","Converter")) %>%  # Remove mixed/NA
  
    distinct(SYMBOL, CHROM,REF,ALT, POS , .keep_all = TRUE)  # Unique variant/sample 
    
    # Filter based on the genes in the brain, brain filter_ntpm, and additional gene sets
    expanded_df <- expanded_df %>%
        mutate(brain = ifelse(SYMBOL %in% brain_gene_consensus_filtered_consensus_no_pitular$Gene.name, 1, 0),
            brain_filter_ntpm = ifelse(SYMBOL %in% brain_gene_consensus_ntm_consensus_no_pitular$Gene.name, 1, 0),
            schema_pval = ifelse(SYMBOL %in% genes_schema_pval, 1, 0),
            bipolar = ifelse(SYMBOL %in% genes_bipolar, 1, 0),
            sfari_non_syndromic = ifelse(SYMBOL %in% genes_sfari_non_syndromic, 1, 0),
            schema_or = ifelse(SYMBOL %in% genes_schema_or, 1, 0))

    variant_df <- variant_df %>%
            mutate(brain = ifelse(SYMBOL %in% brain_gene_consensus_filtered_consensus_no_pitular$Gene.name, 1, 0),
            brain_filter_ntpm = ifelse(SYMBOL %in% brain_gene_consensus_ntm_consensus_no_pitular$Gene.name, 1, 0),
            schema_pval = ifelse(SYMBOL %in% genes_schema_pval, 1, 0),
            bipolar = ifelse(SYMBOL %in% genes_bipolar, 1, 0),
            sfari_non_syndromic = ifelse(SYMBOL %in% genes_sfari_non_syndromic, 1, 0),
            schema_or = ifelse(SYMBOL %in% genes_schema_or, 1, 0))
    
    num_non_pathogenic_pLi_variants <- variant_df %>%
        group_by(group) %>%
        summarize(num_variants = n(), .groups = "drop")
    
    num_individuals_unique_with_non_pathogenic_pLi <- expanded_df %>%
        group_by(group) %>%
        summarize(unique_individuals = n_distinct(sample_id), .groups = "drop")

    num_individuals_with_non_pathogenic_pLi <- expanded_df %>%
        group_by(group) %>%  
        summarize(num_individuals = n(), .groups = "drop")
    
    num_brain_variants <- variant_df %>%
        filter(brain == 1) %>%
        group_by(group) %>%
        summarize(num_brain_variants = n(), .groups = "drop")
    
    num_individuals_unique_with_brain_variants <- expanded_df %>%
        filter(brain == 1) %>%
        group_by(group) %>%
        summarize(unique_individuals_brain = n_distinct(sample_id), .groups = "drop")

    num_individuals_with_brain_variants <- expanded_df %>%
        filter(brain == 1) %>%
        group_by(group) %>%
        summarize(num_individuals_brain = n(), .groups = "drop")  
    
    num_brain_filter_ntpm_variants <- variant_df %>%
        filter(brain_filter_ntpm == 1) %>%
        group_by(group) %>%
        summarize(num_brain_filter_ntpm_variants = n(), .groups = "drop")
    
    num_individuals_unique_with_brain_filter_ntpm_variants <- expanded_df %>%
        filter(brain_filter_ntpm == 1) %>%
        group_by(group) %>%
        summarize(unique_individuals_brain_filter_ntpm = n_distinct(sample_id), .groups = "drop")

    num_individuals_with_brain_filter_ntpm_variants <- expanded_df %>%
    filter(brain_filter_ntpm == 1) %>%
        group_by(group) %>%
        summarize(num_individuals_brain_filter_ntpm = n(), .groups = "drop")
    
    num_schema_pval_variants <- variant_df %>%
        filter(schema_pval == 1) %>%
        group_by(group) %>%
        summarize(num_schema_pval_variants = n(), .groups = "drop")
    
    num_individuals_unique_with_schema_pval_variants <- expanded_df %>%
        filter(schema_pval == 1) %>%
        group_by(group) %>%
        summarize(unique_individuals_schema_pval = n_distinct(sample_id), .groups = "drop")

    num_individuals_with_schema_pval_variants <- expanded_df %>%
        filter(schema_pval == 1) %>%
        group_by(group) %>%
        summarize(num_individuals_schema_pval = n(), .groups = "drop")  
    
    num_bipolar_variants <- variant_df %>%
        filter(bipolar == 1) %>%
        group_by(group) %>%
        summarize(num_bipolar_variants = n(), .groups = "drop")
    
    num_individuals_unique_with_bipolar_variants <- expanded_df %>%
        filter(bipolar == 1) %>%
        group_by(group) %>%
        summarize(unique_individuals_bipolar = n_distinct(sample_id), .groups = "drop")
    
    num_individuals_with_bipolar_variants <- expanded_df %>%
        filter(bipolar == 1) %>%
        group_by(group) %>%
        summarize(num_individuals_bipolar = n(), .groups = "drop")
    
    num_sfari_non_syndromic_variants <- variant_df %>%
        filter(sfari_non_syndromic == 1) %>%
        group_by(group) %>%
        summarize(num_sfari_non_syndromic_variants = n(), .groups = "drop")
    
    num_individuals_unique_with_sfari_non_syndromic_variants <- expanded_df %>%
        filter(sfari_non_syndromic == 1) %>%
        group_by(group) %>%
        summarize(unique_individuals_sfari_non_syndromic = n_distinct(sample_id), .groups = "drop")

    num_individuals_with_sfari_non_syndromic_variants <- expanded_df %>%
        filter(sfari_non_syndromic == 1) %>%
        group_by(group) %>%
        summarize(num_individuals_sfari_non_syndromic = n(), .groups = "drop")  
    
    num_schema_or_variants <- variant_df %>%
        filter(schema_or == 1) %>%
        group_by(group) %>%
        summarize(num_schema_or_variants = n(), .groups = "drop")
    
    num_individuals_unique_with_schema_or_variants <- expanded_df %>%
        filter(schema_or == 1) %>%
        group_by(group) %>%
        summarize(unique_individuals_schema_or = n_distinct(sample_id), .groups = "drop")
    
    num_individuals_with_schema_or_variants <- expanded_df %>%
        filter(schema_or == 1) %>%
        group_by(group) %>%
        summarize(num_individuals_schema_or = n(), .groups = "drop")
    
    # Ensure unique genes are counted
    num_genes <- expanded_df %>%
        group_by(group) %>%
        summarize(num_genes = n_distinct(SYMBOL), .groups = "drop")
    
    num_brain_genes <- expanded_df %>%
        filter(brain == 1) %>%
        group_by(group) %>%
        summarize(num_brain_genes = n_distinct(SYMBOL), .groups = "drop")
    
    num_brain_filter_ntpm_genes <- expanded_df %>%
        filter(brain_filter_ntpm == 1) %>%
        group_by(group) %>%
        summarize(num_brain_filter_ntpm_genes = n_distinct(SYMBOL), .groups = "drop")
    
    num_schema_pval_genes <- expanded_df %>%
        filter(schema_pval == 1) %>%
        group_by(group) %>%
        summarize(num_schema_pval_genes = n_distinct(SYMBOL), .groups = "drop")
    
    num_bipolar_genes <- expanded_df %>%
        filter(bipolar == 1) %>%
        group_by(group) %>%
        summarize(num_bipolar_genes = n_distinct(SYMBOL), .groups = "drop")
    
    num_sfari_non_syndromic_genes <- expanded_df %>%
        filter(sfari_non_syndromic == 1) %>%
        group_by(group) %>%
        summarize(num_sfari_non_syndromic_genes = n_distinct(SYMBOL), .groups = "drop")
    
    num_schema_or_genes <- expanded_df %>%
        filter(schema_or == 1) %>%
        group_by(group) %>%
        summarize(num_schema_or_genes = n_distinct(SYMBOL), .groups = "drop")

    # Ensure all groups have the same number of rows by filling missing values with zeros
    groups <- c("SCZ", "Non_Converter","BD","Converter")
    num_non_pathogenic_pLi_variants <- num_non_pathogenic_pLi_variants %>% complete(group = groups, fill = list(num_variants = 0))
    num_individuals_with_non_pathogenic_pLi <- num_individuals_with_non_pathogenic_pLi %>% complete(group = groups, fill = list(unique_individuals = 0))
    num_brain_variants <- num_brain_variants %>% complete(group = groups, fill = list(num_brain_variants = 0))
    num_individuals_with_brain_variants <- num_individuals_with_brain_variants %>% complete(group = groups, fill = list(unique_individuals_brain = 0))
    num_brain_filter_ntpm_variants <- num_brain_filter_ntpm_variants %>% complete(group = groups, fill = list(num_brain_filter_ntpm_variants = 0))
    num_individuals_with_brain_filter_ntpm_variants <- num_individuals_with_brain_filter_ntpm_variants %>% complete(group = groups, fill = list(unique_individuals_brain_filter_ntpm = 0))
    num_schema_pval_variants <- num_schema_pval_variants %>% complete(group = groups, fill = list(num_schema_pval_variants = 0))
    num_individuals_with_schema_pval_variants <- num_individuals_with_schema_pval_variants %>% complete(group = groups, fill = list(unique_individuals_schema_pval = 0))
    num_bipolar_variants <- num_bipolar_variants %>% complete(group = groups, fill = list(num_bipolar_variants = 0))
    num_individuals_with_bipolar_variants <- num_individuals_with_bipolar_variants %>% complete(group = groups, fill = list(unique_individuals_bipolar = 0))
    num_sfari_non_syndromic_variants <- num_sfari_non_syndromic_variants %>% complete(group = groups, fill = list(num_sfari_non_syndromic_variants = 0))
    num_individuals_with_sfari_non_syndromic_variants <- num_individuals_with_sfari_non_syndromic_variants %>% complete(group = groups, fill = list(unique_individuals_sfari_non_syndromic = 0))
    num_schema_or_variants <- num_schema_or_variants %>% complete(group = groups, fill = list(num_schema_or_variants = 0))
    num_individuals_with_schema_or_variants <- num_individuals_with_schema_or_variants %>% complete(group = groups, fill = list(unique_individuals_schema_or = 0))
    num_genes <- num_genes %>% complete(group = groups, fill = list(num_genes = 0))
    num_brain_genes <- num_brain_genes %>% complete(group = groups, fill = list(num_brain_genes = 0))
    num_brain_filter_ntpm_genes <- num_brain_filter_ntpm_genes %>% complete(group = groups, fill = list(num_brain_filter_ntpm_genes = 0))
    num_schema_pval_genes <- num_schema_pval_genes %>% complete(group = groups, fill = list(num_schema_pval_genes = 0))
    num_bipolar_genes <- num_bipolar_genes %>% complete(group = groups, fill = list(num_bipolar_genes = 0))
    num_sfari_non_syndromic_genes <- num_sfari_non_syndromic_genes %>% complete(group = groups, fill = list(num_sfari_non_syndromic_genes = 0))
    num_schema_or_genes <- num_schema_or_genes %>% complete(group = groups, fill = list(num_schema_or_genes = 0))

    num_individuals_unique_with_non_pathogenic_pLi <- num_individuals_unique_with_non_pathogenic_pLi %>% complete(group = groups, fill = list(unique_individuals = 0))
    num_individuals_unique_with_brain_variants <- num_individuals_unique_with_brain_variants %>% complete(group = groups, fill = list(unique_individuals_brain = 0))
    num_individuals_unique_with_brain_filter_ntpm_variants <- num_individuals_unique_with_brain_filter_ntpm_variants %>% complete(group = groups, fill = list(unique_individuals_brain_filter_ntpm = 0))
    num_individuals_unique_with_schema_pval_variants <- num_individuals_unique_with_schema_pval_variants %>% complete(group = groups, fill = list(unique_individuals_schema_pval = 0))
    num_individuals_unique_with_bipolar_variants <- num_individuals_unique_with_bipolar_variants %>% complete(group = groups, fill = list(unique_individuals_bipolar = 0))
    num_individuals_unique_with_sfari_non_syndromic_variants <- num_individuals_unique_with_sfari_non_syndromic_variants %>% complete(group = groups, fill = list(unique_individuals_sfari_non_syndromic = 0))
    num_individuals_unique_with_schema_or_variants <- num_individuals_unique_with_schema_or_variants %>% complete(group = groups, fill = list(unique_individuals_schema_or = 0))


    # Define the row names based on the dataframe's name
    row_names <- c(
        paste("Number of", name ,"variants"),
        paste("Number of individuals with",name ,"variants"),
        paste("Number of individuals unique with",name ,"variants"),
        paste("Number of brain", name ,"variants"),
        paste("Number of individuals with brain",name ,"variants"),
        paste("Number of individuals unique with brain",name ,"variants"),
        paste("Number of brain filter_ntpm", name ,"variants"),
        paste("Number of individuals with brain filter_ntpm",name ,"variants"),
        paste("Number of individuals unique with brain filter_ntpm",name ,"variants"),
        paste("Number of schema_pval", name ,"variants"),
        paste("Number of individuals with schema_pval",name ,"variants"),
        paste("Number of individuals unique with schema_pval",name ,"variants"),
        paste("Number of bipolar", name ,"variants"),
        paste("Number of individuals with bipolar",name ,"variants"),
        paste("Number of individuals unique with bipolar",name ,"variants"),
        paste("Number of sfari_non_syndromic", name ,"variants"),
        paste("Number of individuals with sfari_non_syndromic",name ,"variants"),
        paste("Number of individuals unique with sfari_non_syndromic",name ,"variants"),
        paste("Number of schema_or", name ,"variants"),
        paste("Number of individuals with schema_or",name ,"variants"),
        paste("Number of individuals unique with schema_or",name ,"variants"),
        paste("Number of genes in", name),
        paste("Number of brain genes in", name),
        paste("Number of brain filter_ntpm genes in", name),
        paste("Number of schema_pval genes in", name),
        paste("Number of bipolar genes in", name),
        paste("Number of sfari_non_syndromic genes in", name),
        paste("Number of schema_or genes in", name)
    )
    
    # Create the result dataframe
    result_df <- data.frame(
        row_names = row_names,
        
        Converter = c(
            num_non_pathogenic_pLi_variants$num_variants[num_non_pathogenic_pLi_variants$group == "Converter"],
            num_individuals_with_non_pathogenic_pLi$num_individuals[num_individuals_with_non_pathogenic_pLi$group == "Converter"],  # Corrected column
            num_individuals_unique_with_non_pathogenic_pLi$unique_individuals[num_individuals_unique_with_non_pathogenic_pLi$group == "Converter"],
            num_brain_variants$num_brain_variants[num_brain_variants$group == "Converter"],
            num_individuals_with_brain_variants$num_individuals_brain[num_individuals_with_brain_variants$group == "Converter"],  # Corrected column
            num_individuals_unique_with_brain_variants$unique_individuals_brain[num_individuals_unique_with_brain_variants$group == "Converter"],
            num_brain_filter_ntpm_variants$num_brain_filter_ntpm_variants[num_brain_filter_ntpm_variants$group == "Converter"],
            num_individuals_with_brain_filter_ntpm_variants$num_individuals_brain_filter_ntpm[num_individuals_with_brain_filter_ntpm_variants$group == "Converter"],  # Corrected column
            num_individuals_unique_with_brain_filter_ntpm_variants$unique_individuals_brain_filter_ntpm[num_individuals_unique_with_brain_filter_ntpm_variants$group == "Converter"],
            num_schema_pval_variants$num_schema_pval_variants[num_schema_pval_variants$group == "Converter"],
            num_individuals_with_schema_pval_variants$num_individuals_schema_pval[num_individuals_with_schema_pval_variants$group == "Converter"],  # Corrected column
            num_individuals_unique_with_schema_pval_variants$unique_individuals_schema_pval[num_individuals_unique_with_schema_pval_variants$group == "Converter"],
            num_bipolar_variants$num_bipolar_variants[num_bipolar_variants$group == "Converter"],
            num_individuals_with_bipolar_variants$num_individuals_bipolar[num_individuals_with_bipolar_variants$group == "Converter"],  # Corrected column
            num_individuals_unique_with_bipolar_variants$unique_individuals_bipolar[num_individuals_unique_with_bipolar_variants$group == "Converter"],
            num_sfari_non_syndromic_variants$num_sfari_non_syndromic_variants[num_sfari_non_syndromic_variants$group == "Converter"],
            num_individuals_with_sfari_non_syndromic_variants$num_individuals_sfari_non_syndromic[num_individuals_with_sfari_non_syndromic_variants$group == "Converter"],  # Corrected column
            num_individuals_unique_with_sfari_non_syndromic_variants$unique_individuals_sfari_non_syndromic[num_individuals_unique_with_sfari_non_syndromic_variants$group == "Converter"],
            num_schema_or_variants$num_schema_or_variants[num_schema_or_variants$group == "Converter"],
            num_individuals_with_schema_or_variants$num_individuals_schema_or[num_individuals_with_schema_or_variants$group == "Converter"],  # Corrected column
            num_individuals_unique_with_schema_or_variants$unique_individuals_schema_or[num_individuals_unique_with_schema_or_variants$group == "Converter"],
            num_genes$num_genes[num_genes$group == "Converter"],
            num_brain_genes$num_brain_genes[num_brain_genes$group == "Converter"],
            num_brain_filter_ntpm_genes$num_brain_filter_ntpm_genes[num_brain_filter_ntpm_genes$group == "Converter"],
            num_schema_pval_genes$num_schema_pval_genes[num_schema_pval_genes$group == "Converter"],
            num_bipolar_genes$num_bipolar_genes[num_bipolar_genes$group == "Converter"],
            num_sfari_non_syndromic_genes$num_sfari_non_syndromic_genes[num_sfari_non_syndromic_genes$group == "Converter"],
            num_schema_or_genes$num_schema_or_genes[num_schema_or_genes$group == "Converter"]
        ),
        
        Non_Converter = c(
            num_non_pathogenic_pLi_variants$num_variants[num_non_pathogenic_pLi_variants$group == "Non_Converter"],
            num_individuals_with_non_pathogenic_pLi$num_individuals[num_individuals_with_non_pathogenic_pLi$group == "Non_Converter"],  # Corrected column
            num_individuals_unique_with_non_pathogenic_pLi$unique_individuals[num_individuals_unique_with_non_pathogenic_pLi$group == "Non_Converter"],
            num_brain_variants$num_brain_variants[num_brain_variants$group == "Non_Converter"],
            num_individuals_with_brain_variants$num_individuals_brain[num_individuals_with_brain_variants$group == "Non_Converter"],  # Corrected column
            num_individuals_unique_with_brain_variants$unique_individuals_brain[num_individuals_unique_with_brain_variants$group == "Non_Converter"],
            num_brain_filter_ntpm_variants$num_brain_filter_ntpm_variants[num_brain_filter_ntpm_variants$group == "Non_Converter"],
            num_individuals_with_brain_filter_ntpm_variants$num_individuals_brain_filter_ntpm[num_individuals_with_brain_filter_ntpm_variants$group == "Non_Converter"],  # Corrected column
            num_individuals_unique_with_brain_filter_ntpm_variants$unique_individuals_brain_filter_ntpm[num_individuals_unique_with_brain_filter_ntpm_variants$group == "Non_Converter"],
            num_schema_pval_variants$num_schema_pval_variants[num_schema_pval_variants$group == "Non_Converter"],
            num_individuals_with_schema_pval_variants$num_individuals_schema_pval[num_individuals_with_schema_pval_variants$group == "Non_Converter"],  # Corrected column
            num_individuals_unique_with_schema_pval_variants$unique_individuals_schema_pval[num_individuals_unique_with_schema_pval_variants$group == "Non_Converter"],
            num_bipolar_variants$num_bipolar_variants[num_bipolar_variants$group == "Non_Converter"],
            num_individuals_with_bipolar_variants$num_individuals_bipolar[num_individuals_with_bipolar_variants$group == "Non_Converter"],  # Corrected column
            num_individuals_unique_with_bipolar_variants$unique_individuals_bipolar[num_individuals_unique_with_bipolar_variants$group == "Non_Converter"],
            num_sfari_non_syndromic_variants$num_sfari_non_syndromic_variants[num_sfari_non_syndromic_variants$group == "Non_Converter"],
            num_individuals_with_sfari_non_syndromic_variants$num_individuals_sfari_non_syndromic[num_individuals_with_sfari_non_syndromic_variants$group == "Non_Converter"],  # Corrected column
            num_individuals_unique_with_sfari_non_syndromic_variants$unique_individuals_sfari_non_syndromic[num_individuals_unique_with_sfari_non_syndromic_variants$group == "Non_Converter"],
            num_schema_or_variants$num_schema_or_variants[num_schema_or_variants$group == "Non_Converter"],
            num_individuals_with_schema_or_variants$num_individuals_schema_or[num_individuals_with_schema_or_variants$group == "Non_Converter"],  # Corrected column
            num_individuals_unique_with_schema_or_variants$unique_individuals_schema_or[num_individuals_unique_with_schema_or_variants$group == "Non_Converter"],
            num_genes$num_genes[num_genes$group == "Non_Converter"],
            num_brain_genes$num_brain_genes[num_brain_genes$group == "Non_Converter"],
            num_brain_filter_ntpm_genes$num_brain_filter_ntpm_genes[num_brain_filter_ntpm_genes$group == "Non_Converter"],
            num_schema_pval_genes$num_schema_pval_genes[num_schema_pval_genes$group == "Non_Converter"],
            num_bipolar_genes$num_bipolar_genes[num_bipolar_genes$group == "Non_Converter"],
            num_sfari_non_syndromic_genes$num_sfari_non_syndromic_genes[num_sfari_non_syndromic_genes$group == "Non_Converter"],
            num_schema_or_genes$num_schema_or_genes[num_schema_or_genes$group == "Non_Converter"]

        ),
            SCZ = c(
            num_non_pathogenic_pLi_variants$num_variants[num_non_pathogenic_pLi_variants$group == "SCZ"],
            num_individuals_with_non_pathogenic_pLi$num_individuals[num_individuals_with_non_pathogenic_pLi$group == "SCZ"],  # Corrected column
            num_individuals_unique_with_non_pathogenic_pLi$unique_individuals[num_individuals_unique_with_non_pathogenic_pLi$group == "SCZ"],
            num_brain_variants$num_brain_variants[num_brain_variants$group == "SCZ"],
            num_individuals_with_brain_variants$num_individuals_brain[num_individuals_with_brain_variants$group == "SCZ"],  # Corrected column
            num_individuals_unique_with_brain_variants$unique_individuals_brain[num_individuals_unique_with_brain_variants$group == "SCZ"],
            num_brain_filter_ntpm_variants$num_brain_filter_ntpm_variants[num_brain_filter_ntpm_variants$group == "SCZ"],
            num_individuals_with_brain_filter_ntpm_variants$num_individuals_brain_filter_ntpm[num_individuals_with_brain_filter_ntpm_variants$group == "SCZ"],  # Corrected column
            num_individuals_unique_with_brain_filter_ntpm_variants$unique_individuals_brain_filter_ntpm[num_individuals_unique_with_brain_filter_ntpm_variants$group == "SCZ"],
            num_schema_pval_variants$num_schema_pval_variants[num_schema_pval_variants$group == "SCZ"],
            num_individuals_with_schema_pval_variants$num_individuals_schema_pval[num_individuals_with_schema_pval_variants$group == "SCZ"],  # Corrected column
            num_individuals_unique_with_schema_pval_variants$unique_individuals_schema_pval[num_individuals_unique_with_schema_pval_variants$group == "SCZ"],
            num_bipolar_variants$num_bipolar_variants[num_bipolar_variants$group == "SCZ"],
            num_individuals_with_bipolar_variants$num_individuals_bipolar[num_individuals_with_bipolar_variants$group == "SCZ"],  # Corrected column
            num_individuals_unique_with_bipolar_variants$unique_individuals_bipolar[num_individuals_unique_with_bipolar_variants$group == "SCZ"],
            num_sfari_non_syndromic_variants$num_sfari_non_syndromic_variants[num_sfari_non_syndromic_variants$group == "SCZ"],
            num_individuals_with_sfari_non_syndromic_variants$num_individuals_sfari_non_syndromic[num_individuals_with_sfari_non_syndromic_variants$group == "SCZ"],  # Corrected column
            num_individuals_unique_with_sfari_non_syndromic_variants$unique_individuals_sfari_non_syndromic[num_individuals_unique_with_sfari_non_syndromic_variants$group == "SCZ"],
            num_schema_or_variants$num_schema_or_variants[num_schema_or_variants$group == "SCZ"],
            num_individuals_with_schema_or_variants$num_individuals_schema_or[num_individuals_with_schema_or_variants$group == "SCZ"],  # Corrected column
            num_individuals_unique_with_schema_or_variants$unique_individuals_schema_or[num_individuals_unique_with_schema_or_variants$group == "SCZ"],
            num_genes$num_genes[num_genes$group == "SCZ"],
            num_brain_genes$num_brain_genes[num_brain_genes$group == "SCZ"],
            num_brain_filter_ntpm_genes$num_brain_filter_ntpm_genes[num_brain_filter_ntpm_genes$group == "SCZ"],
            num_schema_pval_genes$num_schema_pval_genes[num_schema_pval_genes$group == "SCZ"],
            num_bipolar_genes$num_bipolar_genes[num_bipolar_genes$group == "SCZ"],
            num_sfari_non_syndromic_genes$num_sfari_non_syndromic_genes[num_sfari_non_syndromic_genes$group == "SCZ"],
            num_schema_or_genes$num_schema_or_genes[num_schema_or_genes$group == "SCZ"]  

                ),


            BD = c(
            num_non_pathogenic_pLi_variants$num_variants[num_non_pathogenic_pLi_variants$group == "BD"],
            num_individuals_with_non_pathogenic_pLi$num_individuals[num_individuals_with_non_pathogenic_pLi$group == "BD"],  # Corrected column
            num_individuals_unique_with_non_pathogenic_pLi$unique_individuals[num_individuals_unique_with_non_pathogenic_pLi$group == "BD"],
            num_brain_variants$num_brain_variants[num_brain_variants$group == "BD"],
            num_individuals_with_brain_variants$num_individuals_brain[num_individuals_with_brain_variants$group == "BD"],  # Corrected column
            num_individuals_unique_with_brain_variants$unique_individuals_brain[num_individuals_unique_with_brain_variants$group == "BD"],
            num_brain_filter_ntpm_variants$num_brain_filter_ntpm_variants[num_brain_filter_ntpm_variants$group == "BD"],
            num_individuals_with_brain_filter_ntpm_variants$num_individuals_brain_filter_ntpm[num_individuals_with_brain_filter_ntpm_variants$group == "BD"],  # Corrected column
            num_individuals_unique_with_brain_filter_ntpm_variants$unique_individuals_brain_filter_ntpm[num_individuals_unique_with_brain_filter_ntpm_variants$group == "BD"],
            num_schema_pval_variants$num_schema_pval_variants[num_schema_pval_variants$group == "BD"],
            num_individuals_with_schema_pval_variants$num_individuals_schema_pval[num_individuals_with_schema_pval_variants$group == "BD"],  # Corrected column
            num_individuals_unique_with_schema_pval_variants$unique_individuals_schema_pval[num_individuals_unique_with_schema_pval_variants$group == "BD"],
            num_bipolar_variants$num_bipolar_variants[num_bipolar_variants$group == "BD"],
            num_individuals_with_bipolar_variants$num_individuals_bipolar[num_individuals_with_bipolar_variants$group == "BD"],  # Corrected column
            num_individuals_unique_with_bipolar_variants$unique_individuals_bipolar[num_individuals_unique_with_bipolar_variants$group == "BD"],
            num_sfari_non_syndromic_variants$num_sfari_non_syndromic_variants[num_sfari_non_syndromic_variants$group == "BD"],
            num_individuals_with_sfari_non_syndromic_variants$num_individuals_sfari_non_syndromic[num_individuals_with_sfari_non_syndromic_variants$group == "BD"],  # Corrected column
            num_individuals_unique_with_sfari_non_syndromic_variants$unique_individuals_sfari_non_syndromic[num_individuals_unique_with_sfari_non_syndromic_variants$group == "BD"],
            num_schema_or_variants$num_schema_or_variants[num_schema_or_variants$group == "BD"],
            num_individuals_with_schema_or_variants$num_individuals_schema_or[num_individuals_with_schema_or_variants$group == "BD"],  # Corrected column
            num_individuals_unique_with_schema_or_variants$unique_individuals_schema_or[num_individuals_unique_with_schema_or_variants$group == "BD"],
            num_genes$num_genes[num_genes$group == "BD"],
            num_brain_genes$num_brain_genes[num_brain_genes$group == "BD"],
            num_brain_filter_ntpm_genes$num_brain_filter_ntpm_genes[num_brain_filter_ntpm_genes$group == "BD"],
            num_schema_pval_genes$num_schema_pval_genes[num_schema_pval_genes$group == "BD"],
            num_bipolar_genes$num_bipolar_genes[num_bipolar_genes$group == "BD"],
            num_sfari_non_syndromic_genes$num_sfari_non_syndromic_genes[num_sfari_non_syndromic_genes$group == "BD"],
            num_schema_or_genes$num_schema_or_genes[num_schema_or_genes$group == "BD"]



        )
    )
    
    return(result_df)
    }


 # Process each dataframe
    result_df_list <- list(
    process_df(Missense, "Missense", brain_gene_consensus_filtered_consensus_no_pitular, brain_gene_consensus_ntm_consensus_no_pitular, genes_schema_pval, genes_bipolar, genes_sfari_non_syndromic, genes_schema_or),
    process_df(Missense_canonical, "Missense_canonical", brain_gene_consensus_filtered_consensus_no_pitular, brain_gene_consensus_ntm_consensus_no_pitular, genes_schema_pval, genes_bipolar, genes_sfari_non_syndromic, genes_schema_or),
    process_df(MPC_only, "MPC_only", brain_gene_consensus_filtered_consensus_no_pitular, brain_gene_consensus_ntm_consensus_no_pitular, genes_schema_pval, genes_bipolar, genes_sfari_non_syndromic, genes_schema_or),
    process_df(PTV_canonic, "PTV_canonic", brain_gene_consensus_filtered_consensus_no_pitular, brain_gene_consensus_ntm_consensus_no_pitular, genes_schema_pval, genes_bipolar, genes_sfari_non_syndromic, genes_schema_or),
    process_df(PTV, "PTV", brain_gene_consensus_filtered_consensus_no_pitular, brain_gene_consensus_ntm_consensus_no_pitular, genes_schema_pval, genes_bipolar, genes_sfari_non_syndromic, genes_schema_or)
    )

    # Combine all result dataframes into one
    final_result_df <- bind_rows(result_df_list)

    # View the final result
    print(final_result_df)
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


#when there is an NA put 0 
final_result_df[is.na(final_result_df)] <-0
result_df_Missense[is.na(result_df_Missense)] <- 0
result_df_Missense_canonical[is.na(result_df_Missense_canonical)] <- 0
result_df_MPC_only[is.na(result_df_MPC_only)] <- 0
result_df_PTV_canonic[is.na(result_df_PTV_canonic)] <- 0
result_df_PTV[is.na(result_df_PTV)] <- 0


# pairwise fisher individuals 
# Function to perform pairwise Fisher's exact tests between groups
perform_pairwise_fisher_per_ind <- function(group1_count, group2_count, total_group1, total_group2) {
    # Input validation
    if (group1_count > total_group1 || group2_count > total_group2) {
        stop("Carrier counts exceed group totals")
    }
    
    # Build contingency table
    a <- group1_count  # Group1 with variant
    b <- total_group1 - a  # Group1 without variant
    c <- group2_count  # Group2 with variant
    d <- total_group2 - c  # Group2 without variant
    
    contingency_table <- matrix(
        c(a, b, c, d), 
        nrow = 2,
        byrow = TRUE,
        dimnames = list(
            Group = c("Group1", "Group2"),
            Variant = c("Carrier", "NonCarrier")
        )
    )

    # Perform Fisher's exact test
    test <- fisher.test(contingency_table)
    
    return(list(
        p_value = test$p.value,
        odds_ratio = test$estimate,
        ci_low = test$conf.int[1],
        ci_high = test$conf.int[2],
        table = contingency_table
    ))
}

# Function to perform pairwise tests for all group combinations
perform_pairwise_tests_per_ind <- function(result_df) {
    # Define rows to test and their corresponding names
    rows <- c(3, 6, 9, 12, 15, 18, 21)
    row_names <- c("filtered", "brain_variants", "brain_filter_ntpm_variants",
                   "schema_pval_variants", "bipolar_variants", 
                   "sfari_non_syndromic_variants", "schema_or_variants")
    
    # Define group totals
    group_totals <- c(
        "SCZ" = 222,
        "BD" = 33,
        "Converter" = 47,
        "Non_Converter" = 75
    )
    
    # Generate all possible group pairs
    groups <- names(group_totals)
    group_pairs <- combn(groups, 2, simplify = FALSE)
    
    # Perform tests for each row and group pair
    results <- lapply(seq_along(rows), function(i) {
        row <- rows[i]
        row_name <- row_names[i]
        
        pair_results <- lapply(group_pairs, function(pair) {
            group1 <- pair[1]
            group2 <- pair[2]
            
            # Get counts and totals for this pair
            group1_count <- result_df[row, group1]
            group2_count <- result_df[row, group2]
            total_group1 <- group_totals[group1]
            total_group2 <- group_totals[group2]
            
            # Skip if any count is NA
            if (is.na(group1_count)) return(NULL)
            if (is.na(group2_count)) return(NULL)
            
            # Perform test
            res <- tryCatch({
                perform_pairwise_fisher_per_ind(
                    group1_count, group2_count,
                    total_group1, total_group2
                )
            }, error = function(e) {
                return(NULL)
            })
            
            if (is.null(res)) return(NULL)
            
            data.frame(
                Row = row_name,
                Group1 = group1,
                Group2 = group2,
                Group1_Carriers = res$table[1, 1],
                Group1_NonCarriers = res$table[1, 2],
                Group2_Carriers = res$table[2, 1],
                Group2_NonCarriers = res$table[2, 2],
                OR = round(res$odds_ratio, 2),
                CI = paste0("[", round(res$ci_low, 2), ", ", round(res$ci_high, 2), "]"),
                p_value = res$p_value,
                stringsAsFactors = FALSE
            )
        })
        
        # Remove NULL results and combine
        Filter(Negate(is.null), pair_results) %>% bind_rows()
    })
    
    bind_rows(results)
}

# Apply the function to each result dataframe and compile results
pairwise_results_ind <- list(
    Missense = perform_pairwise_tests_per_ind(result_df_Missense),
    Missense_canonical = perform_pairwise_tests_per_ind(result_df_Missense_canonical),
    MPC_only = perform_pairwise_tests_per_ind(result_df_MPC_only),
    PTV_canonic = perform_pairwise_tests_per_ind(result_df_PTV_canonic),  
    PTV = perform_pairwise_tests_per_ind(result_df_PTV)
)

# Combine all results into a single dataframe
pairwise_results_ind_df <- bind_rows(pairwise_results_ind, .id = "Test")

# Apply multiple testing corrections
# Format p-values
# Apply multiple testing corrections
pairwise_results_ind_df <- pairwise_results_ind_df %>%
    group_by(Test, Row) %>%  # Group by test type and variant category
    mutate(p_adj_per_test = p.adjust(p_value, method = "fdr")) %>%
    ungroup() %>%
    mutate(p_adj_global = p.adjust(p_value, method = "fdr"))

# View final results
print(pairwise_results_ind_df)


#rate comparison

# Function to calculate carrier rates for all groups
rate_confront <- function(result_df, row) {
    # Define group totals (from your specifications)
    group_totals <- c(
        SCZ = 222,
        BD = 33,
        Converter = 47,
        Non_Converter = 75
    )
    
    # Calculate rates for all groups
    rates <- lapply(names(group_totals), function(group) {
        result_df[row, group] / group_totals[group]
    })
    
    names(rates) <- paste0(names(group_totals), "_Rate")
    return(rates)
}

# Function to apply the calculation to specified rows
perform_tests_for_all_rows <- function(result_df) {
    # Define rows and their names
    rows <- c(1, 4, 7, 10, 13, 16, 19)
    row_names <- c("total", "brain_variants", "brain_filter_ntpm_variants", 
                 "schema_pval_variants", "bipolar_variants", 
                 "sfari_non_syndromic_variants", "schema_or_variants")

    # Calculate rates for all rows
    results <- lapply(seq_along(rows), function(i) {
        data.frame(
            Row = row_names[i],
            as.data.frame(rate_confront(result_df, rows[i])))
    })
    
    do.call(rbind, results)
}

# Run the analysis on all datasets
rate_confront <- list(
    Missense = perform_tests_for_all_rows(result_df_Missense),
    Missense_canonical = perform_tests_for_all_rows(result_df_Missense_canonical),
    MPC_only = perform_tests_for_all_rows(result_df_MPC_only),
    PTV_canonic = perform_tests_for_all_rows(result_df_PTV_canonic),
    PTV = perform_tests_for_all_rows(result_df_PTV)
)

# Combine results into a single dataframe
rate_confront_df <- bind_rows(rate_confront, .id = "Test") %>%
    select(Test, Row, SCZ_Rate, BD_Rate, Converter_Rate, Non_Converter_Rate)

# Print formatted results
print(rate_confront_df)


#visaluzation 
volcano_plot_ind <- pairwise_results_ind_df %>%
    mutate(log_OR = log2(OR),
           log_p = -log10(p_adj_per_test)) %>%
    ggplot(aes(x = log_OR, y = log_p)) +
    geom_point(aes(color = interaction(Group1, Group2), 
                   shape = p_adj_per_test < 0.05), 
               size = 3, alpha = 0.7) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), color = "red") +
    ggrepel::geom_text_repel(
        data = . %>% filter(p_adj_per_test < 0.05),
        aes(label = paste(Row, "\n", Group1, "vs", Group2)),
        max.overlaps = 20, size = 3
    ) +
    facet_grid(Test ~ ., scales = "free") +
    scale_color_viridis_d(option = "turbo") +
    labs(title = "Pairwise Variant Associations",
         x = "log2(Odds Ratio)",
         y = "-log10( FDR p-value)",
         color = "Group Comparison") +
    theme_bw() +
    theme(legend.position = "bottom")


volcano_plot_ind_all <- pairwise_results_ind_df %>%
  mutate(
    Comparison = paste(Group1, "vs", Group2),
    # Only label significant points
    Label = ifelse(p_adj_per_test < 0.05, Comparison, "")
  ) %>%
  ggplot(aes(x = log2(OR), y = -log10(p_adj_per_test))) +
  geom_point(aes(color = Row, shape = Test), size = 3, alpha = 0.7) +
  geom_label_repel(aes(label = Label),  # Auto-avoid overlaps
                   box.padding = 0.5,
                   max.overlaps = 20,
                   size = 2.5,
                   segment.color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  labs(title = "Volcano Plot: Effect Size vs. Significance",
       x = "log2(Odds Ratio)",
       y = "-log10(FDR-adjusted p-value)") +
  theme_bw()

matrix_plot_ind_OR <- pairwise_results_ind_df %>%
    mutate(Significance = case_when(
        p_adj_global < 0.001 ~ "***",
        p_adj_global < 0.01 ~ "**",
        p_adj_global < 0.05 ~ "*",
        TRUE ~ ""
    )) %>%
    ggplot(aes(x = Group1, y = Group2)) +
    geom_tile(aes(fill = log2(OR)), color = "white") +
    geom_text(aes(label = paste0(round(OR, 2), Significance)), 
              size = 3, color = "black") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                         midpoint = 0, limits = c(-2, 2)) +
    facet_grid(Test ~ Row, scales = "free", space = "free") +
    labs(title = "Odds Ratio Matrix by Group Comparison",
         fill = "log2(OR)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.spacing = unit(0.5, "lines"))  

parallel_plot_ind <- pairwise_results_ind_df %>%
    mutate(Comparison = paste(Group1, "vs", Group2),
           log_OR = log2(OR)) %>%
    ggplot(aes(x = Comparison, y = log_OR)) +
    geom_line(aes(group = interaction(Test, Row), color = Row), 
              alpha = 0.3) +
    geom_point(aes(color = Row, size = -log10(p_adj_global)), 
               alpha = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_boxplot(width = 0.2, outlier.shape = NA) +
    facet_wrap(~Test, ncol = 2) +
    scale_size_continuous(range = c(1, 4)) +
    labs(title = "Effect Size Patterns Across Comparisons",
         y = "log2(Odds Ratio)",
         size = "-log10(p-value)") +
    coord_flip() +
    theme_bw()

matrix_plot_ind_p <- pairwise_results_ind_df %>%
  mutate(
    Significance = case_when(
      p_adj_per_test < 0.001 ~ "***",
      p_adj_per_test < 0.01 ~ "**",
      p_adj_per_test < 0.05 ~ "*",
      TRUE ~ ""
    ),
    # Format p-values for display
    p_formatted = ifelse(p_adj_per_test < 0.001,
                         formatC(p_adj_per_test, format = "e", digits = 1),
                         round(p_adj_per_test, 3))
  ) %>%
  ggplot(aes(x = Group1, y = Group2)) +
  geom_tile(aes(fill = -log10(p_adj_per_test)), color = "white") +
  geom_text(aes(label = paste0(p_formatted, Significance)), 
            size = 3, color = "black") +
  scale_fill_gradient(low = "white", high = "red",
                      name = "-log10(FDR p-value)",
                      limits = c(0, 3)) +  # Adjust limits based on your data
  facet_grid(Test ~ Row, scales = "free", space = "free") +
  labs(title = "P-value Matrix by Group Comparison") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.spacing = unit(0.5, "lines")
  )


dot_plot_ind <- pairwise_results_ind_df %>%
  mutate(Comparison = paste(Group1, "vs", Group2)) %>%
  ggplot(aes(x = -log10(p_adj_per_test), y = Row)) +
  geom_point(aes(color = Comparison, shape = p_adj_per_test < 0.05), size = 3) +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "red") +
  scale_shape_manual(values = c(16, 17)) +  # 16 = non-significant, 17 = significant
  labs(
    title = "Significance by datset & Test Type",
    x = "-log10(FDR-adjusted p-value)",
    y = "Region",
    color = "Comparison",  # Changed from "Test Type" to match aesthetic
    shape = "Significant (FDR < 0.05)"
  ) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 8))


manhattan_plot_ind <- pairwise_results_ind_df %>%
  ggplot(aes(x = interaction(Group1, Group2), y = -log10(p_adj_per_test))) +
  geom_point(aes(color = Test), size = 2) +  # Color by test type
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  labs(
    title = "Manhattan Plot: Significance Across Comparisons",
    x = "Group Comparisons",
    y = "-log10(FDR-adjusted p-value)"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#view plots

volcano_plot_ind
volcano_plot_ind_all
dot_plot_ind
manhattan_plot_ind
matrix_plot_ind_OR
matrix_plot_ind_p


#set working directry
setwd("/home/rachele/SNVs/results_pasteur/group4/private")


# save plots 

ggsave("Fisher_volcano_plot_ind_SCZ.png", volcano_plot_ind, width = 8, height = 6)
ggsave("Fisher_volcano_plot_ind_all_SCZ.png", volcano_plot_ind_all, width = 8, height = 6)
ggsave("Fisher_matrix_plot_ind_OR_SCZ.png", matrix_plot_ind_OR, width = 8, height = 6)
ggsave("Fisher_matrix_plot_ind_p_SCZ.png", matrix_plot_ind_p, width = 8, height = 6)
ggsave("Fisher_parallel_plot_ind_SCZ.png", parallel_plot_ind, width = 8, height = 6)
ggsave("Fisher_dot_plot_ind_SCZ.png", dot_plot_ind, width = 8, height = 6)
ggsave("Fisher_manhattan_plot_ind_SCZ.png", manhattan_plot_ind, width = 8, height = 6)

#save final results
write.csv(pairwise_results_ind_df, "fisher_individuals_SCZ.csv", row.names = FALSE)
write.csv(rate_confront_df, "rate_confront_SCZ.csv", row.names = FALSE)





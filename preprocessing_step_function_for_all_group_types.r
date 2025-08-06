# Load necessary libraries for data manipulation and analysis
library(readr)
library(data.table)
library(dplyr)
library(gprofiler2)
library(httr)
library(readxl)
library(tidyr)
library(stringr)
library(ggplot2)

# ! LOAD GENESETS
# Load my data from CSV files
    SCHEMA <- read_csv("/home/rachele/Documents/old/geneset_confrontation/SCHEMA.csv")
    GWAS_120 <- read.delim("~/GWAS_120.csv")       
    BipEx_Bipolar <- read_csv("/home/rachele/Documents/old/geneset_confrontation/Bipolar_Disorder.csv")
    SFARI <- read_csv("/home/rachele/Documents/old/geneset_confrontation/SFARI.csv")

# Clean GENESETS
    # Convert gene names to a standard format for easier comparison
    # BipEx_Bipolar
    convert_col <- gconvert(query = BipEx_Bipolar$Gene, organism = "hsapiens", target = "ENSG", mthreshold = Inf, filter_na = TRUE)
    BipEx_Bipolar <- merge(BipEx_Bipolar, convert_col[, c("target", "name")], by.x = "Gene", by.y = "target", all.x = TRUE)
    BipEx_Bipolar <- BipEx_Bipolar %>% select(Gene, name, everything())

    # SCHEMA
    convert_col <- gconvert(query = SCHEMA$Gene, organism = "hsapiens", target = "ENSG", mthreshold = Inf, filter_na = TRUE)
    SCHEMA <- merge(SCHEMA, convert_col[, c("target", "name")], by.x = "Gene", by.y = "target", all.x = TRUE)
    SCHEMA <- SCHEMA %>% select(Gene, name, everything())

    # General clean-up of data
    # Remove unnecessary columns, rename columns, and clean up gene names
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

    # Combine all the dataframes into a single list of unique genes
    genes_schema_pval <- unique(SCHEMA_pVal$Gene)
    genes_bipolar <- unique(c(BipEx_Bipolar_p_val_PTV$Gene, BipEx_Bipolar_p_val_Missense$Gene))
    genes_sfari_non_syndromic <- unique(SFARI_non_syndromic_lower_3$Gene)
    genes_gwas <- unique(GWAS_120$Gene)
    genes_schema_or <- unique(SCHEMA_OR$Gene)
    genes_sfari_syndromic <- unique(SFARI_syndromic$Gene)

# ! load data 
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


# !preprocess variant data 
process_variant_data <- function(group_type, manifest, dfs) {
  # Clean SAMPLES column to remove UHR_NA/UHR-NA
  clean_samples <- function(df) {
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

  # Recode manifest status based on group type
  manifest <- manifest %>%
    mutate(Status = case_when(
      group_type == "case_control" & Status %in% c("FEP-SCZ", "FEP-BD","Converter") ~ "Case",
      group_type == "case_control" & Status == "Non_Converter" ~ "Control",
      group_type == "fep_converter" ~ gsub("FEP-SCZ|FEP-BD", "FEP", Status),
      group_type == "scz_bd_converter" ~ gsub("FEP-SCZ", "SCZ", gsub("FEP-BD", "BD", Status)),
      TRUE ~ Status
    ))
  
  # Helper functions
  get_outlier_labels <- function(outlier_value, manifest_df) {
    outlier_values <- strsplit(outlier_value, ",")[[1]]
    labels <- sapply(outlier_values, function(val) {
      if (val %in% manifest_df$Sequencing_number) {
        manifest_df$Status[manifest_df$Sequencing_number == val]
      } else {
        NA
      }
    })
    paste(labels, collapse = ", ")
  }

  check_outlier_label <- function(outlier_value, manifest_df) {
    outlier_values <- strsplit(outlier_value, ",")[[1]]
    labels <- sapply(outlier_values, function(val) {
      if (val %in% manifest_df$Sequencing_number) {
        manifest_df$Status[manifest_df$Sequencing_number == val]
      } else {
        NA
      }
    })
    unique_labels <- unique(labels)
    if (length(unique_labels) > 1) "mixed" 
    else if (length(unique_labels) == 1) unique_labels
    else NA
  }

  clean_na_samples <- function(df) {
    df %>%
      filter(!is.na(sample_label_2) & sample_label_2 != "NA") %>%
      mutate(
        sample_label_2 = gsub("\\bNA\\b", "", sample_label_2),
        sample_label_2 = gsub(",+", ",", sample_label_2),
        sample_label_2 = gsub("^,|,$", "", sample_label_2),
        sample_label_2 = trimws(sample_label_2),
        sample_label_2 = ifelse(sample_label_2 == "", NA, sample_label_2)
      ) %>%
      filter(!is.na(sample_label_2))
  }


  derive_group_label <- function(sample_label_2) {
    if (is.na(sample_label_2)) return(NA)
    labels <- trimws(strsplit(sample_label_2, ",")[[1]])
    allowed <- switch(group_type,
      "case_control" = c("Case", "Control"),
      "fep_converter" = c("FEP", "Converter", "Non_Converter"),
      "scz_bd_converter" = c("SCZ", "BD", "Converter", "Non_Converter")
    )
    if (!all(labels %in% allowed)) return(NA)
    unique_labels <- unique(labels)
    if (length(unique_labels) == 1) unique_labels else "mixed"
  }

  # Process each dataframe
  result_dfs <- lapply(dfs, function(df) {
    # Clean SAMPLES

    
    # Add labels
    df$sample_label_2 <- sapply(df$SAMPLES, get_outlier_labels, manifest)
    df$sample_label <- sapply(df$SAMPLES, check_outlier_label, manifest)
    
    # Apply cleaning and grouping
    df <- clean_samples(df)
    df <- clean_na_samples(df)
    df$sample_label_3 <- sapply(df$sample_label_2, derive_group_label)
    df$count <- sapply(strsplit(df$sample_label_2, ","), length)
    
    df
  })
  
  result_dfs
}


dfs_to_process <- list(
  MPC_only = MPC_only,
  Missense = Missense,
  Missense_canonical = Missense_canonical,
  PTV = PTV,
  PTV_canonic = PTV_canonic
)

# Process data for SCZ/BD/Converter groups
processed_data <- process_variant_data(
  group_type = "scz_bd_converter",
  manifest = manifest_correct,
  dfs = dfs_to_process
)

processed_data <- process_variant_data(
  group_type = "case_control",
  manifest = manifest_correct,
  dfs = dfs_to_process
)

processed_data <- process_variant_data(
  group_type = "fep_converter",
  manifest = manifest_correct,
  dfs = dfs_to_process
)

MPC_only <- processed_data$MPC_only
Missense <- processed_data$Missense
Missense_canonical <- processed_data$Missense_canonical
PTV <- processed_data$PTV
PTV_canonic <- processed_data$PTV_canonic

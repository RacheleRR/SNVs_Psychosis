# Load necessary libraries for data manipulation and analysis
library(readr)
library(data.table)
library(dplyr)
library(gprofiler2)
library(httr)
library(readxl)
library(tidyr)
library(stringr)
library(openxlsx)

# Load my data from CSV files
    SCHEMA <- read_csv("/home/rachele/Documents/old/geneset_confrontation/SCHEMA.csv")
    GWAS_120 <- read.delim("~/GWAS_120.csv")       
    BipEx_Bipolar <- read_csv("/home/rachele/Documents/old/geneset_confrontation/Bipolar_Disorder.csv")
    SFARI <- read_csv("/home/rachele/Documents/old/geneset_confrontation/SFARI.csv")

# !Clean GENESETS
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
PTV_canonic <- `Rubiu-GRCh38.deepvariant.splitted.norm.vep.merged.filtered_gnomad_mpc2.filtered_rarity_stringent_AF.PTV_HC.vcf.gz`
rm('Rubiu-GRCh38.deepvariant.splitted.norm.vep.merged.filtered_gnomad_mpc2.filtered_rarity_stringent_AF.missense_mpc_AM.vcf.gz')
rm('Rubiu-GRCh38.deepvariant.splitted.norm.vep.merged.filtered_gnomad_mpc2.filtered_rarity_stringent_AF.missense_mpc_AM_CANONIC.vcf.gz')
rm('Rubiu-GRCh38.deepvariant.splitted.norm.vep.merged.filtered_gnomad_mpc2.filtered_rarity_stringent_AF.missense_mpc.vcf.gz')
PTV_canonic <- PTV_canonic %>% filter(CANONICAL == "YES")

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


setwd("/home/rachele/SNVs/results_pasteur/group2/genes")
# Separate into case and control
# Load packages
library(dplyr)
library(openxlsx)

# --- PARAMETERS ---------------------------------------------------------

group_column <- "sample_label_3"  # Column to split by (e.g., "sample_label_3")
gene_column <- "SYMBOL"           # Column that stores gene names
dedup_columns <- c("SYMBOL", "CHROM", "REF", "ALT", "POS", "SAMPLES")

# --- INPUT: Replace these with your actual datasets ---------------------

datasets <- list(
  Missense = Missense,
  Missense_canonical = Missense_canonical,
  MPC_only = MPC_only,
  PTV_canonic = PTV_canonic,
  PTV = PTV
)

gene_lists <- list(
  schema_pval = unique(SCHEMA_pVal$Gene),
  bipolar = unique(c(BipEx_Bipolar_p_val_PTV$Gene, BipEx_Bipolar_p_val_Missense$Gene)),
  sfari_non_syndromic = unique(SFARI_non_syndromic_lower_3$Gene),
  gwas = unique(GWAS_120$Gene),
  schema_or = unique(SCHEMA_OR$Gene),
  sfari_syndromic = unique(SFARI_syndromic$Gene),
  brain_ntm = unique(brain_gene_consensus_ntm_consensus_no_pitular$Gene.name),
  brain_filt = unique(brain_gene_consensus_filtered_consensus_no_pitular$Gene.name)
)

# --- FUNCTION: Split dataset by arbitrary group labels ------------------

split_by_group <- function(df, group_col = "sample_label_3") {
  df %>%
    split(.[[group_col]]) %>%
    lapply(function(x) distinct(x, across(all_of(dedup_columns)), .keep_all = TRUE))
}


split_datasets <- lapply(datasets, split_by_group)

# --- FUNCTION: Save unique genes from each group to TSV ------------------

save_unique_genes_all <- function(split_datasets, output_prefix = "") {
  for (data_name in names(split_datasets)) {
    for (group in names(split_datasets[[data_name]])) {
      df <- split_datasets[[data_name]][[group]]
      unique_df <- df %>% distinct(across(all_of(gene_column)), .keep_all = TRUE)
      file_name <- paste0(output_prefix, tolower(group), "_", data_name, ".tsv")
      write.table(unique_df, file = file_name, sep = "\t", row.names = FALSE, quote = FALSE)
    }
  }
}

save_unique_genes_all(split_datasets)

# --- FUNCTION: Filter group datasets by gene lists -----------------------

filter_by_gene_lists <- function(split_datasets, gene_lists, gene_col = "SYMBOL") {
  results <- list()
  for (dataset_name in names(split_datasets)) {
    for (group in names(split_datasets[[dataset_name]])) {
      df <- split_datasets[[dataset_name]][[group]]
      for (list_name in names(gene_lists)) {
        filtered <- df %>%
          filter(!!sym(gene_col) %in% gene_lists[[list_name]]) %>%
          mutate(type = list_name, group = group)
        key <- paste(dataset_name, group, list_name, sep = "_")
        results[[key]] <- filtered
      }
    }
  }
  return(results)
}

filtered_results <- filter_by_gene_lists(split_datasets, gene_lists)

# --- FUNCTION: Combine filtered results by dataset -----------------------

combine_by_dataset <- function(filtered_results, dataset_prefix) {
  keys <- grep(paste0("^", dataset_prefix, "_"), names(filtered_results), value = TRUE)
  bind_rows(filtered_results[keys])
}

combined_filtered <- list()
for (name in names(datasets)) {
  combined_filtered[[paste0("combined_", name)]] <- combine_by_dataset(filtered_results, name)
}

# --- FUNCTION: Write all results to Excel files --------------------------

write_list_to_excel <- function(data_list, filename) {
  write.xlsx(data_list, file = filename)
}


#!write list excel but with summary 
# --- FUNCTION: Write Excel with summary sheet ----------------------------
write_list_to_excel <- function(data_list, filename) {
  library(openxlsx)

  wb <- createWorkbook()

  original_names <- names(data_list)

  # Generate safe sheet names: truncate to 31 chars and make unique
  safe_names <- substr(original_names, 1, 31)
  safe_names <- make.unique(safe_names, sep = "_")

  # Create mapping: original name → safe sheet name
  name_mapping <- data.frame(
    Original_Name = original_names,
    Sheet_Name = safe_names,
    Total_Rows = sapply(data_list, nrow),
    Unique_Genes = sapply(data_list, function(df) n_distinct(df[[gene_column]])),
    stringsAsFactors = FALSE
  )

  # Write summary sheet
  addWorksheet(wb, "Summary")
  writeData(wb, "Summary", name_mapping)

  # Write each dataframe with safe sheet name
  for (i in seq_along(data_list)) {
    addWorksheet(wb, safe_names[i])
    writeData(wb, safe_names[i], data_list[[i]])
  }

  saveWorkbook(wb, filename, overwrite = TRUE)
}


# Flat group datasets
flat_split_list <- unlist(split_datasets, recursive = FALSE)
write_list_to_excel(flat_split_list, "grouped_datasets.xlsx")

# Unique genes per group
filter_unique_genes <- function(df) df %>% distinct(across(all_of(gene_column)), .keep_all = TRUE)
unique_split_list <- lapply(flat_split_list, filter_unique_genes)
write_list_to_excel(unique_split_list, "unique_grouped_datasets.xlsx")

# Combined filtered by gene list
write_list_to_excel(combined_filtered, "combined_filtered_by_gene_lists.xlsx")

# Unique version of combined filtered
unique_combined_filtered <- lapply(combined_filtered, filter_unique_genes)
write_list_to_excel(unique_combined_filtered, "unique_combined_filtered_by_gene_lists.xlsx")






# --- FUNCTION: Optional intersection analysis ----------------------------

# Example: find genes shared between brain_filtered and schema_or/gwas in one group
perform_overlap_analysis <- function(group_name = "case_1", dataset_name = "PTV") {
  dataset <- split_datasets[[dataset_name]][[group_name]]
  gene_set <- unique(dataset[[gene_column]])
  brain_genes <- gene_lists$brain_filt
  schema_or_gwas <- union(gene_lists$schema_pval, gene_lists$gwas)
  intersected <- intersect(intersect(gene_set, brain_genes), schema_or_gwas)
  df <- data.frame(Gene = intersected)
  write.table(df, paste0("overlap_", dataset_name, "_", group_name, ".tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
}

# Run example
perform_overlap_analysis("case_1", "PTV")

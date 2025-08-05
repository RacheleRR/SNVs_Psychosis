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
rm(list = c("SCHEMA", "SCHEMA_OR", "SCHEMA_pVal", "SFARI", "SFARI_non_syndromic_lower_3", "SFARI_syndromic", "GWAS_120", "BipEx_Bipolar", "BipEx_Bipolar_p_val_Missense", "BipEx_Bipolar_p_val_PTV", "convert_col"))
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


#separte into group status 
# Function to split dataframe into case and control
split_case_control <- function(df) {
case_df <- df[df$sample_label_3 == "Case", ]
control_df <- df[df$sample_label_3 == "Control", ]
case_df <- case_df %>% distinct(SYMBOL, CHROM,REF,ALT, POS, ID, .keep_all = TRUE) 
control_df <- control_df %>% distinct(SYMBOL, CHROM,REF,ALT, POS,ID,.keep_all=TRUE)
return(list(case = case_df, control = control_df))
}

# Applying the function to each dataset
split_Missense <- split_case_control(Missense)
split_Missense_canonical <- split_case_control(Missense_canonical)
split_MPC_only <- split_case_control(MPC_only)
split_PTV_canonic <- split_case_control(PTV_canonic)
split_PTV <- split_case_control(PTV)

# Accessing the case and control dataframes
case_Missense <- split_Missense$case
control_Missense <- split_Missense$control
case_Missense_canonical <- split_Missense_canonical$case
control_Missense_canonical <- split_Missense_canonical$control
case_MPC_only <- split_MPC_only$case
control_MPC_only <- split_MPC_only$control
case_PTV_canonic <- split_PTV_canonic$case
control_PTV_canonic <- split_PTV_canonic$control
case_PTV <- split_PTV$case
control_PTV <- split_PTV$control



# set working directory 
 setwd("/home/rachele/SNVs/results_pasteur/genes")
# Function to filter unique genes and save as TSV
save_unique_genes <- function(df, filename) {
unique_df <- df[!duplicated(df$SYMBOL), ]  # Remove duplicates based on GENE column
write.table(unique_df, file = filename, sep = "\t", row.names = FALSE, quote = FALSE)
}

# Apply to case and control datasets
save_unique_genes(case_Missense, "case_Missense.tsv")
save_unique_genes(control_Missense, "control_Missense.tsv")
save_unique_genes(case_Missense_canonical, "case_Missense_canonical.tsv")
save_unique_genes(control_Missense_canonical, "control_Missense_canonical.tsv")
save_unique_genes(case_MPC_only, "case_MPC_only.tsv")
save_unique_genes(control_MPC_only, "control_MPC_only.tsv")
save_unique_genes(case_PTV_canonic, "case_PTV_canonic.tsv")
save_unique_genes(control_PTV_canonic, "control_PTV_canonic.tsv")
save_unique_genes(case_PTV, "case_PTV.tsv")
save_unique_genes(control_PTV, "control_PTV.tsv")


# GENES IN SCHIZO ETC 
# Define the unique gene lists
genes_schema_pval <- unique(SCHEMA_pVal$Gene)
genes_bipolar <- unique(c(BipEx_Bipolar_p_val_PTV$Gene, BipEx_Bipolar_p_val_Missense$Gene))
genes_sfari_non_syndromic <- unique(SFARI_non_syndromic_lower_3$Gene)
genes_gwas <- unique(GWAS_120$Gene)
genes_schema_or <- unique(SCHEMA_OR$Gene)
genes_sfari_syndromic <- unique(SFARI_syndromic$Gene)
genes_brain_ntm <- unique(brain_gene_consensus_ntm_consensus_no_pitular$Gene.name)
genes_brain_filt <- unique(brain_gene_consensus_filtered_consensus_no_pitular$Gene.name)

# Function to filter genes that are present in both the case dataset and a given gene list
filter_genes <- function(case_df, gene_list, type) {
  filtered_df <- case_df %>% filter(SYMBOL %in% gene_list) %>% mutate(type = type)
  return(filtered_df)
}

# Filter the case datasets based on the gene lists
filtered_case_Missense_canonical_schema_pval <- filter_genes(case_Missense_canonical, genes_schema_pval, "schema_pval")
filtered_case_Missense_canonical_bipolar <- filter_genes(case_Missense_canonical, genes_bipolar, "bipolar")
filtered_case_Missense_canonical_sfari_non_syndromic <- filter_genes(case_Missense_canonical, genes_sfari_non_syndromic, "sfari_non_syndromic")
filtered_case_Missense_canonical_gwas <- filter_genes(case_Missense_canonical, genes_gwas, "gwas")
filtered_case_Missense_canonical_schema_or <- filter_genes(case_Missense_canonical, genes_schema_or, "schema_or")
filtered_case_Missense_canonical_sfari_syndromic <- filter_genes(case_Missense_canonical, genes_sfari_syndromic, "sfari_syndromic")
filtered_case_Missense_canonical_brain_ntm <- filter_genes(case_Missense_canonical, genes_brain_ntm, "brain_ntm")
filtered_case_Missense_canonical_brain_filt <- filter_genes(case_Missense_canonical, genes_brain_filt, "brain_filt")


filtered_case_Missense_schema_pval <- filter_genes(case_Missense, genes_schema_pval, "schema_pval")
filtered_case_Missense_bipolar <- filter_genes(case_Missense, genes_bipolar, "bipolar")
filtered_case_Missense_sfari_non_syndromic <- filter_genes(case_Missense, genes_sfari_non_syndromic, "sfari_non_syndromic")
filtered_case_Missense_gwas <- filter_genes(case_Missense, genes_gwas, "gwas")
filtered_case_Missense_schema_or <- filter_genes(case_Missense, genes_schema_or, "schema_or")
filtered_case_Missense_sfari_syndromic <- filter_genes(case_Missense, genes_sfari_syndromic, "sfari_syndromic")
filtered_case_Missense_brain_ntm <- filter_genes(case_Missense, genes_brain_ntm, "brain_ntm")
filtered_case_Missense_brain_filt <- filter_genes(case_Missense, genes_brain_filt, "brain_filt")



filtered_case_PTV_schema_pval <- filter_genes(case_PTV, genes_schema_pval, "schema_pval")
filtered_case_PTV_bipolar <- filter_genes(case_PTV, genes_bipolar, "bipolar")
filtered_case_PTV_sfari_non_syndromic <- filter_genes(case_PTV, genes_sfari_non_syndromic, "sfari_non_syndromic")
filtered_case_PTV_gwas <- filter_genes(case_PTV, genes_gwas, "gwas")
filtered_case_PTV_schema_or <- filter_genes(case_PTV, genes_schema_or, "schema_or")
filtered_case_PTV_sfari_syndromic <- filter_genes(case_PTV, genes_sfari_syndromic, "sfari_syndromic")
filtered_case_PTV_brain_ntm <- filter_genes(case_PTV, genes_brain_ntm, "brain_ntm")
filtered_case_PTV_brain_filt <- filter_genes(case_PTV, genes_brain_filt, "brain_filt")

filtered_case_PTV_canonic_schema_pval <- filter_genes(case_PTV_canonic, genes_schema_pval, "schema_pval")
filtered_case_PTV_canonic_bipolar <- filter_genes(case_PTV_canonic, genes_bipolar, "bipolar")
filtered_case_PTV_canonic_sfari_non_syndromic <- filter_genes(case_PTV_canonic, genes_sfari_non_syndromic, "sfari_non_syndromic")
filtered_case_PTV_canonic_gwas <- filter_genes(case_PTV_canonic, genes_gwas, "gwas")
filtered_case_PTV_canonic_schema_or <- filter_genes(case_PTV_canonic, genes_schema_or, "schema_or")
filtered_case_PTV_canonic_sfari_syndromic <- filter_genes(case_PTV_canonic, genes_sfari_syndromic, "sfari_syndromic")
filtered_case_PTV_canonic_brain_ntm <- filter_genes(case_PTV_canonic, genes_brain_ntm, "brain_ntm")
filtered_case_PTV_canonic_brain_filt <- filter_genes(case_PTV_canonic, genes_brain_filt, "brain_filt")

filtered_case_MPC_only_schema_pval <- filter_genes(case_MPC_only, genes_schema_pval, "schema_pval")
filtered_case_MPC_only_bipolar <- filter_genes(case_MPC_only, genes_bipolar, "bipolar")
filtered_case_MPC_only_sfari_non_syndromic <- filter_genes(case_MPC_only, genes_sfari_non_syndromic, "sfari_non_syndromic")
filtered_case_MPC_only_gwas <- filter_genes(case_MPC_only, genes_gwas, "gwas")
filtered_case_MPC_only_schema_or <- filter_genes(case_MPC_only, genes_schema_or, "schema_or")
filtered_case_MPC_only_sfari_syndromic <- filter_genes(case_MPC_only, genes_sfari_syndromic, "sfari_syndromic")
filtered_case_MPC_only_brain_ntm <- filter_genes(case_MPC_only, genes_brain_ntm, "brain_ntm")
filtered_case_MPC_only_brain_filt <- filter_genes(case_MPC_only, genes_brain_filt, "brain_filt")

# Combine the filtered dataframes into a single dataframe for each case dataset
combined_case_Missense_canonical <- bind_rows(
  filtered_case_Missense_canonical_schema_pval,
  filtered_case_Missense_canonical_bipolar,
  filtered_case_Missense_canonical_sfari_non_syndromic,
  filtered_case_Missense_canonical_gwas,
  filtered_case_Missense_canonical_schema_or,
  filtered_case_Missense_canonical_sfari_syndromic,
  filtered_case_Missense_canonical_brain_ntm,
  filtered_case_Missense_canonical_brain_filt
)

combined_case_Missense <- bind_rows(
  filt_case_Missense_schema_pval,
  filt_case_Missense_bipolar,
  filt_case_Missense_sfari_non_syndromic,
  filt_case_Missense_gwas,
  filt_case_Missense_schema_or,
  filt_case_Missense_sfari_syndromic,
  filt_case_Missense_brain_ntm,
  filt_case_Missense_brain_filt
)


combined_case_PTV <- bind_rows(
  filtered_case_PTV_schema_pval,
  filtered_case_PTV_bipolar,
  filtered_case_PTV_sfari_non_syndromic,
  filtered_case_PTV_gwas,
  filtered_case_PTV_schema_or,
  filtered_case_PTV_sfari_syndromic,
  filtered_case_PTV_brain_ntm,
  filtered_case_PTV_brain_filt
)

combined_case_PTV_canonic <- bind_rows(
  filtered_case_PTV_canonic_schema_pval,
  filtered_case_PTV_canonic_bipolar,
  filtered_case_PTV_canonic_sfari_non_syndromic,
  filtered_case_PTV_canonic_gwas,
  filtered_case_PTV_canonic_schema_or,
  filtered_case_PTV_canonic_sfari_syndromic,
  filtered_case_PTV_canonic_brain_ntm,
  filtered_case_PTV_canonic_brain_filt
)

combined_case_MPC_only <- bind_rows(
  filtered_case_MPC_only_schema_pval,
  filtered_case_MPC_only_bipolar,
  filtered_case_MPC_only_sfari_non_syndromic,
  filtered_case_MPC_only_gwas,
  filtered_case_MPC_only_schema_or,
  filtered_case_MPC_only_sfari_syndromic,
  filtered_case_MPC_only_brain_ntm,
  filtered_case_MPC_only_brain_filt
)

# Save the combined dataframes to new TSV files
write.table(combined_case_Missense_canonical, "combined_case_Missense_canonical.tsv", sep="\t", row.names=FALSE, quote=FALSE)
write.table(combined_case_Missense, "combined_case_Missense.tsv", sep="\t", row.names=FALSE, quote=FALSE)
write.table(combined_case_PTV, "combined_case_PTV.tsv", sep="\t", row.names=FALSE, quote=FALSE)
write.table(combined_case_PTV_canonic, "combined_case_PTV_canonic.tsv", sep="\t", row.names=FALSE, quote=FALSE)
write.table(combined_case_MPC_only, "combined_case_MPC_only.tsv", sep="\t", row.names=FALSE, quote=FALSE)


# FIND GENES WITHIN MULTIPLES
# Define the unique gene lists
genes_schema_pval <- unique(SCHEMA_pVal$Gene)
genes_gwas <- unique(GWAS_120$Gene)
genes_brain_filt <- unique(brain_gene_consensus_filtered_consensus_no_pitular$Gene)
genes_case_PTV <- unique(case_PTV$SYMBOL)

# Find the genes that are in case_PTV and brain_filt
genes_case_PTV_brain <- intersect(genes_case_PTV, genes_brain_filt)

# Find the genes that are either in schema_pval or gwas
genes_schema_or_gwas <- union(genes_schema_pval, genes_gwas)

# Find the genes that are in case_PTV_brain and either in schema_pval or gwas
common_genes <- intersect(genes_case_PTV_brain, genes_schema_or_gwas)

# Create a dataframe with the unique genes
common_genes_df <- data.frame(Gene = common_genes)

# Save the dataframe to a TSV file
write.table(common_genes_df, "common_genes_case_PTV_brain_schema_or_gwas.tsv", sep="\t", row.names=FALSE, quote=FALSE)



# !Excel document 
# case_control genes
  #Create a list of all case and control datasets
  case_control_list <- list(
    case_Missense = case_Missense,
    control_Missense = control_Missense,
    case_Missense_canonical = case_Missense_canonical,
    control_Missense_canonical = control_Missense_canonical,
    case_MPC_only = case_MPC_only,
    control_MPC_only = control_MPC_only,
    case_PTV_canonic = case_PTV_canonic,
    control_PTV_canonic = control_PTV_canonic,
    case_PTV = case_PTV,
    control_PTV = control_PTV
  )

  #Save the list to an Excel file
  write.xlsx(case_control_list, file = "case_control_datasets.xlsx")

# case control genes unqiue 
  #Create a function to filter for unique genes based on the SYMBOL column
  filter_unique_genes <- function(df) {
    df %>% distinct(SYMBOL, .keep_all = TRUE)
  }

  #Apply the function to each dataset in the list
  unique_case_control_list <- lapply(case_control_list, filter_unique_genes)

  #Save the unique datasets to another Excel file
  write.xlsx(unique_case_control_list, file = "unique_case_control_datasets.xlsx")

#  case control with disaes and brain 
#Create a list of all combined and filtered datasets
#Combine the filtered dataframes into a single dataframe for each case dataset
combined_case_Missense_canonical <- bind_rows(
  filtered_case_Missense_canonical_schema_pval,
  filtered_case_Missense_canonical_bipolar,
  filtered_case_Missense_canonical_sfari_non_syndromic,
  filtered_case_Missense_canonical_gwas,
  filtered_case_Missense_canonical_schema_or,
  filtered_case_Missense_canonical_sfari_syndromic
)

combined_case_Missense <- bind_rows(
  filt_case_Missense_schema_pval,
  filt_case_Missense_bipolar,
  filt_case_Missense_sfari_non_syndromic,
  filt_case_Missense_gwas,
  filt_case_Missense_schema_or,
  filt_case_Missense_sfari_syndromic
)

combined_case_PTV <- bind_rows(
  filtered_case_PTV_schema_pval,
  filtered_case_PTV_bipolar,
  filtered_case_PTV_sfari_non_syndromic,
  filtered_case_PTV_gwas,
  filtered_case_PTV_schema_or,
  filtered_case_PTV_sfari_syndromic
)

combined_case_PTV_canonic <- bind_rows(
  filtered_case_PTV_canonic_schema_pval,
  filtered_case_PTV_canonic_bipolar,
  filtered_case_PTV_canonic_sfari_non_syndromic,
  filtered_case_PTV_canonic_gwas,
  filtered_case_PTV_canonic_schema_or,
  filtered_case_PTV_canonic_sfari_syndromic
)

combined_case_MPC_only <- bind_rows(
  filtered_case_MPC_only_schema_pval,
  filtered_case_MPC_only_bipolar,
  filtered_case_MPC_only_sfari_non_syndromic,
  filtered_case_MPC_only_gwas,
  filtered_case_MPC_only_schema_or,
  filtered_case_MPC_only_sfari_syndromic
)

combined_filtered_list <- list(
  combined_case_Missense_canonical = combined_case_Missense_canonical,
  filt_case_Missense_canonical_brainntm = filtered_case_Missense_canonical_brain_ntm,
  filt_case_Missense_canonical_brainfil = filtered_case_Missense_canonical_brain_filt,
  combined_case_PTV = combined_case_PTV,
  combined_case_PTV_canonic = combined_case_PTV_canonic,
  combined_case_MPC_only = combined_case_MPC_only,
  combined_case_Missense = combined_case_Missense,
  filt_case_Missense_brainntm = filtered_case_Missense_brain_ntm,
  filt_case_Missense_brainfil = filtered_case_Missense_brain_filt,
  filt_case_PTV_brain_ntm = filtered_case_PTV_brain_ntm,
  filt_case_PTV_brain_filt = filtered_case_PTV_brain_filt,
  filt_case_PTV_canonic_brain_ntm = filtered_case_PTV_canonic_brain_ntm,
  filt_case_PTV_canonic_brain_filt = filtered_case_PTV_canonic_brain_filt,
  filt_case_MPC_only_brainntm = filtered_case_MPC_only_brain_ntm,
  filt_case_MPC_only_brainfil = filtered_case_MPC_only_brain_filt
)

#Save the list to an Excel file
write.xlsx(combined_filtered_list, file = "combined_filtered_datasets.xlsx")

# case control with disease and brain unique 
#Apply the unique filter function to each dataset in the combined_filtered_list
unique_combined_filtered_list <- lapply(combined_filtered_list, filter_unique_genes)

#Save the unique datasets to another Excel file
write.xlsx(unique_combined_filtered_list, file = "unique_combined_filtered_datasets.xlsx")


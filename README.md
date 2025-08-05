---

# 🧬 Welcome to SNVs\_Psychosis!

This project guides you through processing, annotating, and statistically analyzing Single Nucleotide Variants (SNVs) using a structured pipeline. The workflow covers filtering, transformation, annotation, and statistical testing — leading up to biological interpretation through enrichment and network analysis.

---

## 🌟 Overview

This repository walks you through how to:

* ✅ Filter variants based on biological scores and annotations
* ✅ Convert VCF files into TSVs for R-friendly processing
* ✅ Annotate SNVs using public databases
* ✅ Run statistical, enrichment, and network analyses

---

## 🎯 Objectives

* ✔️ Apply variant filtering with score and canonic thresholds
* ✔️ Convert VCF files to TSV format for visualization and analysis
* ✔️ Annotate SNVs using external databases
* ✔️ Run statistical tests for group comparisons
* ✔️ Perform enrichment and network-based biological interpretation

---

## 🛠️ Prerequisites

Ensure you have the following tools installed:

* `bcftools` 🔧
* `R` with required packages (e.g., `data.table`, `ggplot2`, `stats`, etc.) 📊

---

## 🔍 SNV Processing Steps

### 1️⃣ Filter Variants

Use `filters_snvs_pasteur.sh` to retain only relevant variants:

* ✅ **MPC only**: MPC ≥ 2.6
* ✅ **Missense**: MPC ≥ 2.6 and AM ≥ 0.545
* ✅ **Missense Canonic**: same as Missense + canonic = YES
* ✅ **PTV Canonic**: high-confidence LoF + canonic = YES
* 🧬 **PTVs**:  high-confidence LoF 

---

### 2️⃣ Convert VCF to TSV

Use `transform_VCF.sh` to make data R-friendly.

**Example usage:**

```bash
./transform_VCF.sh -i path/to/your_file.vcf.gz
```

**Example file:**

```bash
./transform_VCF.sh -i /media/rachele/One\ Touch/Rubiu-GRCh38.deepvariant.splitted.norm.vep.merged.filtered_gnomad_mpc2.filtered_rarity_stringent_AF.PTV_HC.vcf.gz
```

---

### 3️⃣ Statistical Tests

#### For 2-Group Comparison:

* `FISHER_SNVs_PASTEUR_2GROUPS.r`
* `Wilcoxon_SNVs_PASTEUR_2GROUPS.r`
* `Regression_and_poisson_2_groups.r`

#### For 3-Group Comparison:

* `FISHER_SNVs_PASTEUR_3GROUPS.r`
* `Kruskal_SNVs_PASTEUR_3GROUPS.r`
* `Regression_and_poisson_3_groups.r`

#### For 4-Group Comparison:

* `FISHER_SNVs_PASTEUR_4GROUPS.r`
* `Kruskal_SNVs_PASTEUR_4GROUPS.r`
* `Regression_and_poisson_4_groups.r`

---

### 4️⃣ Gene Lists

Extract gene names based on filters, variant categories, and public databases.

Run:

```bash
Rscript Create_gene_list.r
```

---

### 5️⃣ Enrichment Analysis

Use tools like:

* [Enrichr](https://maayanlab.cloud/Enrichr/)
* Cytoscape + EnrichmentMap
* GSEA

To identify enriched biological processes, pathways, and gene ontologies.

---

### 6️⃣ Network Analysis

Visualize gene relationships using:

* [GeneMANIA](https://genemania.org/)
* [STRING](https://string-db.org/)
* [HumanBase](https://hb.flatironinstitute.org/)

---

## 💻 Command Summary

| Task                         | Command / Script                                            |
| ---------------------------- | ----------------------------------------------------------- |
| Filter SNVs                  | `./filters_snvs_pasteur.sh`                                 |
| Convert VCF to TSV           | `./transform_VCF.sh -i your_file.vcf.gz`                    |
| Statistical Test (2 Groups)  | `FISHER_SNVs_PASTEUR_2GROUPS.r`, etc.                       |
| Statistical Test (3 Groups)  | `FISHER_SNVs_PASTEUR_3GROUPS.r`, etc.                       |
| Statistical Test (4 Groups)  | `FISHER_SNVs_PASTEUR_4GROUPS.r`, etc.                       |
| Preprocessing for Enrichment | `Rscript preprocessing_step_function_for_all_group_types.r` |
| Create Gene List             | `Rscript Create_gene_list.r`                                |

---


## 🧾 License

MIT License — free to use and modify. Please cite if used in publications.

---


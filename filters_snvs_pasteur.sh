
#! FILTER for allele frequency

bcftools view -e 'INFO/AF_genomes > 0.01' \
  -Oz \
  -o "/media/rachele/One Touch/Rubiu-GRCh38.deepvariant.splitted.norm.vep.merged.filtered_gnomad_mpc2.filtered_rarity_stringent_AF.vcf.gz" \
  Rubiu-GRCh38.deepvariant.splitted.norm.vep.merged.filtered_gnomad_mpc2.vcf.gz


#! PATHOGENISTY
bcftools view -i 'MPC >= 2.6' Rubiu-GRCh38.deepvariant.splitted.norm.vep.merged.filtered_gnomad_mpc2.filtered_rarity_stringent_AF.vcf.gz -Oz -o Rubiu-GRCh38.deepvariant.splitted.norm.vep.merged.filtered_gnomad_mpc2.filtered_rarity_stringent_AF.missense_mpc.vcf.gz

filter_vep \
  -i  Rubiu-GRCh38.deepvariant.splitted.norm.vep.merged.filtered_gnomad_mpc2.filtered_rarity_stringent_AF.missense_mpc.vcf.gz \
  -o  Rubiu-GRCh38.deepvariant.splitted.norm.vep.merged.filtered_gnomad_mpc2.filtered_rarity_stringent_AF.missense_mpc_AM_98.vcf.gz \
  --only_matched \
  --filter "am_pathogenicity >= 0.98"

filter_vep \
  -i  Rubiu-GRCh38.deepvariant.splitted.norm.vep.merged.filtered_gnomad_mpc2.filtered_rarity_stringent_AF.missense_mpc.vcf.gz \
  -o  Rubiu-GRCh38.deepvariant.splitted.norm.vep.merged.filtered_gnomad_mpc2.filtered_rarity_stringent_AF.missense_mpc_AM.vcf.gz \
  --only_matched \
  --filter "am_pathogenicity >= 0.565"

filter_vep \
  -i  Rubiu-GRCh38.deepvariant.splitted.norm.vep.merged.filtered_gnomad_mpc2.filtered_rarity_stringent_AF.missense_mpc.vcf.gz \
  -o  Rubiu-GRCh38.deepvariant.splitted.norm.vep.merged.filtered_gnomad_mpc2.filtered_rarity_stringent_AF.missense_mpc_AM_CANONIC.vcf.gz \
  --only_matched \
  --filter "am_pathogenicity >= 0.565 and CANONICAL is YES"

# PTVs 
filter_vep \
  -i Rubiu-GRCh38.deepvariant.splitted.norm.vep.merged.filtered_gnomad_mpc2.filtered_rarity_stringent_AF.vcf.gz  \
  -o Rubiu-GRCh38.deepvariant.splitted.norm.vep.merged.filtered_gnomad_mpc2.filtered_rarity_stringent_AF.PTV_HC.vcf.gz  \
  --only_matched  \
  --filter "LoF is HC"

filter_vep \
  -i Rubiu-GRCh38.deepvariant.splitted.norm.vep.merged.filtered_gnomad_mpc2.filtered_rarity_stringent_AF.vcf.gz \
  -o Rubiu-GRCh38.deepvariant.splitted.norm.vep.merged.filtered_gnomad_mpc2.filtered_rarity_stringent_AF.PTV_HC3.vcf.gz \
  --only_matched \
  --filter "LoF is HC and CANONICAL is YES"



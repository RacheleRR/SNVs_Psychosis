
bcftools view -H Rubiu-GRCh38.deepvariant.splitted.norm.vep.merged.filtered_gnomad_mpc2.vcf.gz | wc -l

#! FILTER 1

bcftools view -e 'INFO/AF_genomes > 0.01' \
  -Oz \
  -o "/media/rachele/One Touch/Rubiu-GRCh38.deepvariant.splitted.norm.vep.merged.filtered_gnomad_mpc2.filtered_rarity_stringent_AF.vcf.gz" \
  Rubiu-GRCh38.deepvariant.splitted.norm.vep.merged.filtered_gnomad_mpc2.vcf.gz

bcftools view -H Rubiu-GRCh38.deepvariant.splitted.norm.vep.merged.filtered_gnomad_mpc2.filtered_rarity_stringent_AF.vcf.gz | wc -l



# PATHOGENISTY
bcftools view -i 'MPC >= 2.6' Rubiu-GRCh38.deepvariant.splitted.norm.vep.merged.filtered_gnomad_mpc2.filtered_rarity_stringent_AF.vcf.gz -Oz -o Rubiu-GRCh38.deepvariant.splitted.norm.vep.merged.filtered_gnomad_mpc2.filtered_rarity_stringent_AF.missense_mpc.vcf.gz

bcftools view -H Rubiu-GRCh38.deepvariant.splitted.norm.vep.merged.filtered_gnomad_mpc2.filtered_rarity_stringent_AF.missense_mpc.vcf.gz | wc -l


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

bcftools view -H Rubiu-GRCh38.deepvariant.splitted.norm.vep.merged.filtered_gnomad_mpc2.filtered_rarity_stringent_AF.missense_mpc_AM.vcf.gz | wc -l

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

filter_vep \
  -i Rubiu-GRCh38.deepvariant.splitted.norm.vep.merged.filtered_gnomad_mpc2.filtered_rarity_stringent_AF.vcf.gz \
  -o Rubiu-GRCh38.deepvariant.splitted.norm.vep.merged.filtered_gnomad_mpc2.filtered_rarity_stringent_AF.PTV_HC6.vcf.gz \
  --only_matched \
  --filter "LoF is HC and CANONICAL is YES" 

bcftools view -H Rubiu-GRCh38.deepvariant.splitted.norm.vep.merged.filtered_gnomad_mpc2.filtered_rarity_stringent_AF.PTV_HC.vcf.gz | wc -l


#OPTIONALS 
# PTVs with LoF not defined
filter_vep \
  -i Rubiu-GRCh38.deepvariant.splitted.norm.vep.merged.filtered_gnomad_mpc2.filtered_rarity_stringent_AF.vcf.gz  \
  -o Rubiu-GRCh38.deepvariant.splitted.norm.vep.merged.filtered_gnomad_mpc2.filtered_rarity_stringent_AF.PTV.vcf.gz  \
  --only_matched  \
  --filter "(
    Consequence is 'frameshift_variant' or
    Consequence is 'stop_gained' or
    Consequence is 'splice_acceptor_variant' or
    Consequence is 'splice_donor_variant' or
    Consequence is 'start_lost'
  ) and (LoF is 'HC' or not LoF)"

screen -r vep_filter


#
bcftools view -e 'INFO/nhomalt_joint != 0 || INFO/AC_genomes >= 10 || INFO/genomes_filters != "PASS"' \
  -Oz \
  -o "/media/rachele/One Touch/Rubiu-GRCh38.deepvariant.splitted.norm.vep.merged.filtered_gnomad_mpc2.filtered_rarity_stringent_INFO.vcf.gz"\
  "/media/rachele/One Touch/Rubiu-GRCh38.deepvariant.splitted.norm.vep.merged.filtered_gnomad_mpc2.filtered_rarity_stringent_AF.vcf.gz"

bcftools view -e 'INFO/nhomalt_joint <= 1  || INFO/AC_genomes >= 50|| INFO/genomes_filters != "PASS"' \
  -Oz \
  -o "/media/rachele/One Touch/Rubiu-GRCh38.deepvariant.splitted.norm.vep.merged.filtered_gnomad_mpc2.filtered_rarity_leniant_INFO.vcf.gz"\
  "/media/rachele/One Touch/Rubiu-GRCh38.deepvariant.splitted.norm.vep.merged.filtered_gnomad_mpc2.filtered_rarity_stringent_AF.vcf.gz"

bcftools view -e 'INFO/nhomalt_joint != 0 ' \
  -Oz \
  -o "/media/rachele/One Touch/Rubiu-GRCh38.deepvariant.splitted.norm.vep.merged.filtered_gnomad_mpc2.filtered_rarity_stringent_nH.vcf.gz"\
  "/media/rachele/One Touch/Rubiu-GRCh38.deepvariant.splitted.norm.vep.merged.filtered_gnomad_mpc2.filtered_rarity_stringent_AF.vcf.gz"

bcftools view -e 'INFO/nhomalt_joint <= 1 ' \
  -Oz \
  -o "/media/rachele/One Touch/Rubiu-GRCh38.deepvariant.splitted.norm.vep.merged.filtered_gnomad_mpc2.filtered_rarity_leniant_nH.vcf.gz"\
  "/media/rachele/One Touch/Rubiu-GRCh38.deepvariant.splitted.norm.vep.merged.filtered_gnomad_mpc2.filtered_rarity_stringent_AF.vcf.gz"

bcftools view -e 'INFO/AC_genomes >= 10 ' \
  -Oz \
  -o "/media/rachele/One Touch/Rubiu-GRCh38.deepvariant.splitted.norm.vep.merged.filtered_gnomad_mpc2.filtered_rarity_stringent_AC.vcf.gz"\
  "/media/rachele/One Touch/Rubiu-GRCh38.deepvariant.splitted.norm.vep.merged.filtered_gnomad_mpc2.filtered_rarity_stringent_AF.vcf.gz"

bcftools view -e 'INFO/AC_genomes >= 50 ' \
  -Oz \
  -o "/media/rachele/One Touch/Rubiu-GRCh38.deepvariant.splitted.norm.vep.merged.filtered_gnomad_mpc2.filtered_rarity_leniant_AC.vcf.gz"\
  "/media/rachele/One Touch/Rubiu-GRCh38.deepvariant.splitted.norm.vep.merged.filtered_gnomad_mpc2.filtered_rarity_stringent_AF.vcf.gz"


bcftools view -e 'INFO/genomes_filters != "PASS"' \
  -Oz \
  -o "/media/rachele/One Touch/Rubiu-GRCh38.deepvariant.splitted.norm.vep.merged.filtered_gnomad_mpc2.filtered_rarity_stringent_PASS.vcf.gz"\
  "/media/rachele/One Touch/Rubiu-GRCh38.deepvariant.splitted.norm.vep.merged.filtered_gnomad_mpc2.filtered_rarity_stringent_AF.vcf.gz"




#code to see all variant amount 
#!/bin/bash

echo "Filename,Variant_Count"

for file in *.vcf.gz; do
    count=$(bcftools view -H "$file" | wc -l)
    echo "$file,$count"
done


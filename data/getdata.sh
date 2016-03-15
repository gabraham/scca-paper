#!/bin/bash

set -e

# See ftp://ftp.ncbi.nlm.nih.gov/hapmap/genotypes/2009-01_phaseIII/00README.txt
# for information about the genotyping QC procedures

# Download the genotype data (not the consensus data):
#wget ftp://ftp.ncbi.nlm.nih.gov/hapmap/genotypes/2010-05_phaseIII/plink_format/hapmap3_r3_b36_fwd.qc.poly.tar.gz
#wget ftp://ftp.ncbi.nlm.nih.gov/hapmap/genotypes/2010-05_phaseIII/relationships_w_pops_041510.txt

#tar xzvf hapmap3_r3_b36_fwd.qc.poly.tar.gz

# Typo in directory name from HapMap?
DIR=hapmap3_r2_b36_fwd.qc.poly

#pushd $DIR
#for f in *.map
#do
#   plink --file ${f/.map/} \
#      --make-bed \
#      --out ${f/.map}
#
#   plink --bfile ${f/.map} \
#      --maf 0.01 \
#      --geno 0.1 \
#      --mind 0.1 \
#      --hwe 1e-6 \
#      --autosome \
#      --filter-founders \
#      --make-bed \
#      --out ${f/.map/.filtered}
#done
#
## Only look at the 8 populations with gene expression data
#awk '{print $2}' \
#   hapmap3_r3_b36_fwd.{CEU,CHB,GIH,JPT,LWK,MEX,MKK,YRI}.qc.poly.bim \
#   | LC_ALL=C sort | uniq -c | awk '$1 == 8 {print $2}' \
#   > orig_consensus_snps.txt
#
#awk '{print $2}' \
#   hapmap3_r3_b36_fwd.{CEU,CHB,GIH,JPT,LWK,MEX,MKK,YRI}.qc.poly.filtered.bim \
#   | LC_ALL=C sort | uniq -c | awk '$1 == 8 {print $2}' \
#   > extrafiltering_consensus_snps.txt
#
#rm -f merge_list.txt
#for p in CEU CHB GIH JPT LWK MEX MKK YRI
#do
#   echo hapmap3_r3_b36_fwd.${p}.qc.poly.filtered.{bed,bim,fam} \
#      >> merge_list.txt
#done
#
#plink --merge-list merge_list.txt \
#   --make-bed \
#   --out tmp
#
#plink --bfile tmp \
#   --extract extrafiltering_consensus_snps.txt \
#   --make-bed \
#   --out hapmap3_r3_b36_fwd.qc.poly.filtered
#popd
#
## Download the gene expression data
#wget http://www.ebi.ac.uk/arrayexpress/files/E-MTAB-264/E-MTAB-264.processed.1.zip
#wget http://www.ebi.ac.uk/arrayexpress/files/E-MTAB-198/E-MTAB-198.processed.1.zip
#
#unzip E-MTAB-264.processed.1.zip
#unzip E-MTAB-198.processed.1.zip
#
## Rename CEU for consistency (ignore CEU RNA-seq data)
#cp normalized_array_data.109.txt CEU_p3_expression.txt
#
## Get the sample IDs for the gene expression data
#grep Normalization *_p3_expression.txt \
#   | cut -f 2- -d$'\t' | tr '\t' '\n' \
#   > expression_samples.txt
#
#wc -l expression_samples.txt
#
#awk 'NR == FNR {a[$1]; next} $2 in a {print $1, $2}' \
#   expression_samples.txt \
#   $DIR/hapmap3_r3_b36_fwd.qc.poly.filtered.fam \
#   > common_samples.txt
#
#plink \
#   --bfile $DIR/hapmap3_r3_b36_fwd.qc.poly.filtered \
#   --keep common_samples.txt \
#   --make-bed \
#   --out $DIR/hapmap3_r3_b36_fwd.qc.poly.filtered.wexpr
#
#plink \
#   --bfile $DIR/hapmap3_r3_b36_fwd.qc.poly.filtered.wexpr \
#   --chr 1 \
#   --make-bed \
#   --out $DIR/hapmap3_r3_b36_fwd.qc.poly.filtered.wexpr.chr1

# Also extract the samples with gene expression, for each population
for p in CEU CHB GIH JPT LWK MEX MKK YRI
do
   plink \
      --bfile $DIR/hapmap3_r3_b36_fwd.${p}.qc.poly.filtered \
      --keep common_samples.txt \
      --make-bed \
      --out $DIR/hapmap3_r3_b36_fwd.${p}.qc.poly.filtered.wexpr
done


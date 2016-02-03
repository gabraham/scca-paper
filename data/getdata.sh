#!/bin/bash

# Download the genotype data:
wget ftp://ftp.ncbi.nlm.nih.gov/hapmap/genotypes/2009-01_phaseIII/plink_format/hapmap3_r2_b36_fwd.consensus.qc.poly.ped.bz2
wget ftp://ftp.ncbi.nlm.nih.gov/hapmap/genotypes/2009-01_phaseIII/plink_format/hapmap3_r2_b36_fwd.consensus.qc.poly.map.bz2

bzip2 -d hapmap3_r2_b36_fwd.consensus.qc.poly.ped.bz2
bzip2 -d hapmap3_r2_b36_fwd.consensus.qc.poly.map.bz2

plink \
   --file hapmap3_r2_b36_fwd.consensus.qc.poly \
   --make-bed \
   --maf 0.05 \
   --hwe 5e-6 \
   --mind 0.01 \
   --geno 0.01 \
   --filter-founders \
   --out hapmap3_r2_b36_fwd.consensus.qc.poly_founders_filtered

plink \
   --bfile hapmap3_r2_b36_fwd.consensus.qc.poly_founders_filtered \
   --chr 1 \
   --make-bed \
   --out hapmap3_r2_b36_fwd.consensus.qc.poly_founders_filtered


rm hapmap3_r2_b36_fwd.consensus.qc.poly.{ped,map}

# Download the gene expression data from
wget http://www.ebi.ac.uk/arrayexpress/files/E-MTAB-264/E-MTAB-264.processed.1.zip
unzip E-MTAB-264.processed.1.zip


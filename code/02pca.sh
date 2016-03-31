#!/bin/bash

set -e

PLINK=plink
FPCA=flashpca

type $PLINK >/dev/null 2>&1 || {
   echo >&2 "plink is not available; set PLINK variable in 02pca.sh script or fix path"
   exit 1
}
type $FPCA >/dev/null 2>&1 || { 
   echo >&2 "flashpca is not available; set FPCA variable in 02pca.sh script or fix path"
   exit 1
}

mkdir PCA
pushd PCA

# Perform PCA within each population
for pop in CHB GIH JPT LWK MEX MKK YRI
do
   ROOT=../../data/hapmap3_r2_b36_fwd.qc.poly/hapmap3_r3_b36_fwd.${pop}.qc.poly.filtered.wexpr

   eval $PLINK \
      --bfile $ROOT \
      --indep-pairwise 1000 50 0.1 \
      --exclude range ~/Code/flashpca/exclusion_regions_hg19.txt \
      --out ${pop}

   eval $PLINK \
      --bfile $ROOT \
      --extract ${pop}.prune.in \
      --make-bed \
      --out $(basename $ROOT)_thinned

   eval $FPCA \
      --bfile $(basename $ROOT)_thinned \
      --mem low \
      --ndim 10 \
      --nextra 100 \
      --tol 1e-3 \
      --suffix _${pop}.txt \
      --v 2>&1 | tee $(basename $ROOT)_thinned.log
done

popd


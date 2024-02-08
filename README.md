# Merothon

Merothon is a collection of scripts designed for genetic data analysis, focusing on linkage disequilibrium calculations and genotype visualizations.

## Table of Contents

- [Installation](#installation)
- [Scripts](#scripts)
  - [Calculating R2 All SNPs, 2 VCFS](#calculating-r2-all-snps-2-vcfs)
  - [Plot Genotypes from VCF](#plot-genotypes-from-vcf)

## Installation

Most of the basic dependencies you likely already have, but it is easily installable with a fresh environment and setup.py:

```
git clone https://github.com/merondun/merothon.git
mamba create -n merothon python=3.8
mamba activate merothon
#cd merothon #or wherever the git repo was downloaded
pip install -e .
```

## Scripts

### Calculating R2 All SNPs, 2 VCFS

Calculates LD (R2) between the SNP genotypes in 2 VCF files. This is useful for estimating LD for e.g. a mtDNA variant and the rest of the autosomal SNPs. 

NOTE: This only works correctly for biallelic SNPs. It works for variable ploidy (same output as plink --ld-window 999999999 --ld-window-kb 100000000 --ld-window-r2 0), but it won't output meaingful results for 3N, 4N sites. 

Example from /examples/ directory: `calculate_r2 --vcf1 chr_MT_Biallelic_SNPs.vcf.gz --vcf2 chr_MT_Target_SNP.vcf.gz --out chr_MT_LD.txt`

If you have any invariant or constant sites in your VCF, you will get a warning "ConstantInputWarning: An input array is constant; the correlation coefficient is not defined.", but it does not affect calculations for other sites (output will be nan). 

| chrVCF1 | posVCF1 | chrVCF2 | posVCF2 | num_missing_genotypesVCF1 | num_missing_genotypesVCF2 | R2                         |
|---------|---------|---------|---------|---------------------------|---------------------------|----------------------------|
| chr_MT  | 44      | chr_MT  | 4270    | 0                         | 0                         | 0.3048275862068961         |
| chr_MT  | 192     | chr_MT  | 4270    | 0                         | 0                         | 0.20816783216783216        |
| chr_MT  | 196     | chr_MT  | 4270    | 0                         | 0                         | 0.19325217121588079        |

### Plot Genotypes from VCF

Plots color-coded genotypes for SNP positions. Provide a VCF file, a list of SNP positions (tab separated e.g. chr_MT	4270) a metadata file where $ID matches the sample names, and a phenotype of interest to distinguish individuals, as well as an output .png name.

Example from /examples/ directory: `plot_genotypes --vcf chr_MT_Biallelic_SNPs.vcf.gz --metadata Egg_Metadata.txt --pos Genotype_Inspect_Positions.txt --phenotype EggType --out Eggtype.png --size 50`

This will output a png showing genotypes via color, ordered by the EggType column in the metadata sheet. 

![Example Plot](examples/Eggtype.png)


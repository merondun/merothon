import argparse
import pysam
import numpy as np
from scipy.stats import pearsonr

def extract_genotypes(vcf_path, positions):
    vcf = pysam.VariantFile(vcf_path)
    records = []
    for chrom, pos in positions:
        for rec in vcf.fetch(str(chrom), pos-1, pos):
            # Check for biallelic sites
            if len(rec.alleles) != 2:
                # This is a non-biallelic site, mark as 'na'
                records.append((rec.chrom, rec.pos, 'na', 'na'))
                continue
            genotypes = [s['GT'] for s in rec.samples.values()]
            num_missing = sum(1 for g in genotypes if None in g)
            records.append((rec.chrom, rec.pos, genotypes, num_missing))
    return records

def calculate_r2(genotypes1, genotypes2):
    """
    Calculate R^2 (squared Pearson correlation coefficient) between two sets of genotypes.
    """

    # If either site is non-biallelic, return 'na'
    if genotypes1 == 'na' or genotypes2 == 'na':
        return 'na'
    
    # Convert genotypes to numeric arrays (0, 1, 2) and handle missing data
    g1 = np.array([sum(g) if None not in g else np.nan for g in genotypes1])
    g2 = np.array([sum(g) if None not in g else np.nan for g in genotypes2])
    
       # Check if either array is constant
    if len(set(g1[~np.isnan(g1)])) <= 1 or len(set(g2[~np.isnan(g2)])) <= 1:
        return np.nan  # Return NaN or another indicator value
    
    # Calculate Pearson correlation coefficient
    try:
        r, _ = pearsonr(g1[~np.isnan(g1) & ~np.isnan(g2)], g2[~np.isnan(g1) & ~np.isnan(g2)])
        return r**2
    except ConstantInputWarning:
        return np.nan  # Handle the warning gracefully

def read_genotypes(vcf_path):
    """
    Read genotypes from a VCF file.
    """
    vcf = pysam.VariantFile(vcf_path)
    snps = []
    for rec in vcf.fetch():
        genotypes = [s['GT'] for s in rec.samples.values()]
        num_missing = sum(1 for g in genotypes if None in g)
        snps.append((rec.chrom, rec.pos, genotypes, num_missing))
    return snps

def main(vcf1, vcf2, out):
    """
    Main function to calculate LD between SNPs in two VCF files.
    """
    snps_vcf1 = read_genotypes(vcf1)
    snps_vcf2 = read_genotypes(vcf2)
    
    with open(out, 'w') as out_f:
        out_f.write("chrVCF1\tposVCF1\tchrVCF2\tposVCF2\tnum_missing_genotypesVCF1\tnum_missing_genotypesVCF2\tR2\n")
        for snp1 in snps_vcf1:
            for snp2 in snps_vcf2:
                r2 = calculate_r2(snp1[2], snp2[2])
                out_f.write(f"{snp1[0]}\t{snp1[1]}\t{snp2[0]}\t{snp2[1]}\t{snp1[3]}\t{snp2[3]}\t{r2}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Calculate LD (R2) between all SNPs in two VCF files.')
    parser.add_argument('--vcf1', help='Path to the first VCF file',required=True)
    parser.add_argument('--vcf2', help='Path to the second VCF file',required=True)
    parser.add_argument('--out', help='Path to the output file',required=True)
    
    args = parser.parse_args()
    
    main(args.vcf1, args.vcf2, args.out)


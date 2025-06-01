import argparse
import pysam
import numpy as np
import sys

def extract_genotypes(vcf_path, positions):
    vcf = pysam.VariantFile(vcf_path)
    records = []
    for chrom, pos in positions:
        for rec in vcf.fetch(str(chrom), pos - 1, pos):
            if len(rec.alleles) != 2:
                records.append((rec.chrom, rec.pos, 'na', 'na'))
                continue
            genotypes = [s['GT'] for s in rec.samples.values()]
            num_missing = sum(1 for g in genotypes if None in g)
            records.append((rec.chrom, rec.pos, genotypes, num_missing))
    return records

def calculate_r2(genotypes1, genotypes2):
    if genotypes1 == 'na' or genotypes2 == 'na':
        return 'na'
    g1 = np.array([sum(g) if None not in g else np.nan for g in genotypes1])
    g2 = np.array([sum(g) if None not in g else np.nan for g in genotypes2])
    if len(set(g1[~np.isnan(g1)])) <= 1 or len(set(g2[~np.isnan(g2)])) <= 1:
        return np.nan
    mask = ~np.isnan(g1) & ~np.isnan(g2)
    g1, g2 = g1[mask], g2[mask]
    if len(g1) < 2 or np.std(g1) == 0 or np.std(g2) == 0:
        return np.nan
    r = np.corrcoef(g1, g2)[0, 1]
    return r**2

def read_genotypes(vcf_path):
    try:
        vcf = pysam.VariantFile(vcf_path)
    except ValueError as e:
        sys.stderr.write(f"Error reading VCF: {e}\n")
        sys.exit(1)
    snps = []
    for rec in vcf.fetch():
        genotypes = [s['GT'] for s in rec.samples.values()]
        num_missing = sum(1 for g in genotypes if None in g)
        snps.append((rec.chrom, rec.pos, genotypes, num_missing))
    return snps

def main():
    parser = argparse.ArgumentParser(description='Calculate LD (R²) between all SNPs in two VCF files.')
    parser.add_argument('--vcf1', required=True, help='Path to the first VCF file')
    parser.add_argument('--vcf2', required=True, help='Path to the second VCF file')
    parser.add_argument('--out', required=True, help='Path to the output file')
    args = parser.parse_args()

    snps_vcf1 = read_genotypes(args.vcf1)
    snps_vcf2 = read_genotypes(args.vcf2)

    with open(args.out, 'w') as out_f:
        out_f.write("chrVCF1\tposVCF1\tchrVCF2\tposVCF2\tnum_missing_genotypesVCF1\tnum_missing_genotypesVCF2\tR2\n")
        for snp1 in snps_vcf1:
            for snp2 in snps_vcf2:
                r2 = calculate_r2(snp1[2], snp2[2])
                out_f.write(f"{snp1[0]}\t{snp1[1]}\t{snp2[0]}\t{snp2[1]}\t{snp1[3]}\t{snp2[3]}\t{r2}\n")

    print(f"Successfully calculated LD between {args.vcf1} and {args.vcf2}")

if __name__ == "__main__":
    main()


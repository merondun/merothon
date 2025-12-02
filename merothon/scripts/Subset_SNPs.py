import argparse
import subprocess
import tempfile
import os
import sys

def get_total_snps(vcf_path, threads=1):
    try:
        result = subprocess.run(
            ["bcftools", "index", "-n", vcf_path],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=True
        )
        return int(result.stdout.strip())
    except subprocess.CalledProcessError:
        sys.stderr.write("Index not found, creating index...\n")
        subprocess.run(["bcftools", "index", "--threads", str(threads), vcf_path], check=True)
        result = subprocess.run(
            ["bcftools", "index", "-n", vcf_path],
            stdout=subprocess.PIPE,
            check=True
        )
        return int(result.stdout.strip())

def run_command(cmd):
    result = subprocess.run(cmd, shell=True, stderr=subprocess.PIPE)
    if result.returncode != 0:
        sys.stderr.write(result.stderr.decode())
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(description='Subset SNPs from a VCF file.')
    parser.add_argument('--vcf', required=True, help='Input VCF file (can be .vcf or .vcf.gz)')
    parser.add_argument('--out', required=True, help='Output VCF.GZ file path')
    parser.add_argument('--threads', type=int, default=1, help='Number of threads for bcftools (default: 1)')

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--prop', type=float, help='Proportion of SNPs to retain (0 < prop <= 1.0)')
    group.add_argument('--num', type=int, help='Number of SNPs to retain (integer > 0)')

    args = parser.parse_args()

    total_snps = get_total_snps(args.vcf, threads=args.threads)

    if args.prop is not None:
        if not (0 < args.prop <= 1.0):
            sys.stderr.write("--prop must be > 0 and <= 1.0\n")
            sys.exit(1)
        n_subset = int(total_snps * args.prop)
    else:
        if args.num <= 0:
            sys.stderr.write("--num must be > 0\n")
            sys.exit(1)
        if args.num > total_snps:
            sys.stderr.write(f"--num ({args.num}) exceeds total SNPs ({total_snps})\n")
            sys.exit(1)
        n_subset = args.num

    site_file = os.path.join(os.getcwd(), "subset_sites.txt")

    cmd = (
        f"bcftools query -f '%CHROM\\t%POS\\n' {args.vcf} | "
        f"shuf -n {n_subset} > {site_file}"
    )
    run_command(cmd)

    run_command(f"bcftools view --threads {args.threads} -T {site_file} -Oz -o {args.out} {args.vcf}")
    os.remove(site_file)

    print(f"Subsetted VCF from: {args.vcf}")
    print(f"Output VCF: {args.out}")
    if args.prop is not None:
        print(f"Proportion requested: {args.prop:.2f}")
    else:
        print(f"Number requested: {args.num}")
    print(f"Original SNP count: {total_snps}")
    print(f"Subset SNP count: {n_subset}")

if __name__ == "__main__":
    main()

import argparse
import subprocess
import tempfile
import os
import sys

def get_total_snps(vcf_path):
    cmd = ["bcftools", "view", "-v", "snps", vcf_path]
    p1 = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    p2 = subprocess.Popen(["bcftools", "view", "-H"], stdin=p1.stdout, stdout=subprocess.PIPE)
    p3 = subprocess.Popen(["wc", "-l"], stdin=p2.stdout, stdout=subprocess.PIPE)
    out, _ = p3.communicate()
    return int(out.strip())

def run_command(cmd):
    result = subprocess.run(cmd, shell=True, stderr=subprocess.PIPE)
    if result.returncode != 0:
        sys.stderr.write(result.stderr.decode())
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(description='Subset a proportion of SNPs from a VCF file.')
    parser.add_argument('--vcf', required=True, help='Input VCF file (can be .vcf or .vcf.gz)')
    parser.add_argument('--prop', type=float, required=True, help='Proportion of SNPs to retain (0 < prop <= 1.0)')
    parser.add_argument('--out', required=True, help='Output VCF.GZ file path')

    args = parser.parse_args()

    if not (0 < args.prop <= 1.0):
        sys.stderr.write("--prop must be > 0 and <= 1.0\n")
        sys.exit(1)

    total_snps = get_total_snps(args.vcf)
    n_subset = int(total_snps * args.prop)

    with tempfile.NamedTemporaryFile(mode='w+', delete=False) as temp_sites:
        site_file = temp_sites.name
        cmd = (
            f"bcftools view -v snps {args.vcf} | "
            f"bcftools view -H | shuf -n {n_subset} | "
            f"awk '{{print $1 \"\\t\" $2}}' > {site_file}"
        )
        run_command(cmd)

    run_command(f"bcftools view -T {site_file} -Oz -o {args.out} {args.vcf}")
    os.remove(site_file)

    print(f"Subsetted VCF from: {args.vcf}")
    print(f"Output VCF: {args.out}")
    print(f"Proportion requested: {args.prop:.2f}")
    print(f"Original SNP count: {total_snps}")
    print(f"Subset SNP count: {n_subset}")

if __name__ == "__main__":
    main()


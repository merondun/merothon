import argparse
import pysam

def parse_arguments():
    parser = argparse.ArgumentParser(description="Update VCF files with ancestral allele information.")
    parser.add_argument("--outgroups", required=True, help="File containing outgroup sample names.")
    parser.add_argument("--vcf", required=True, help="Input VCF file to be processed.")
    parser.add_argument("--out", required=True, help="Prefix for the output VCF file.")
    return parser.parse_args()

def load_outgroup_samples(outgroups):
    with open(outgroups, 'r') as file:
        outgroup_samples = [line.strip() for line in file]
    return outgroup_samples

def process_vcf(vcf, outgroup_samples, out):
    with pysam.VariantFile(vcf) as vcf_in:
        vcf_in.header.info.add('AA', '1', 'String', 'Ancestral Allele')
        with pysam.VariantFile(f"{out}.vcf.gz", 'w', header=vcf_in.header) as vcf_out:
            for record in vcf_in:
                alleles, ancestral_allele = analyze_record(record, outgroup_samples)
                record.info['AA'] = ancestral_allele
                vcf_out.write(record)

def analyze_record(record, outgroup_samples):
    alleles = []
    for sample in outgroup_samples:
        genotype = record.samples[sample]['GT']
        #check if genotype is not None and then proceed
        if genotype is None or len(genotype) == 0:
            # skip site genotype is None or empty
            continue
        # Adjust handling for haploid genotypes
        if len(genotype) == 1 or genotype[1] is None:  # Haploid or missing second allele
            allele_index = genotype[0]
        else:
            # if the genotype is missing or heterozygous in diploid, skip this sample
            if genotype.count(None) > 1 or set(genotype) == {0, 1}:
                alleles = []
                break
            allele_index = genotype[0]
        # Extract the allele using the determined index
        if allele_index is not None:  # Ensure allele_index is not None
            allele = record.alleles[allele_index]
            alleles.append(allele)
        else:
            # Handle cases where allele_index is None (shouldn't happen with the above checks, but just in case)
            continue

    if len(set(alleles)) > 1 or len(alleles) == 0:
        ancestral_allele = 'U'
    else:
        ancestral_allele = alleles[0]
    return alleles, ancestral_allele

def main():
    args = parse_arguments()
    outgroup_samples = load_outgroup_samples(args.outgroups)
    process_vcf(args.vcf, outgroup_samples, args.out)

if __name__ == "__main__":
    main()

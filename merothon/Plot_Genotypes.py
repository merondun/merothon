import argparse
import pysam
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Set up argument parsing
parser = argparse.ArgumentParser(description='Plot genotypes by position and phenotype.')
parser.add_argument('--vcf', help='Path to the VCF file', required=True)
parser.add_argument('--metadata', help='Path to the metadata file', required=True)
parser.add_argument('--pos', help='Path to the file containing positions of interest (tab-separated, chr and pos)', required=True)
parser.add_argument('--phenotype', help='Phenotype column name in the metadata file', required=True)
parser.add_argument('--out', help='Output PNG file name', required=True)
parser.add_argument('--size', help='Size of the markers in the plot [default 100]', type=int, default=100)

# Parse arguments
args = parser.parse_args()

# Load positions of interest from file
positions_of_interest = pd.read_csv(args.pos, sep='\t', header=None, names=['chrom', 'pos'])
positions_of_interest = [tuple(x) for x in positions_of_interest.to_numpy()]

# Read metadata
metadata = pd.read_csv(args.metadata, sep='\t')  # Assuming metadata is also tab-separated

# Function to extract genotypes for positions of interest
def extract_genotypes(vcf, positions):
    vcf = pysam.VariantFile(vcf)
    records = []
    for chrom, pos in positions:
        pos = int(pos)  # Ensure pos is an integer
        for rec in vcf.fetch(str(chrom), pos-1, pos):
            for sample in rec.samples:
                genotype = rec.samples[sample]['GT']
                genotype_str = '|'.join([str(allele) if allele is not None else '.' for allele in genotype])
                records.append({'chrom': chrom, 'pos': pos, 'ID': sample, 'Genotype': genotype_str})
    return pd.DataFrame(records)

# Extract genotypes
genotypes = extract_genotypes(args.vcf, positions_of_interest)

# Merge genotypes with metadata
merged_data = pd.merge(genotypes, metadata, on="ID")
merged_data['site'] = merged_data['chrom'] + ':' + merged_data['pos'].astype(str)
merged_data = merged_data.sort_values(by=args.phenotype)
# Plotting
plt.figure(figsize=(10, 6))
sns.set(style="whitegrid")
plot = sns.scatterplot(
    x='ID',
    y='site',
    hue='Genotype',
    style=args.phenotype,
    palette='Set2',
    data=merged_data,
    s=args.size
)

plt.xticks(rotation=90)
plt.xlabel("Individuals")
plt.ylabel("Position")
plt.title(f"Genotypes by Position and {args.phenotype}")
plt.legend(title='Genotype', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()

# Save the plot to a file
plt.savefig(args.out)
plt.close()

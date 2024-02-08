import argparse
import pysam
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Define file paths and positions of interest
vcf_path = "/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/immaculate_hunting/vcfs_raw/chr_MT.vcf.gz"
metadata_path = "/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/Cuckoo_Full_Metadata_2023OCT3.txt"
positions_of_interest = [('chr_MT', 4646), ('chr_MT', 4270)] # List of tuples (chromosome, position)
phenotype = 'Immaculate'
output_figure = 'Output.png'

# Read metadata
metadata = pd.read_csv(metadata_path,delimiter='\t')

# Function to extract genotypes for positions of interest
def extract_genotypes(vcf_path, positions):
    vcf = pysam.VariantFile(vcf_path)
    records = []
    for chrom, pos in positions:
        for rec in vcf.fetch(str(chrom), pos-1, pos):
            for sample in rec.samples:
                genotype = rec.samples[sample]['GT']
                genotype_str = '|'.join([str(allele) if allele is not None else '.' for allele in genotype])
                records.append({'chrom': chrom, 'pos': pos, 'ID': sample, 'Genotype': genotype_str})
    return pd.DataFrame(records)

# Extract genotypes
genotypes = extract_genotypes(vcf_path, positions_of_interest)
print(genotypes)

# Merge genotypes with metadata
merged_data = pd.merge(genotypes, metadata, on="ID")

# Prepare the data for plotting
merged_data['site'] = merged_data['chrom'].astype(str) + ':' + merged_data['pos'].astype(str)
merged_data['Genotype'] = merged_data['Genotype'].astype('category')
merged_data = merged_data.sort_values(by=phenotype)
print(merged_data)

# Plotting
plt.figure(figsize=(10, 6))
sns.set(style="whitegrid")

plot = sns.scatterplot(
    x='ID',
    y='site',
    hue='Genotype',
    style=phenotype, # Assumes 'phenotype' column in metadata
    palette='Set2',
    data=merged_data,
    s=100
)

plt.xticks(rotation=90)
plt.xlabel("Individuals")
plt.ylabel("Position")
plt.title("Genotypes by Position and Phenotype")
plt.legend(title='Genotype', bbox_to_anchor=(1.05, 1), loc='upper left')

plt.tight_layout()

# Save the plot to a PNG file
plt.savefig(output_figure)

# Optionally, after saving you can also close the plot if you don't need it displayed in your environment
plt.close()


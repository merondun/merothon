import argparse
import pysam
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch
from matplotlib.colors import ListedColormap

def main():
    parser = argparse.ArgumentParser(description='Tile plot of genotypes by position and phenotype.')
    parser.add_argument('--vcf', help='Path to the VCF file', required=True)
    parser.add_argument('--metadata', help='Path to the metadata file', required=True)
    parser.add_argument('--pos', help="Position of interest, formatted as 'chr:pos'", required=True)
    parser.add_argument('--phenotype', help='Phenotype column name in the metadata file', required=True)
    parser.add_argument('--out', help='Output PNG file name', required=True)

    args = parser.parse_args()

    positions_of_interest = [(item.split(':')[0], int(item.split(':')[1])) for item in args.pos.split(',')]
    metadata = pd.read_csv(args.metadata, sep='\t')
    genotypes = extract_genotypes(args.vcf, positions_of_interest)
    merged_data = pd.merge(genotypes, metadata, on="ID")
    merged_data['Phenotype.ID'] = merged_data[args.phenotype] + '.' + merged_data['ID']

    plot_data = prepare_plot_data(merged_data)

    plot_genotype_counts(plot_data, merged_data, args.out)

def extract_genotypes(vcf, positions):
    vcf = pysam.VariantFile(vcf)
    records = []
    for chrom, pos in positions:
        for rec in vcf.fetch(str(chrom), pos-1, pos):
            for sample in rec.samples:
                genotype = rec.samples[sample]['GT']
                genotype_str = '|'.join([str(allele) for allele in genotype]) if genotype else 'Missing'
                records.append({'chrom': chrom, 'pos': str(pos), 'ID': sample, 'Genotype': genotype_str})
    return pd.DataFrame(records)

def prepare_plot_data(merged_data):
    # Convert 'Genotype' to categorical and map to codes, including 'Missing' handling
    merged_data['GenotypeValue'] = pd.Categorical(merged_data['Genotype']).codes
    pivot_table = merged_data.pivot_table(index='Phenotype.ID', columns='pos', values='GenotypeValue', aggfunc='first', fill_value=-1)
    return pivot_table

def plot_genotype_counts(plot_data, merged_data, output_file):
    # Identify unique genotypes and assign each a discrete color from 'viridis'
    unique_genotypes = np.unique(merged_data['Genotype'])
    n_genotypes = len(unique_genotypes)
    
    # Generate a colormap with enough colors for each genotype
    viridis_colors = plt.cm.viridis(np.linspace(0, 1, n_genotypes))
    genotype_to_color = {genotype: viridis_colors[i] for i, genotype in enumerate(unique_genotypes)}

    # Apply colormap to genotypes in plot_data
    plot_data_colored = plot_data.replace({code: genotype_to_color[genotype] for genotype, code in zip(unique_genotypes, range(n_genotypes))})

    fig, ax = plt.subplots(figsize=(12, 8))
    # Plot each cell with the color corresponding to its genotype
    for (i, j), genotype_code in np.ndenumerate(plot_data):
        ax.fill_betweenx([i-0.5, i+0.5], j-0.5, j+0.5, color=genotype_to_color.get(unique_genotypes[genotype_code], 'white'))

    # Set axis ticks
    ax.set_xticks(range(len(plot_data.columns)))
    ax.set_xticklabels(plot_data.columns, rotation=90)
    ax.set_yticks(range(len(plot_data.index)))
    ax.set_yticklabels(plot_data.index)

    # Create a legend
    legend_elements = [Patch(facecolor=genotype_to_color[geno], label=geno) for geno in unique_genotypes]
    ax.legend(handles=legend_elements, bbox_to_anchor=(1.05, 1), loc='upper left', title="Genotype")

    plt.xlabel('Position')
    plt.ylabel('Individual.Phenotype')
    plt.title('Genotypes by Position and Phenotype')
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()

if __name__ == "__main__":
    main()

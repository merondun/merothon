import argparse
import allel
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def parse_arguments():
    parser = argparse.ArgumentParser(description='Perform PCA on VCF data, plot colored by phenotype, and optionally label points and write PCA results.')
    parser.add_argument('--vcf', required=True, help='Path to the gzipped VCF file.')
    parser.add_argument('--metadata', required=True, help='Path to the metadata file (tab-separated, with $ID and additional columns).')
    parser.add_argument('--phenotype', required=True, help='The column in the metadata file to use for coloring PCA points.')
    parser.add_argument('--out', required=True, help='Output file path (supports .png or .pdf formats).')
    parser.add_argument('--label', action='store_true', help='Optionally label PCA points with individual IDs.')
    parser.add_argument('--write', action='store_true', help='Optionally write PCA eigenvectors and eigenvalues to files.')
    return parser.parse_args()

def read_population_map(metadata_file, phenotype_col):
    print(f"Reading population map from {metadata_file}...")
    df = pd.read_csv(metadata_file, sep='\t')
    df[phenotype_col] = df[phenotype_col].astype('category')
    return df

def perform_pca(vcf_path):
    print(f"Reading VCF data from {vcf_path}...")
    callset = allel.read_vcf(vcf_path)
    samples = callset['samples']
    print("Performing PCA analysis...")
    gt = allel.GenotypeArray(callset['calldata/GT'])
    ac = gt.count_alleles()
    flt = ac.is_segregating() & (ac.max_allele() == 1)
    gf = gt.compress(flt, axis=0).to_n_alt()
    coords, model = allel.pca(gf, n_components=4, scaler='patterson')  #compute 4 components, modify if you want 
    return coords, samples, model

def plot_pca(coords, metadata, phenotype_col, samples, output_file, label_ids, model):
    print("Plotting PCA results...")
    fig, axs = plt.subplots(1, 2, figsize=(8, 5))  # Adjusted figure size

    #Ensure metadata is filtered to match samples and phenotype categories
    metadata_filtered = metadata[metadata['ID'].isin(samples)]
    metadata_filtered = metadata_filtered.set_index('ID').loc[samples].reset_index()
    metadata_filtered[phenotype_col] = pd.Categorical(metadata_filtered[phenotype_col])

    # Remove unused categories after filtering
    metadata_filtered[phenotype_col] = metadata_filtered[phenotype_col].cat.remove_unused_categories()

    categories = metadata_filtered[phenotype_col].cat.categories
    cmap = plt.get_cmap('viridis', len(categories))

    # Plotting PCA results for the first four components in two subplots
    for ax, pcs in zip(axs, [(0, 1), (2, 3)]):
        for i, category in enumerate(categories):
            indices = metadata_filtered[phenotype_col] == category
            ax.scatter(coords[indices, pcs[0]], coords[indices, pcs[1]], label=category, color=cmap(i))
            if label_ids:
                for txt, x, y in zip(metadata_filtered.loc[indices, 'ID'], coords[indices, pcs[0]], coords[indices, pcs[1]]):
                    ax.annotate(txt, (x, y))
        ax.set_xlabel(f'PC{pcs[0]+1}: {model.explained_variance_ratio_[pcs[0]]*100:.2f}%')
        ax.set_ylabel(f'PC{pcs[1]+1}: {model.explained_variance_ratio_[pcs[1]]*100:.2f}%')

    #Adjust legend placement and figure layout
    handles, labels = axs[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 0.1), fancybox=True, shadow=True, ncol=3, fontsize='small')
    plt.suptitle('PCA by Phenotype')
    
    #adjust layout to make space for the legend outside the subplots
    plt.tight_layout(rect=[0, 0.1, 1, 0.95])  # Adjust the bottom parameter to accommodate the legend
    plt.savefig(output_file, format=output_file.split('.')[-1])

def write_pca_results(coords, metadata, phenotype_col, samples, prefix):
    # Filter and reorder metadata to match the PCA results
    metadata_filtered = metadata[metadata['ID'].isin(samples)]
    metadata_filtered = metadata_filtered.set_index('ID').loc[samples].reset_index()

    # Create a DataFrame with PCA results
    pca_results = pd.DataFrame({
        'ID': metadata_filtered['ID'],
        phenotype_col: metadata_filtered[phenotype_col],
        'PC1': coords[:, 0],
        'PC2': coords[:, 1],
        'PC3': coords[:, 2],
        'PC4': coords[:, 3]
    })

    #Write the DataFrame to a txt file 
    pca_results.to_csv(f"{prefix}_PCA_results.txt",sep='\t', index=False)
    print(f"PCA results written to {prefix}_PCA_results.txt.")

def main():
    args = parse_arguments()
    metadata = read_population_map(args.metadata, args.phenotype)
    coords, samples, model = perform_pca(args.vcf)
    plot_pca(coords, metadata, args.phenotype, samples, args.out, args.label, model)
    if args.write:
        write_pca_results(coords, metadata, args.phenotype, samples, args.out.rsplit('.', 1)[0])

if __name__ == '__main__':
    main()

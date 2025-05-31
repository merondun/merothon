import argparse
import allel
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def parse_arguments():
    parser = argparse.ArgumentParser(description='Perform PCA on VCF data, plot colored by phenotype if provided, and always write PCA results.')
    parser.add_argument('--vcf', required=True, help='Path to the gzipped VCF file.')
    parser.add_argument('--metadata', help='Path to the metadata file (tab-separated, with $ID and additional columns). Optional.')
    parser.add_argument('--phenotype', help='The column in the metadata file to use for coloring PCA points. Required if --metadata is provided.')
    parser.add_argument('--out', required=True, help='Output file path (supports .png or .pdf formats). Prefix for output files.')
    parser.add_argument('--label', action='store_true', help='Optionally label PCA points with individual IDs.')
    return parser.parse_args()

def read_population_map(metadata_file, phenotype_col):
    if metadata_file and phenotype_col:
        print(f"Reading population map from {metadata_file}...")
        df = pd.read_csv(metadata_file, sep='\t')
        df[phenotype_col] = df[phenotype_col].astype('category')
        return df
    return None

def perform_pca(vcf_path):
    print(f"Reading VCF data from {vcf_path}...")
    callset = allel.read_vcf(vcf_path)
    samples = callset['samples']
    print("Performing PCA analysis...")
    gt = allel.GenotypeArray(callset['calldata/GT'])
    ac = gt.count_alleles()
    flt = ac.is_segregating() & (ac.max_allele() == 1)
    gf = gt.compress(flt, axis=0).to_n_alt()
    coords, model = allel.pca(gf, n_components=4, scaler='patterson')
    return coords, samples, model

def plot_pca(coords, metadata, phenotype_col, samples, output_file, label_ids, model):
    print("Plotting PCA results...")
    fig, axs = plt.subplots(1, 2, figsize=(8, 5))

    if metadata is not None and phenotype_col:
        # Ensure metadata is filtered to match samples and phenotype categories
        metadata_filtered = metadata[metadata['ID'].isin(samples)]
        metadata_filtered = metadata_filtered.set_index('ID').loc[samples].reset_index()
        metadata_filtered[phenotype_col] = pd.Categorical(metadata_filtered[phenotype_col])
        metadata_filtered[phenotype_col] = metadata_filtered[phenotype_col].cat.remove_unused_categories()

        categories = metadata_filtered[phenotype_col].cat.categories
        cmap = plt.get_cmap('viridis', len(categories))

        for ax, pcs in zip(axs, [(0, 1), (2, 3)]):
            for i, category in enumerate(categories):
                indices = metadata_filtered[phenotype_col] == category
                ax.scatter(coords[indices, pcs[0]], coords[indices, pcs[1]], label=category, color=cmap(i))
                if label_ids:
                    for txt, x, y in zip(metadata_filtered.loc[indices, 'ID'], coords[indices, pcs[0]], coords[indices, pcs[1]]):
                        ax.annotate(txt, (x, y))
            ax.set_xlabel(f'PC{pcs[0]+1}: {model.explained_variance_ratio_[pcs[0]]*100:.2f}%')
            ax.set_ylabel(f'PC{pcs[1]+1}: {model.explained_variance_ratio_[pcs[1]]*100:.2f}%')

        handles, labels = axs[0].get_legend_handles_labels()
        fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 0.1), fancybox=True, shadow=True, ncol=3, fontsize='small')
    else:
        for ax, pcs in zip(axs, [(0, 1), (2, 3)]):
            ax.scatter(coords[:, pcs[0]], coords[:, pcs[1]])
            ax.set_xlabel(f'PC{pcs[0]+1}: {model.explained_variance_ratio_[pcs[0]]*100:.2f}%')
            ax.set_ylabel(f'PC{pcs[1]+1}: {model.explained_variance_ratio_[pcs[1]]*100:.2f}%')

    plt.suptitle('PCA by Phenotype' if phenotype_col else 'PCA Results')
    plt.tight_layout(rect=[0, 0.1, 1, 0.95])
    plt.savefig(output_file, format=output_file.split('.')[-1])

def write_pca_results(coords, samples, model, prefix):
    pca_scores = pd.DataFrame(coords, columns=['PC1', 'PC2', 'PC3', 'PC4'], index=samples)
    pca_scores.to_csv(f"{prefix}_scores.txt", sep='\t', index_label='SampleID')

    explained_variance = pd.Series(model.explained_variance_ratio_, index=['PC1', 'PC2', 'PC3', 'PC4'])
    explained_variance.to_csv(f"{prefix}_values.txt", sep='\t', index_label='PC',header=['Explained_Variance'])

def main():
    args = parse_arguments()
    if args.metadata and not args.phenotype:
        raise ValueError("--phenotype must be provided if --metadata is specified.")
    metadata = read_population_map(args.metadata, args.phenotype) if args.metadata else None
    coords, samples, model = perform_pca(args.vcf)
    plot_pca(coords, metadata, args.phenotype, samples, args.out, args.label, model)
    output_prefix = args.out.rsplit('.', 1)[0]
    write_pca_results(coords, samples, model, output_prefix)

if __name__ == '__main__':
    main()

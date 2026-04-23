import argparse
from types import SimpleNamespace

from cyvcf2 import VCF
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


def parse_arguments():
    parser = argparse.ArgumentParser(
        description='Perform PCA on VCF data, plot colored by phenotype if provided, and always write PCA results.'
    )
    parser.add_argument('--vcf', required=True, help='Path to the gzipped VCF file.')
    parser.add_argument('--metadata', help='Path to the metadata file (tab-separated, with ID and additional columns). Optional.')
    parser.add_argument('--phenotype', help='The column in the metadata file to use for coloring PCA points. Required if --metadata is provided.')
    parser.add_argument('--out', required=True, help='Output file path (supports .png or .pdf formats). Prefix for output files.')
    parser.add_argument('--label', action='store_true', help='Optionally label PCA points with individual IDs.')
    return parser.parse_args()


def read_population_map(metadata_file, phenotype_col):
    if metadata_file and phenotype_col:
        print(f"Reading population map from {metadata_file}...")
        df = pd.read_csv(metadata_file, sep='\t')
        if 'ID' not in df.columns:
            raise ValueError("Metadata file must contain an 'ID' column.")
        if phenotype_col not in df.columns:
            raise ValueError(f"Phenotype column '{phenotype_col}' not found in metadata.")
        df[phenotype_col] = df[phenotype_col].astype('category')
        return df
    return None


def _patterson_scale(g):
    # g is variants x samples, values in {0,1,2,np.nan}
    p = np.nanmean(g, axis=1) / 2.0
    keep = np.isfinite(p) & (p > 0.0) & (p < 1.0)
    g = g[keep]
    p = p[keep]

    if g.shape[0] == 0:
        raise ValueError("No segregating biallelic SNPs remain after filtering.")

    # mean-impute missing genotypes by expected genotype count 2p
    means = 2.0 * p
    nan_mask = np.isnan(g)
    if nan_mask.any():
        g = g.copy()
        g[nan_mask] = np.take(means, np.where(nan_mask)[0])

    denom = np.sqrt(2.0 * p * (1.0 - p))
    keep2 = denom > 0
    g = g[keep2]
    means = means[keep2]
    denom = denom[keep2]

    x = (g - means[:, None]) / denom[:, None]
    return x


def _svd_pca(gf, n_components=4):
    # gf is variants x samples
    x = _patterson_scale(gf)

    # sample x variants for PCA coordinates by sample
    x = x.T

    n_samples = x.shape[0]
    if n_samples < 2:
        raise ValueError("Need at least two samples for PCA.")

    max_components = min(n_components, x.shape[0], x.shape[1])
    if max_components < n_components:
        print(f"Warning: reducing n_components from {n_components} to {max_components}.")
        n_components = max_components

    u, s, vh = np.linalg.svd(x, full_matrices=False)
    coords = u[:, :n_components] * s[:n_components]

    # explained variance from singular values
    explained_variance = (s ** 2) / (n_samples - 1)
    total_variance = explained_variance.sum()
    explained_variance_ratio = explained_variance[:n_components] / total_variance

    model = SimpleNamespace(
        explained_variance_ratio_=explained_variance_ratio
    )
    return coords, model


def perform_pca(vcf_path):
    print(f"Reading VCF data from {vcf_path}...")

    vcf = VCF(vcf_path, gts012=True, strict_gt=True)
    samples = np.array(vcf.samples)

    gt_rows = []
    n_total = 0
    n_kept = 0

    for variant in vcf:
        n_total += 1

        if not variant.is_snp:
            continue
        if len(variant.REF) != 1:
            continue
        if len(variant.ALT) != 1:
            continue
        if variant.ALT[0] is None or len(variant.ALT[0]) != 1:
            continue

        gt = variant.gt_types.astype(float)  # 0,1,2,3 where 3 is UNKNOWN
        gt[gt == 3] = np.nan

        # keep segregating biallelic SNPs
        if np.all(np.isnan(gt)):
            continue

        finite = gt[np.isfinite(gt)]
        if finite.size == 0:
            continue
        if np.nanmin(gt) == np.nanmax(gt):
            continue

        gt_rows.append(gt)
        n_kept += 1

    if n_kept == 0:
        raise ValueError("No usable biallelic segregating SNPs found in VCF.")

    gf = np.vstack(gt_rows)  # variants x samples

    print(f"Variants scanned: {n_total}")
    print(f"Biallelic segregating SNPs kept: {n_kept}")
    print("Performing PCA analysis...")

    coords, model = _svd_pca(gf, n_components=4)
    return coords, samples, model


def plot_pca(coords, metadata, phenotype_col, samples, output_file, label_ids, model):
    print("Plotting PCA results...")
    n_axes = 2 if coords.shape[1] >= 4 else 1
    fig, axs = plt.subplots(1, n_axes, figsize=(8 if n_axes == 2 else 4.5, 5))
    if n_axes == 1:
        axs = [axs]

    pc_pairs = [(0, 1)] if coords.shape[1] < 4 else [(0, 1), (2, 3)]

    if metadata is not None and phenotype_col:
        metadata_filtered = metadata[metadata['ID'].isin(samples)].copy()
        metadata_filtered = metadata_filtered.set_index('ID').reindex(samples).reset_index()

        if metadata_filtered[phenotype_col].isna().any():
            missing = metadata_filtered.loc[metadata_filtered[phenotype_col].isna(), 'ID'].tolist()
            raise ValueError(
                f"Some VCF samples are missing phenotype metadata in column '{phenotype_col}': "
                + ", ".join(missing[:10])
                + ("..." if len(missing) > 10 else "")
            )

        metadata_filtered[phenotype_col] = pd.Categorical(metadata_filtered[phenotype_col])
        metadata_filtered[phenotype_col] = metadata_filtered[phenotype_col].cat.remove_unused_categories()

        categories = metadata_filtered[phenotype_col].cat.categories
        cmap = plt.get_cmap('viridis', len(categories))

        for ax, pcs in zip(axs, pc_pairs):
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
        for ax, pcs in zip(axs, pc_pairs):
            ax.scatter(coords[:, pcs[0]], coords[:, pcs[1]])
            ax.set_xlabel(f'PC{pcs[0]+1}: {model.explained_variance_ratio_[pcs[0]]*100:.2f}%')
            ax.set_ylabel(f'PC{pcs[1]+1}: {model.explained_variance_ratio_[pcs[1]]*100:.2f}%')

    plt.suptitle('PCA by Phenotype' if phenotype_col else 'PCA Results')
    plt.tight_layout(rect=[0, 0.1, 1, 0.95] if metadata is not None and phenotype_col else [0, 0, 1, 0.95])
    plt.savefig(output_file, format=output_file.split('.')[-1])


def write_pca_results(coords, samples, model, prefix):
    n_pcs = coords.shape[1]
    pc_names = [f'PC{i+1}' for i in range(n_pcs)]

    pca_scores = pd.DataFrame(coords, columns=pc_names, index=samples)
    pca_scores.to_csv(f"{prefix}_scores.txt", sep='\t', index_label='SampleID')

    explained_variance = pd.Series(model.explained_variance_ratio_, index=pc_names)
    explained_variance.to_csv(
        f"{prefix}_values.txt",
        sep='\t',
        index_label='PC',
        header=['Explained_Variance']
    )


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
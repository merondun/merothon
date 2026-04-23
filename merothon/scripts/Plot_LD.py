import argparse
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np

def bin_ld_data(df, window_size):
    df['BP_A_bin'] = (df['BP_A'] // (window_size * 1000)) * (window_size * 1000)
    df['BP_B_bin'] = (df['BP_B'] // (window_size * 1000)) * (window_size * 1000)
    df = df.groupby(['BP_A_bin', 'BP_B_bin'], as_index=False)['R2'].mean()
    return df

def plot_ld_heatmap(input_file, out, win_size, highlight_region=None):
    df = pd.read_csv(input_file, sep='\s+')
    df['BP_A'] = pd.to_numeric(df['BP_A'])
    df['BP_B'] = pd.to_numeric(df['BP_B'])
    df['R2'] = pd.to_numeric(df['R2'])
    df_binned = bin_ld_data(df, win_size)
    binned_output_file = f"{out}_windows.txt"
    df_binned.to_csv(binned_output_file, sep='\t', index=False)

    ld_matrix = df_binned.pivot(index='BP_B_bin', columns='BP_A_bin', values='R2')
    plt.figure(figsize=(10, 7))
    im = plt.imshow(ld_matrix.values, cmap='YlOrRd', aspect='equal', origin='lower',
                    extent=[ld_matrix.columns.min(), ld_matrix.columns.max(), 
                            ld_matrix.index.min(), ld_matrix.index.max()])
    cbar = plt.colorbar(im)
    cbar.set_label('R²')
    plt.xlabel(f"BP_A ({win_size} KB bins)")
    plt.ylabel(f"BP_B ({win_size} KB bins)")
    plt.title("Linkage Disequilibrium (LD) Heatmap")

    if highlight_region:
        start, end = highlight_region
        rect = patches.Rectangle((start, start), end - start, end - start,
                                 linewidth=2, edgecolor='blue', facecolor='none')
        plt.gca().add_patch(rect)

    if not out.lower().endswith((".png", ".pdf")):
        out += ".png"

    plt.savefig(out, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"LD plotted successfully, output within {out}")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Input LD file")
    parser.add_argument("--out", required=True, help="Output file name (with .png or .pdf extension)")
    parser.add_argument("--win_size", type=int, required=True, help="Window size in KB to average R2 within (recommended 5)")
    parser.add_argument("--highlight", help="Highlight region as start-end (e.g., 700-900)")
    args = parser.parse_args()

    highlight_region = None
    if args.highlight:
        try:
            start, end = map(int, args.highlight.split('-'))
            highlight_region = (start, end)
        except ValueError:
            raise ValueError("Invalid format for --highlight. Use start-end (e.g., 700-900)")

    plot_ld_heatmap(args.input, args.out, args.win_size, highlight_region)

if __name__ == "__main__":
    main()


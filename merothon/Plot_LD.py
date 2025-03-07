import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.patches as patches

def bin_ld_data(df, window_size):
    # Bin BP_A and BP_B into N KB windows
    df['BP_A_bin'] = (df['BP_A'] // (window_size * 1000)) * (window_size * 1000)
    df['BP_B_bin'] = (df['BP_B'] // (window_size * 1000)) * (window_size * 1000)

    # Average R2 values within each window
    df = df.groupby(['BP_A_bin', 'BP_B_bin'], as_index=False)['R2'].mean()
    return df

def plot_ld_heatmap(input_file, out, win_size, highlight_region=None):
    # Read the LD file, skipping the first row (header)
    df = pd.read_csv(input_file, sep='\s+')

    # Convert columns to numeric
    df['BP_A'] = pd.to_numeric(df['BP_A'])
    df['BP_B'] = pd.to_numeric(df['BP_B'])
    df['R2'] = pd.to_numeric(df['R2'])

    # Bin data into N KB windows
    df_binned = bin_ld_data(df, win_size)

    # Save the binned data
    binned_output_file = f"{out}_windows.txt"
    df_binned.to_csv(binned_output_file, sep='\t', index=False)

    # Pivot table for heatmap
    ld_matrix = df_binned.pivot(index='BP_B_bin', columns='BP_A_bin', values='R2')

    # Plot
    plt.figure(figsize=(10, 7))
    ax = sns.heatmap(ld_matrix, cmap='YlOrRd', cbar_kws={'label': 'R²'}, square=True)
    plt.xlabel(f"BP_A ({win_size} KB bins)")
    plt.ylabel(f"BP_B ({win_size} KB bins)")
    plt.title("Linkage Disequilibrium (LD) Heatmap")

    # Highlight a specific region (chr:start-end)
    if highlight_region:
        start, end = highlight_region  # Expecting a tuple (start, end)
        ax = plt.gca()

        # Find closest indices for given positions
        bp_a_values = ld_matrix.columns.to_numpy()
        bp_b_values = ld_matrix.index.to_numpy()

        start_idx = (abs(bp_a_values - start)).argmin()
        end_idx = (abs(bp_a_values - end)).argmin()

        # Draw rectangle outline
        rect = patches.Rectangle((start_idx, start_idx), end_idx - start_idx, end_idx - start_idx,
                                 linewidth=2, edgecolor='blue', facecolor='none')
        ax.add_patch(rect)

    # Determine file extension and save plot accordingly
    if not out.lower().endswith((".png", ".pdf")):
        out += ".png"  # Default to PNG if no valid extension is provided

    plt.savefig(out, dpi=300, bbox_inches='tight')
    plt.close()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Input LD file")
    parser.add_argument("--out", required=True, help="Output file name (with .png or .pdf extension)")
    parser.add_argument("--win_size", type=int, required=True, help="Window size in KB to average R2 within (recommended 5)")
    parser.add_argument("--highlight", help="Highlight region as start-end (e.g., 700-900)")
    args = parser.parse_args()

    # Parse highlight region if provided
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

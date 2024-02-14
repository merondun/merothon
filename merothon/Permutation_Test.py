import argparse
import pandas as pd
import numpy as np

def load_data(all_data_path, regions_path):
    all_data = pd.read_csv(all_data_path, sep='\t', header=None, names=['chr', 'start', 'end', 'value'], dtype={"chr": str, "start": int, "end": int, "value": float})
    regions = pd.read_csv(regions_path, sep='\t', header=None, names=['chr', 'start', 'end', 'name'], dtype={"chr": str, "start": int, "end": int, "name": str})
    return all_data, regions

def check_overlap(row, region):
    return row['chr'] == region['chr'] and not (row['start'] > region['end'] or row['end'] < region['start'])

def perform_permutation_test(all_data, region, permutations, seed=None):
    if seed is not None:
        np.random.seed(seed)

    print(f"Permuting region: {region['name']}")

    target = all_data[all_data.apply(check_overlap, axis=1, args=(region,))]
    if target.empty:
        return None  #return NA if no overlaps 
        print(f"No overlaps for region: {region['name']}")

    background = all_data[~all_data.apply(check_overlap, axis=1, args=(region,))]
    
    observed_mean_diff = target['value'].mean() - background['value'].mean()
    num_target_windows = len(target)

    permutation_diffs = []
    combined_values = pd.concat([target['value'], background['value']]).reset_index(drop=True)
    for _ in range(permutations):
        shuffled = combined_values.sample(frac=1).reset_index(drop=True)
        new_target_mean = shuffled.iloc[:num_target_windows].mean()
        new_background_mean = shuffled.iloc[num_target_windows:].mean()
        mean_diff = new_target_mean - new_background_mean
        permutation_diffs.append(mean_diff)

    return observed_mean_diff, permutation_diffs, num_target_windows

def main():
    parser = argparse.ArgumentParser(description="Perform a permutation test on genomic data.")
    parser.add_argument("--all_data", required=True, help="File path to all genomic data.")
    parser.add_argument("--regions", required=True, help="File path to genomic regions of interest.")
    parser.add_argument("--permutations", type=int, default=1000, help="Number of permutations to perform.")
    parser.add_argument("--out", required=True, help="Output file path for permutation test results.")
    parser.add_argument("--seed", type=int, help="Seed for random number generation to ensure reproducible results.", default=None)

    args = parser.parse_args()

    all_data, regions = load_data(args.all_data, args.regions)
    results = []

    for _, region in regions.iterrows():
        result = perform_permutation_test(all_data, region, args.permutations, args.seed)
        if result is None:
            # Handle the case where no overlaps are found
            print(f"No overlaps for region: {region['name']}")
            results.append({'name': region['name'], 'observed_difference': 'NA', 'permuted_difference': 'NA', 'num_target_windows': 0})
        else:
            observed_mean_diff, permutation_diffs, num_target_windows = result
            for perm_diff in permutation_diffs:
                results.append({'name': region['name'], 'observed_difference': observed_mean_diff, 'permuted_difference': perm_diff, 'num_target_windows': num_target_windows})

    results_df = pd.DataFrame(results)
    results_df.to_csv(args.out, sep='\t', index=False)

if __name__ == "__main__":
    main()
import argparse
import pandas as pd
import numpy as np

def load_data(all_data_path, regions_path):
    # Load all_data and regions without headers, specifying column names
    all_data = pd.read_csv(all_data_path, sep='\t', header=None, names=['chr', 'start', 'end', 'value'], dtype={"chr": str, "start": int, "end": int, "value": float})
    regions = pd.read_csv(regions_path, sep='\t', header=None, names=['chr', 'start', 'end', 'name'], dtype={"chr": str, "start": int, "end": int, "name": str})
    return all_data, regions

def check_overlap(row, region):
    # Check if row overlaps with region
    return row['chr'] == region['chr'] and not (row['start'] > region['end'] or row['end'] < region['start'])

def sample_difference(all_data, region, events):
    valid_differences = 0
    differences = []
    
    while valid_differences < events:
        target = all_data[all_data.apply(check_overlap, axis=1, args=(region,))]
        background = all_data[~all_data.apply(check_overlap, axis=1, args=(region,))]

        # Sample one value from target and one from background
        target_sample = target['value'].sample(1) if not target.empty else pd.Series([np.nan])
        background_sample = background['value'].sample(1) if not background.empty else pd.Series([np.nan])

        # Calculate the difference
        diff = target_sample.values[0] - background_sample.values[0]

        # Check if diff is NaN
        if not np.isnan(diff):
            differences.append(diff)
            valid_differences += 1

    return differences

def bootstraps(all_data, regions, events):
    results = []
    for _, region in regions.iterrows():
        diffs = sample_difference(all_data, region, events)
        for diff in diffs:
            results.append({'name': region['name'], 'difference': diff})

    return pd.DataFrame(results)

def main():
    parser = argparse.ArgumentParser(description="Perform bootstramp sampling on genomic data. Calculates difference of values (within region) - values (outside region) on the same chromosome. ")
    parser.add_argument("--all_data", required=True, help="File path to all genomic data. Must be chr start end value, tab separated ")
    parser.add_argument("--regions", required=True, help="File path to genomic regions of interest. Msut be chr start end region_name, tab separated ")
    parser.add_argument("--events", type=int, default=10000, help="Number of bootstrap sampling events to perform. [default 10000]")
    parser.add_argument("--out", required=True, help="Output file path for bootstrap sampling test results.")
    parser.add_argument("--seed", type=int, help="Seed for random number generation to ensure reproducible results.", default=None)

    args = parser.parse_args()

    if args.seed is not None:
        np.random.seed(args.seed)
        
    all_data, regions = load_data(args.all_data, args.regions)
    results = bootstraps(all_data, regions, args.events)
    results.to_csv(args.out, sep='\t', index=False)

if __name__ == "__main__":
    main()

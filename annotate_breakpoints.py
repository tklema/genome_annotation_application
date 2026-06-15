import argparse

import numpy as np
import pandas as pd
from intervaltree import Interval, IntervalTree


def precompute_coverage_array(repeats_tree, chrom, region_start, region_end):
    """
    Build prefix sum array for O(1) coverage queries.

    Returns: prefix_sum or None if no repeats.
    """
    if chrom not in repeats_tree:
        return None

    length = region_end - region_start
    if length <= 0:
        return None

    # Binary array: 1 where repeats cover
    coverage = np.zeros(length, dtype=np.uint8)
    overlaps = repeats_tree[chrom].overlap(region_start, region_end)
    for iv in overlaps:
        start = max(iv.begin, region_start) - region_start
        end = min(iv.end, region_end) - region_start
        coverage[start:end] = 1

    # Prefix sums for O(1) range queries
    prefix_sum = np.zeros(length + 1, dtype=np.int64)
    prefix_sum[1:] = np.cumsum(coverage)

    return prefix_sum


def fast_coverage_from_prefix(prefix_sum, region_offset, window_start, window_end):
    """O(1) coverage query using precomputed prefix sums."""
    if prefix_sum is None:
        return 0.0

    local_start = window_start - region_offset
    local_end = window_end - region_offset

    if local_start < 0:
        local_start = 0
    if local_end > len(prefix_sum) - 1:
        local_end = len(prefix_sum) - 1

    if local_start >= local_end:
        return 0.0

    covered = prefix_sum[local_end] - prefix_sum[local_start]
    window_length = window_end - window_start

    return (covered / window_length) * 100 if window_length > 0 else 0.0


def generate_random_starts(region_start, max_start, num_samples):
    """Generate random start positions for sampling"""
    return np.random.randint(region_start, max_start + 1, size=num_samples)


def process_sv_pair(chrom, pos1, pos2, genes_tree, repeats_tree, window_size, num_samples=1000):
    """
    Analyze repeat coverage at both breakpoints of structural variation.

    Returns dict with gene, coverage, percentile, and significance for window's start/end.
    """
    # Get windows around each breakpoint
    win1_start, win1_end, reg1, gene1 = get_window_around_breakpoint(
        chrom, pos1, genes_tree, 'start', window_size)
    win2_start, win2_end, reg2, gene2 = get_window_around_breakpoint(
        chrom, pos2, genes_tree, 'end', window_size)

    # Sampling region coordinates for both breakpoints
    sampling_start = min(win1_start, win2_start)
    sampling_end = max(win1_end, win2_end)

    # Precompute coverage for the whole region
    prefix_sum = precompute_coverage_array(repeats_tree, chrom, sampling_start, sampling_end)

    # Real window coverages
    cov1 = fast_coverage_from_prefix(prefix_sum, sampling_start, win1_start, win1_end)
    cov2 = fast_coverage_from_prefix(prefix_sum, sampling_start, win2_start, win2_end)

    results = {}
    for i, (win_start, win_end, cov, reg, gene) in enumerate([
        (win1_start, win1_end, cov1, reg1, gene1),
        (win2_start, win2_end, cov2, reg2, gene2)
    ]):
        suffix = 'start' if i == 0 else 'end'

        win_size = win_end - win_start

        # Calculate percentile
        if win_size <= 0 or sampling_end - sampling_start <= win_size:
            percentile = 50.0
        else:
            # Generate random windows
            max_start = sampling_end - win_size
            random_starts = generate_random_starts(sampling_start, max_start, num_samples)

            # Compute coverage
            random_coverages = np.zeros(num_samples, dtype=np.float64)
            for j in range(num_samples):
                r_start = random_starts[j]
                r_end = r_start + win_size
                random_coverages[j] = fast_coverage_from_prefix(prefix_sum, sampling_start, r_start, r_end)

            percentile = (np.sum(np.array(random_coverages) < cov) / len(random_coverages)) * 100

        results[f'{suffix}_region'] = reg
        results[f'{suffix}_gene'] = gene if gene else 'NA'
        results[f'{suffix}_window_start'] = win_start
        results[f'{suffix}_window_end'] = win_end
        results[f'{suffix}_coverage'] = cov
        results[f'{suffix}_percentile'] = percentile
        results[f'{suffix}_significant'] = percentile > 95

    return results


def build_interval_tree(df):
    """Build mapping (Chromosome -> IntervalTree) for fast overlap queries."""
    trees = {}
    for chrom, group in df.groupby('chrom'):
        tree = IntervalTree()
        intervals = [
            Interval(row['start'], row['end'], row.get('cluster_id'))
            for _, row in group.iterrows()
        ]
        tree.update(intervals)  # Batch insert is faster
        trees[chrom] = tree
    return trees


def get_window_around_breakpoint(chrom, pos, genes_tree, bp_type, window_size=5000):
    """
    Calculate window's coordinates based on gene context.

    If inside gene: from gene boundary to breakpoint.
    If intergenic: to nearest gene or use default window_size.
    """
    if chrom not in genes_tree:
        return pos - window_size, pos + window_size, 'no_gene', None

    # Check if breakpoint is inside a gene
    containing_genes = genes_tree[chrom].at(pos)
    if containing_genes:
        gene = list(containing_genes)[0]
        if bp_type == 'start':
            return gene.begin, pos, 'inside_gene', gene.data
        else:
            return pos, gene.end, 'inside_gene', gene.data

    # Intergenic - find nearest genes
    tree = genes_tree[chrom]

    # Nearest gene on the left
    prev_genes = tree.overlap(0, pos)
    prev_gene = None
    if prev_genes:
        prev_gene = max(prev_genes, key=lambda x: x.end)
        if prev_gene.end >= pos:
            prev_gene = None

    # Nearest gene on the right
    next_genes = tree.overlap(pos, pos + window_size * 10)
    next_gene = None
    for g in sorted(next_genes, key=lambda x: x.begin):
        if g.begin > pos:
            next_gene = g
            break

    if bp_type == 'start':
        if prev_gene:
            return prev_gene.end, pos, 'intergenic', None
        else:
            return pos - window_size, pos, 'intergenic', None
    else:
        if next_gene:
            return pos, next_gene.begin, 'intergenic', None
        else:
            return pos, pos + window_size, 'intergenic', None

def load_genes(gene_file):
    """Load BED gene file (chrom, start, end, gene_id)."""
    print(f"Loading genes from: {gene_file}")
    genes = pd.read_csv(gene_file, sep='\t', header=None,
                        names=['chrom', 'start', 'end', 'gene_id'],
                        dtype={0: str, 1: np.int64, 2: np.int64, 3: str})
    print(f"Loaded {len(genes)} genes")
    return genes

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--syri', required=True)
    parser.add_argument('--repeats1', required=True)
    parser.add_argument('--repeats2', required=True)
    parser.add_argument('--genes1', required=True)
    parser.add_argument('--genes2', required=True)
    parser.add_argument('--window', type=int, default=5000)
    parser.add_argument('--output', default='breakpoints_with_repeats.tsv')
    parser.add_argument('--threads', type=int, default=1)
    parser.add_argument('--samples', type=int, default=1000, help='Number of random samples')
    args = parser.parse_args()

    # Load data
    print("Loading data...")
    syri = pd.read_csv(args.syri, sep='\t')
    print(f"Loaded {len(syri)} SVs from SYRI")

    reps1 = pd.read_csv(args.repeats1, sep='\t', header=None,
                        names=['chrom', 'start', 'end', 'cluster_id'],
                        dtype={'chrom': str, 'start': np.int64, 'end': np.int64})
    reps2 = pd.read_csv(args.repeats2, sep='\t', header=None,
                        names=['chrom', 'start', 'end', 'cluster_id'],
                        dtype={'chrom': str, 'start': np.int64, 'end': np.int64})

    genes1 = load_genes(args.genes1)
    genes2 = load_genes(args.genes2)

    # Build interval trees
    print("Building interval trees...")
    genes_tree1 = build_interval_tree(
        genes1[['chrom', 'start', 'end', 'gene_id']].rename(columns={'gene_id': 'cluster_id'}))
    genes_tree2 = build_interval_tree(
        genes2[['chrom', 'start', 'end', 'gene_id']].rename(columns={'gene_id': 'cluster_id'}))

    repeats_tree1 = build_interval_tree(reps1)
    repeats_tree2 = build_interval_tree(reps2)
    print("Trees built")

    # Process SVs
    print(f"Processing {len(syri)} SVs...")
    results = []

    for idx, sv in syri.iterrows():
        row_data = {
            'sv_id': idx,
            'type': sv['type'],
            'ref_chrom': sv['ref_chrom'],
            'ref_start': sv['ref_start'],
            'ref_end': sv['ref_end'],
            'qry_chrom': sv['qry_chrom'],
            'qry_start': sv['qry_start'],
            'qry_end': sv['qry_end'],
            'confirmed_by_eaglec': sv['confirmed_by_eaglec']
        }

        # Process reference genome
        ref_chrom = str(sv['ref_chrom'])
        if ref_chrom in genes_tree1:
            res1 = process_sv_pair(
                ref_chrom,
                int(sv['ref_start']),
                int(sv['ref_end']),
                genes_tree1, repeats_tree1, args.window
            )
            for key, value in res1.items():
                row_data[f'ref_{key}_genome1'] = value

        # Process query genome
        qry_chrom = str(sv['qry_chrom'])
        if qry_chrom in genes_tree2:
            res2 = process_sv_pair(
                qry_chrom,
                int(sv['qry_start']),
                int(sv['qry_end']),
                genes_tree2, repeats_tree2, args.window
            )
            for key, value in res2.items():
                row_data[f'qry_{key}_genome2'] = value

        # Determine ancestral state
        for bp in ['start', 'end']:
            genome1_sig = row_data.get(f'ref_{bp}_significant_genome1', False)
            genome2_sig = row_data.get(f'qry_{bp}_significant_genome2', False)
            genome1_cov = row_data.get(f'ref_{bp}_coverage_genome1', 0)
            genome2_cov = row_data.get(f'qry_{bp}_coverage_genome2', 0)

            if genome1_sig and not genome2_sig and genome1_cov > genome2_cov:
                row_data[f'ancestral_{bp}'] = 'genome1'
            elif not genome1_sig and genome2_sig and genome2_cov > genome1_cov:
                row_data[f'ancestral_{bp}'] = 'genome2'
            else:
                row_data[f'ancestral_{bp}'] = 'ambiguous'

        results.append(row_data)

        if idx % 100 == 0:
            print(f"Processed {idx}/{len(syri)} SVs")

    pd.DataFrame(results).to_csv(args.output, sep='\t', index=False)
    print(f"Done! Saved to {args.output}")


if __name__ == '__main__':
    main()

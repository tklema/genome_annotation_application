import argparse
import pandas as pd
import numpy as np
import re

def extract_gene_id(attr):
    """Extract gene ID from GFF3 (ID=xxx) or GTF (gene_id "xxx") attribute string."""
    match = re.search(r'ID=([^;]+)', attr) or re.search(r'gene_id "([^"]+)"', attr)
    return match.group(1) if match else None


def build_gene_regions(df):
    """Merge multiple items into single gene regions."""
    return df.groupby('gene_id').agg({
        'chrom': 'first',
        'start': 'min',
        'end': 'max'
    }).reset_index()


def load_genes(gene_file):
    """
    Load and normalize gene annotations from GFF3/GTF.

    Returns BED-like DataFrame with columns: chrom, start, end, gene_id
    """
    genes = pd.read_csv(gene_file, sep='\t', comment='#', header=None,
                        dtype={0: str, 3: np.int64, 4: np.int64})
    genes.columns = ['chrom', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attributes']

    # If file has gene features, use them directly
    if genes['feature'].str.contains('gene', na=False).any():
        genes = genes[genes['feature'].str.contains('gene', na=False)]
        genes['gene_id'] = genes['attributes'].apply(extract_gene_id)
    else:
        # Otherwise, infer gene regions from all feature items
        genes['gene_id'] = genes['attributes'].apply(extract_gene_id)
        genes = build_gene_regions(genes)

    # Clean up and sort
    genes = genes[['chrom', 'start', 'end', 'gene_id']]
    genes = genes[genes['start'] < genes['end']]
    genes = genes.sort_values(['chrom', 'start'])
    print(f"  Loaded {len(genes)} genes")
    return genes

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--genes', required=True)
    args = parser.parse_args()

    genes = load_genes(args.genes)

    genes.to_csv("genes.bed", sep='\t', header=False, index=False)


if __name__ == "__main__":
    main()

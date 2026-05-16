import argparse
import pandas as pd
import numpy as np
from intervaltree import Interval, IntervalTree
import re

def extract_gene_id(attr):
    match = re.search(r'ID=([^;]+)', attr) or re.search(r'gene_id "([^"]+)"', attr)
    return match.group(1) if match else None

def build_gene_regions(df):
    return df.groupby('gene_id').agg({
        'chrom': 'first',
        'start': 'min',
        'end': 'max'
    }).reset_index()

def load_genes(gene_file):
    genes = pd.read_csv(gene_file, sep='\t', comment='#', header=None,
                        dtype={0: str, 3: np.int64, 4: np.int64})
    genes.columns = ['chrom', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attributes']

    # if 'gene' in genes['feature'].values:
    if genes['feature'].str.contains('gene', na=False).any():
        #genes = genes[genes['feature'] == 'gene']
        genes = genes[genes['feature'].str.contains('gene', na=False)]
        genes['gene_id'] = genes['attributes'].apply(extract_gene_id)
    else:
        genes['gene_id'] = genes['attributes'].apply(extract_gene_id)
        genes = build_gene_regions(genes)

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

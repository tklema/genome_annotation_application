#!/usr/bin/env python3
import pandas as pd
import sys

def split_breakpoints(tsv_file):
    df = pd.read_csv(tsv_file, sep='\t')

    with open('breakpoints_genome1.bed', 'w') as f:
        for _, row in df.iterrows():
            f.write(f"{row['ref_chrom']}\t{row['ref_start']}\t{row['ref_end']}\t{row['sv_id']}\n")

    with open('breakpoints_genome2.bed', 'w') as f:
        for _, row in df.iterrows():
            f.write(f"{row['qry_chrom']}\t{row['qry_start']}\t{row['qry_end']}\t{row['sv_id']}\n")


if __name__ == "__main__":
    split_breakpoints(sys.argv[1])

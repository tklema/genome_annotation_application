import argparse
import sys
from collections import defaultdict
import pandas as pd

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--syri', required=True)
    parser.add_argument('--eaglec2', required=True)
    args = parser.parse_args()

    PAD = 100000

    SV_TYPES = {
        'INV', 'TRANS', 'DUP', 'INVTR', 'INVDUP'
    }

    eaglec_by_chrom = defaultdict(list)
    with open(args.eaglec2) as f:
        next(f)
        for eaglec_id, line in enumerate(f, 1):
            parts = line.strip().split()
            if len(parts) < 4:
                continue
            chrom, pos1, pos2 = parts[0], int(parts[1]), int(parts[3])
            x, y = sorted([pos1, pos2])
            eaglec_by_chrom[chrom].append((x, y, eaglec_id))

    matches = []
    syri_count = 0
    matched_syri_indices = []

    try:
        df_syri = pd.read_csv(args.syri, sep='\t')

        for syri_id, row in df_syri.iterrows():
            sv_type = str(row['type'])
            if sv_type not in SV_TYPES:
                continue

            syri_count += 1
            chrom = str(row['ref_chrom'])
            x1, y1 = sorted([int(row['ref_start']), int(row['ref_end'])])

            if chrom not in eaglec_by_chrom:
                continue

            for x2, y2, eaglec_id in eaglec_by_chrom[chrom]:
                dist_x = abs(x1 - x2)
                dist_y = abs(y1 - y2)
                if dist_x < PAD and dist_y < PAD:
                    matched_syri_indices.append(syri_id)
                    matches.append([
                        chrom, x1, y1, x2, y2,
                        dist_x, dist_y,
                        syri_id, eaglec_id, sv_type
                    ])
                    break
    except FileNotFoundError:
        print(f"Error: file SyRI not found: {args.syri}")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading SyRI file: {e}")
        sys.exit(1)

    if matched_syri_indices:
        df_syri_matched = df_syri.iloc[matched_syri_indices]
        df_syri_matched.to_csv("syri_matched.tsv", sep='\t', index=False)
        print(f"Filtered SyRI file saved to: syri_matched.tsv")
        print(f"Saved {len(df_syri_matched)} out of {syri_count} ({len(df_syri_matched) / syri_count * 100:.1f}%)")
    else:
        if syri_count > 0:
            pd.DataFrame(columns=df_syri.columns).to_csv("syri_matched.tsv", sep='\t', index=False)
            print(f"Matches not found. Created empty file: syri_matched.tsv")
        else:
            print("There are no SyRI rows to process")

    matches_columns = [
            'chrom', 'syri_x1', 'syri_y1', 'eaglec_x2', 'eaglec_y2',
            'distance_x', 'distance_y', 'syri_id', 'eaglec_id', 'sv_type'
        ]
    if matches:
        pd.DataFrame(matches, columns=matches_columns).to_csv("matches.tsv", sep='\t', index=False)
        print(f"Matches info saved to: matches.tsv")
    else:
        pd.DataFrame(columns=matches_columns).to_csv("matches.tsv", sep='\t', index=False)
        print("Matches not found. Created empty file: matches.tsv")


if __name__ == "__main__":
    main()

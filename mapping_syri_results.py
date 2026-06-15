import pandas as pd
import sys

def load_chromosome_mapping(mapping_file):
    """
    Load mapping between synthetic and original chromosomes.

    Returns dict: synthetic_chrom -> {genome1: {chroms, lengths, cumulative}, genome2: {...}}.
    """
    df = pd.read_csv(mapping_file, sep='\t')
    mapping = {}

    for _, row in df.iterrows():
        synth_chrom = row['synthetic_chrom']

        # Parse genome1: comma-separated lists
        genome1_chroms = row['genome1_original_chroms'].split(',')
        genome1_lengths = list(map(int, str(row['genome1_lengths']).split(',')))

        # Cumulative offsets for coordinate conversion
        genome1_cumulative = [0]
        for length in genome1_lengths:
            genome1_cumulative.append(genome1_cumulative[-1] + length)

        # Parse genome2 similarly
        genome2_chroms = row['genome2_original_chroms'].split(',')
        genome2_lengths = list(map(int, str(row['genome2_lengths']).split(',')))

        genome2_cumulative = [0]
        for length in genome2_lengths:
            genome2_cumulative.append(genome2_cumulative[-1] + length)

        mapping[synth_chrom] = {
            'genome1': {
                'chroms': genome1_chroms,
                'lengths': genome1_lengths,
                'cumulative': genome1_cumulative
            },
            'genome2': {
                'chroms': genome2_chroms,
                'lengths': genome2_lengths,
                'cumulative': genome2_cumulative
            }
        }

    return mapping


def map_position(synth_chrom, position, mapping, genome):
    """
    Convert synthetic coordinate back to original chromosome and position.

    Returns: (original_chrom, original_position, chrom_length) or (None, None, None).
    """
    if synth_chrom not in mapping:
        return None, None, None

    chrom_info = mapping[synth_chrom][genome]

    # Find which original chromosome contains this position
    for i in range(len(chrom_info['cumulative']) - 1):
        if chrom_info['cumulative'][i] <= position < chrom_info['cumulative'][i + 1]:
            original_chrom = chrom_info['chroms'][i]
            original_pos = position - chrom_info['cumulative'][i]
            return original_chrom, original_pos, chrom_info['lengths'][i]

    return None, None, None


def main():
    if len(sys.argv) != 5:
        print("Usage: python map_syri_results.py <syri_output> <chrom_mapping> <original_genome1> <original_genome2>")
        sys.exit(1)

    syri_output = sys.argv[1]
    chrom_mapping = sys.argv[2]
    original_genome1 = sys.argv[3]
    original_genome2 = sys.argv[4]

    print("Loading chromosome mapping...")
    mapping = load_chromosome_mapping(chrom_mapping)

    print("Loading SyRI results...")
    syri_columns = [
        'refchr', 'refStart', 'refEnd', 'refSequence', 'qrySequence',
        'qrychr', 'qryStart', 'qryEnd', 'annotationID', 'parentID', 'type', 'info'
    ]

    try:
        syri_df = pd.read_csv(syri_output, sep='\t', header=None, names=syri_columns)
    except:
        try:
            syri_df = pd.read_csv(syri_output, sep='\s+', header=None, names=syri_columns)
        except Exception as e:
            print(f"Error reading SyRI output: {e}")
            sys.exit(1)

    print(f"Processing {len(syri_df)} SVs...")
    print(f"Sample row: {syri_df.iloc[0].tolist()}")

    SV_TYPES = {
        'INV', 'TRANS', 'DUP', 'INVTR', 'INVDUP'
    }

    filtered_svs = []
    for idx, row in syri_df.iterrows():
        synth_ref_chrom = str(row.get('refchr', ''))
        ref_start_val = row.get('refStart', 0)
        ref_end_val = row.get('refEnd', 0)
        ref_start = int(ref_start_val) if str(ref_start_val).isdigit() else 0
        ref_end = int(ref_end_val) if str(ref_end_val).isdigit() else 0

        sv_type = row.get('type', 'UNK')
        if sv_type not in SV_TYPES:
            continue

        synth_qry_chrom = str(row.get('qrychr', ''))
        qry_start_val = row.get('qryStart', 0)
        qry_end_val = row.get('qryEnd', 0)
        qry_start = int(qry_start_val) if str(qry_start_val).isdigit() else 0
        qry_end = int(qry_end_val) if str(qry_end_val).isdigit() else 0

        # Skip unmapped or invalid entries
        if synth_qry_chrom == '-' or synth_ref_chrom == '-':
            continue

        # Map coordinates back to original genomes
        ref_chrom, ref_start_mapped, ref_chrom_len = map_position(synth_ref_chrom, ref_start, mapping, 'genome1')
        _, ref_end_mapped, _ = map_position(synth_ref_chrom, ref_end, mapping, 'genome1')

        qry_chrom, qry_start_mapped, qry_chrom_len = map_position(synth_qry_chrom, qry_start, mapping, 'genome2')
        _, qry_end_mapped, _ = map_position(synth_qry_chrom, qry_end, mapping, 'genome2')

        # Only keep if all coordinates mapped successfully
        if ref_chrom and ref_start_mapped is not None and ref_end_mapped is not None and qry_start_mapped is not None and qry_end_mapped is not None:
            mapped_sv = {
                'type': sv_type,
                'ref_chrom': ref_chrom,
                'ref_start': int(ref_start_mapped),
                'ref_end': int(ref_end_mapped),
                'qry_chrom': qry_chrom if qry_chrom else 'unknown',
                'qry_start': int(qry_start_mapped),
                'qry_end': int(qry_end_mapped)
            }
            filtered_svs.append(mapped_sv)

    # Write outputs
    if filtered_svs:
        mapped_df = pd.DataFrame(filtered_svs)
        mapped_df.to_csv('mapped_svs.tsv', sep='\t', index=False)
        print(f"Mapped {len(filtered_svs)} SVs to original chromosomes")

        with open('summary.txt', 'w') as f:
            f.write(f"Total SVs: {len(filtered_svs)}\n")
            f.write(f"Original genome1: {original_genome1}\n")
            f.write(f"Original genome2: {original_genome2}\n")
            f.write(f"\nSV type distribution:\n")
            f.write(str(mapped_df['type'].value_counts()))
    else:
        print("No SVs could be mapped")
        pd.DataFrame().to_csv('mapped_svs.tsv', sep='\t')
        with open('summary.txt', 'w') as f:
            f.write(f"No SVs were mapped\n")


if __name__ == "__main__":
    main()

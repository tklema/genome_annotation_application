#!/usr/bin/env python3
import argparse
import pandas as pd
import numpy as np
from intervaltree import Interval, IntervalTree
import random
import re

def generate_random_windows_in_region(region_start, region_end, window_size, num_samples=1000):
    """
    Генерирует случайные окна заданного размера внутри региона
    """
    windows = []
    region_length = region_end - region_start
    
    if region_length <= window_size:
        return [(region_start, region_end)]
    
    max_start = region_end - window_size
    for _ in range(num_samples):
        start = random.randint(region_start, max_start)
        windows.append((start, start + window_size))
    
    return windows

def calculate_percentile(value, distribution):
    """Calculate percentile of value in distribution"""
    return (np.sum(np.array(distribution) < value) / len(distribution)) * 100

def build_interval_tree(df):
    """Build interval tree for fast overlap queries"""
    trees = {}
    for chrom in df['chrom'].unique():
        chrom_df = df[df['chrom'] == chrom]
        tree = IntervalTree()
        for _, row in chrom_df.iterrows():
            tree.add(Interval(row['start'], row['end'], row.get('cluster_id', None)))
        trees[chrom] = tree
    return trees

def get_window_around_breakpoint(chrom, pos, genes_tree, bp_type, window_size=5000):
    """
    Определяем окно для брейкпоинта на основе генов
    Возвращает: (start, end, region_type, gene_id)
    region_type: 'inside_gene', 'intergenic', 'no_gene'
    """
    if chrom not in genes_tree:
        return (pos - window_size, pos + window_size, 'no_gene', None)
    
    containing_genes = genes_tree[chrom].overlap(pos, pos)
    print("count of containing_genes: ", len(containing_genes))
    print()

    if containing_genes:
        gene = containing_genes[0]
        if bp_type == 'start':
            return (gene.begin, pos, 'inside_gene', gene.data)
        else:
            return (pos, gene.end, 'inside_gene', gene.data)
    
    all_genes = sorted(genes_tree[chrom], key=lambda x: x.begin)
    
    prev_gene = None
    next_gene = None
    
    for gene in all_genes:
        if gene.end < pos:
            prev_gene = gene
        elif gene.begin > pos:
            next_gene = gene
            break
    
    if bp_type == 'start':
        if prev_gene:
            return (prev_gene.end, pos, 'intergenic', None)
        else:
            return (pos - window_size, pos, 'intergenic', None)
    else:
        if next_gene:
            return (pos, next_gene.begin, 'intergenic', None)
        else:
            return (pos, pos + window_size, 'intergenic', None)

def calculate_coverage_percentage(window_start, window_end, repeats_tree):
    """
    Считает процент покрытия окна повторами
    """
    if window_start >= window_end:
        return 0.0
    
    overlaps = repeats_tree.overlap(window_start, window_end)
    
    if not overlaps:
        return 0.0
    
    # Обрезаем повторы по границам окна
    intervals = []
    for iv in overlaps:
        start = max(iv.begin, window_start)
        end = min(iv.end, window_end)
        if start < end:
            intervals.append((start, end))
    
    if not intervals:
        return 0.0

    intervals.sort()
    
    merged = []
    current_start, current_end = intervals[0]
    
    for start, end in intervals[1:]:
        if start <= current_end:
            current_end = max(current_end, end)
        else:
            merged.append((current_start, current_end))
            current_start, current_end = start, end
    merged.append((current_start, current_end))
    
    covered_length = sum(end - start for start, end in merged)
    window_length = window_end - window_start
    
    return (covered_length / window_length) * 100

def process_sv_pair(chrom, pos1, pos2, genes_tree, repeats_tree, window_size):
    """
    Обрабатывает пару брейкпоинтов (start и end) для одного генома
    Возвращает словарь с результатами для обоих брейкпоинтов
    """
    results = {}
    
    # 1. Находим окна для обоих брейкпоинтов
    win1_start, win1_end, reg1, gene1 = get_window_around_breakpoint(
        chrom, pos1, genes_tree, 'start', window_size)
    win2_start, win2_end, reg2, gene2 = get_window_around_breakpoint(
        chrom, pos2, genes_tree, 'end', window_size)
    
    # 2. Считаем покрытия
    cov1 = calculate_coverage_percentage(win1_start, win1_end, repeats_tree[chrom])
    cov2 = calculate_coverage_percentage(win2_start, win2_end, repeats_tree[chrom])
    
    # 3. Определяем регион для сэмплирования (от начала первого окна до конца второго)
    sampling_region = (min(win1_start, win2_start), max(win1_end, win2_end))
    
    # 4. Для каждого брейкпоинта генерируем случайные окна и считаем перцентили
    for i, (win_start, win_end, cov, reg, gene) in enumerate([
        (win1_start, win1_end, cov1, reg1, gene1),
        (win2_start, win2_end, cov2, reg2, gene2)
    ]):
        suffix = 'start' if i == 0 else 'end'
        
        window_size = win_end - win_start
        random_windows = generate_random_windows_in_region(
            sampling_region[0], sampling_region[1], window_size)
        
        random_coverages = []
        for r_start, r_end in random_windows:
            rand_cov = calculate_coverage_percentage(r_start, r_end, repeats_tree[chrom])
            random_coverages.append(rand_cov)
        
        percentile = (np.sum(np.array(random_coverages) < cov) / len(random_coverages)) * 100
        
        results[f'{suffix}_region'] = reg
        results[f'{suffix}_gene'] = gene if gene else 'NA'
        results[f'{suffix}_window_start'] = win_start
        results[f'{suffix}_window_end'] = win_end
        results[f'{suffix}_coverage'] = cov
        results[f'{suffix}_percentile'] = percentile
        results[f'{suffix}_significant'] = percentile > 95 or percentile < 5
    
    return results

def extract_gene_id_gff3(attr):
    match = re.search(r'ID=([^;]+)', attr)
    return match.group(1) if match else None

def extract_gene_id_gtf(attr):
    match = re.search(r'gene_id "([^"]+)"', attr)
    return match.group(1) if match else None

def build_gene_regions(df):
    return df.groupby('gene_id').agg({
        'chrom': 'first',
        'start': 'min',
        'end': 'max'
    }).reset_index()

def load_genes(gene_file):
    print(f"Loading genes from: {gene_file}")

    # Загружаем данные
    genes = pd.read_csv(gene_file, sep='\t', comment='#', header=None)
    genes.columns = ['chrom', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attributes']

    # Обрабатываем в зависимости от формата
    if gene_file.endswith('.gtf'):
        print("  Detected GTF format")

        # Извлекаем gene_id
        genes['gene_id'] = genes['attributes'].apply(extract_gene_id_gtf)

        # Удаляем строки без gene_id
        # genes = genes.dropna(subset=['gene_id'])

        # Строим регионы генов
        genes = build_gene_regions(genes)

    elif gene_file.endswith(('.gff', '.gff3')):
        print("  Detected GFF3 format")

        # Берем только записи с типом "gene"
        genes = genes[genes['feature'] == 'gene']

        # Извлекаем gene_id
        genes['gene_id'] = genes['attributes'].apply(extract_gene_id_gff3)

        # Удаляем строки без gene_id
        # genes = genes.dropna(subset=['gene_id'])

        # Для GFF3 у нас уже есть координаты, просто выбираем нужные колонки
        genes = genes[['chrom', 'start', 'end', 'gene_id']]
    else:
        raise ValueError(f"Unknown file format: {gene_file}. File must end with .gtf, .gff, or .gff3")

    genes = genes[genes['start'] < genes['end']]

    print(f"  Loaded {len(genes)} genes")
    return genes

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--syri', required=True)
    parser.add_argument('--repeats1', required=True)
    parser.add_argument('--repeats2', required=True)
    parser.add_argument('--genes1', required=True, help='GFF3/GTF файл с генами для genome1')
    parser.add_argument('--genes2', required=True, help='GFF3/GTF файл с генами для genome2')
    parser.add_argument('--window', type=int, default=5000, help='Размер окна по умолчанию (если нет генов)')
    parser.add_argument('--output', default='breakpoints_with_repeats.tsv')
    args = parser.parse_args()

    syri = pd.read_csv(args.syri, sep='\t')
    print(f"Loaded {len(syri)} SVs from SYRI")

    dtype_spec = {0: str, 1: int, 2: int, 3: str}
    
    reps1 = pd.read_csv(args.repeats1, sep='\t', header=None,
                        names=['chrom', 'start', 'end', 'cluster_id'])
    print(f"Loaded {len(reps1)} repeats from genome1")
    
    reps2 = pd.read_csv(args.repeats2, sep='\t', header=None,
                        names=['chrom', 'start', 'end', 'cluster_id'])
    print(f"Loaded {len(reps2)} repeats from genome2")

    genes1 = load_genes(args.genes1)
    genes2 = load_genes(args.genes2)

    print("Sample gene_ids from genome1 (first 20):")
    print(genes1['gene_id'].head(20).tolist())
    print("\nSample gene_ids from genome2 (first 20):")
    print(genes2['gene_id'].head(20).tolist())
    print(f"\nTotal unique gene_ids in genome1: {genes1['gene_id'].nunique()}")
    print(f"Total unique gene_ids in genome2: {genes2['gene_id'].nunique()}")

    print(f"Built {len(genes1)} gene regions from genome1")
    print(f"Built {len(genes2)} gene regions from genome2")

    # Строим interval tree для генов
    genes_tree1 = build_interval_tree(genes1[['chrom', 'start', 'end', 'gene_id']].rename(
        columns={'gene_id': 'cluster_id'}))
    genes_tree2 = build_interval_tree(genes2[['chrom', 'start', 'end', 'gene_id']].rename(
        columns={'gene_id': 'cluster_id'}))
    
    print("Building repeat interval trees...")
    repeats_tree1 = build_interval_tree(reps1)
    repeats_tree2 = build_interval_tree(reps2)
    print("Trees built")

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
            'qry_end': sv['qry_end']
        }
        
        # Для genome1
        if str(sv['ref_chrom']) in genes_tree1:
            res1 = process_sv_pair(
                str(sv['ref_chrom']), 
                int(sv['ref_start']), 
                int(sv['ref_end']),
                genes_tree1, repeats_tree1, args.window
            )
            for key, value in res1.items():
                row_data[f'ref_{key}_genome1'] = value
    
        # Для genome2
        if str(sv['qry_chrom']) in genes_tree2:
            res2 = process_sv_pair(
                str(sv['qry_chrom']),
                int(sv['qry_start']),
                int(sv['qry_end']),
                genes_tree2, repeats_tree2, args.window
            )
            for key, value in res2.items():
                row_data[f'qry_{key}_genome2'] = value
    
        results.append(row_data)
        
        if idx % 100 == 0:
            print(f"Processed {idx}/{len(syri)} SVs")

    pd.DataFrame(results).to_csv(args.output, sep='\t', index=False)
    print(f"Done! Saved to {args.output}")

if __name__ == '__main__':
    main()

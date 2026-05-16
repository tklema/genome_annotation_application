import argparse
import pandas as pd
import numpy as np
from intervaltree import Interval, IntervalTree
import re

def precompute_coverage_array(repeats_tree, chrom, region_start, region_end):
    if chrom not in repeats_tree:
        return None, region_start
    
    length = region_end - region_start
    if length <= 0:
        return None, region_start
    
    # Создаем массив покрытия
    coverage = np.zeros(length, dtype=np.uint8)
    
    overlaps = repeats_tree[chrom].overlap(region_start, region_end)
    for iv in overlaps:
        start = max(iv.begin, region_start) - region_start
        end = min(iv.end, region_end) - region_start
        coverage[start:end] = 1
    
    # Префиксные суммы для O(1) запросов
    prefix_sum = np.zeros(length + 1, dtype=np.int64)
    prefix_sum[1:] = np.cumsum(coverage)
    
    return prefix_sum, region_start


def fast_coverage_from_prefix(prefix_sum, region_offset, window_start, window_end):
    if prefix_sum is None:
        return 0.0
    
    # Преобразуем в локальные координаты
    local_start = window_start - region_offset
    local_end = window_end - region_offset
    
    # Проверяем границы
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
    return np.random.randint(region_start, max_start + 1, size=num_samples)


def calculate_coverage_percentage(window_start, window_end, repeats_tree, chrom):
    if window_start >= window_end:
        return 0.0
    
    if chrom not in repeats_tree:
        return 0.0
    
    overlaps = repeats_tree[chrom].overlap(window_start, window_end)
    
    if not overlaps:
        return 0.0
    
    # Собираем интервалы в numpy массив для быстрой обработки
    intervals = []
    for iv in overlaps:
        start = max(iv.begin, window_start)
        end = min(iv.end, window_end)
        if start < end:
            intervals.append((start, end))
    
    if not intervals:
        return 0.0
    
    # Быстрое слияние интервалов
    intervals = np.array(intervals, dtype=np.int64)
    intervals = intervals[intervals[:, 0].argsort()]
    
    covered_length = merge_intervals(intervals)
    window_length = window_end - window_start
    
    return (covered_length / window_length) * 100

def merge_intervals(intervals):
    if len(intervals) == 0:
        return 0
    
    total_covered = 0
    current_start = intervals[0, 0]
    current_end = intervals[0, 1]
    
    for i in range(1, len(intervals)):
        start = intervals[i, 0]
        end = intervals[i, 1]
        
        if start <= current_end:
            if end > current_end:
                current_end = end
        else:
            total_covered += current_end - current_start
            current_start = start
            current_end = end
    
    total_covered += current_end - current_start
    return total_covered


def process_sv_pair(chrom, pos1, pos2, genes_tree, repeats_tree, window_size, num_samples=1000):
    results = {}
    
    # Находим окна для обоих брейкпоинтов
    win1_start, win1_end, reg1, gene1 = get_window_around_breakpoint(
        chrom, pos1, genes_tree, 'start', window_size)
    win2_start, win2_end, reg2, gene2 = get_window_around_breakpoint(
        chrom, pos2, genes_tree, 'end', window_size)
    
    # Определяем общий регион для сэмплирования
    sampling_start = min(win1_start, win2_start)
    sampling_end = max(win1_end, win2_end)
    
    # Предвычисляем покрытие для всего региона
    prefix_sum, region_offset = precompute_coverage_array(
        repeats_tree, chrom, sampling_start, sampling_end)
    
    # Считаем покрытия для реальных окон
    cov1 = fast_coverage_from_prefix(prefix_sum, region_offset, win1_start, win1_end)
    cov2 = fast_coverage_from_prefix(prefix_sum, region_offset, win2_start, win2_end)
    
    # Обрабатываем оба брейкпоинта
    for i, (win_start, win_end, cov, reg, gene) in enumerate([
        (win1_start, win1_end, cov1, reg1, gene1),
        (win2_start, win2_end, cov2, reg2, gene2)
    ]):
        suffix = 'start' if i == 0 else 'end'
        
        win_size = win_end - win_start
        
        if win_size <= 0 or sampling_end - sampling_start <= win_size:
            percentile = 50.0
            random_coverages = np.array([cov])
        else:
            # Генерируем случайные окна
            max_start = sampling_end - win_size
            random_starts = generate_random_starts(sampling_start, max_start, num_samples)
            
            # Векторизованный расчет покрытий
            random_coverages = np.zeros(num_samples, dtype=np.float64)
            for j in range(num_samples):
                r_start = random_starts[j]
                r_end = r_start + win_size
                random_coverages[j] = fast_coverage_from_prefix(
                    prefix_sum, region_offset, r_start, r_end)
            
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
    trees = {}
    
    for chrom, group in df.groupby('chrom'):
        tree = IntervalTree()
        # Batch добавление
        intervals = [
            Interval(row['start'], row['end'], row.get('cluster_id'))
            for _, row in group.iterrows()
        ]
        tree.update(intervals)
        trees[chrom] = tree
    
    return trees


def get_window_around_breakpoint(chrom, pos, genes_tree, bp_type, window_size=5000):
    if chrom not in genes_tree:
        return (pos - window_size, pos + window_size, 'no_gene', None)
    
    containing_genes = genes_tree[chrom].at(pos)
    
    if containing_genes:
        gene = list(containing_genes)[0]  # Преобразуем в list для доступа
        if bp_type == 'start':
            return (gene.begin, pos, 'inside_gene', gene.data)
        else:
            return (pos, gene.end, 'inside_gene', gene.data)
    
    # Используем envelop для быстрого поиска соседних генов
    tree = genes_tree[chrom]
    
    # Ищем ближайший ген слева
    prev_genes = tree.overlap(0, pos)
    prev_gene = None
    if prev_genes:
        prev_gene = max(prev_genes, key=lambda x: x.end)
        if prev_gene.end >= pos:
            prev_gene = None
    
    # Ищем ближайший ген справа
    next_genes = tree.overlap(pos, pos + window_size * 10)
    next_gene = None
    for g in sorted(next_genes, key=lambda x: x.begin):
        if g.begin > pos:
            next_gene = g
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


#def load_genes(gene_file):
#    print(f"Loading genes from: {gene_file}")
#    
#    genes = pd.read_csv(gene_file, sep='\t', comment='#', header=None,
#                        dtype={0: str, 3: np.int64, 4: np.int64})
#    genes.columns = ['chrom', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attributes']
#    
#    if gene_file.endswith('.gtf'):
#        print("  Detected GTF format")
#        genes['gene_id'] = genes['attributes'].apply(extract_gene_id_gtf)
#        genes = build_gene_regions(genes)
#    elif gene_file.endswith(('.gff', '.gff3')):
#        print("  Detected GFF3 format")
#        genes = genes[genes['feature'] == 'gene']
#        genes['gene_id'] = genes['attributes'].apply(extract_gene_id_gff3)
#        genes = genes[['chrom', 'start', 'end', 'gene_id']]
#    else:
#        raise ValueError(f"Unknown file format: {gene_file}")
#    
#    genes = genes[genes['start'] < genes['end']]
#    print(f"  Loaded {len(genes)} genes")
#    return genes

def load_genes(gene_file):
    print(f"Loading genes from: {gene_file}")
    
    genes = pd.read_csv(gene_file, sep='\t', header=None,
                        names=['chrom', 'start', 'end', 'gene_id'],
                        dtype={0: str, 1: np.int64, 2: np.int64, 3: str})
    print(f"  Loaded {len(genes)} genes")
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
    
    # Загрузка данных
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
    
    # Построение деревьев
    print("Building interval trees...")
    genes_tree1 = build_interval_tree(
        genes1[['chrom', 'start', 'end', 'gene_id']].rename(columns={'gene_id': 'cluster_id'}))
    genes_tree2 = build_interval_tree(
        genes2[['chrom', 'start', 'end', 'gene_id']].rename(columns={'gene_id': 'cluster_id'}))
    
    repeats_tree1 = build_interval_tree(reps1)
    repeats_tree2 = build_interval_tree(reps2)
    print("Trees built")
    
    # Обработка
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

        breakpoints = ['start', 'end']
        for bp in breakpoints:
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

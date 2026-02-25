#!/usr/bin/env python3
import argparse
import pandas as pd
import numpy as np

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--annotation', required=True)
    parser.add_argument('--syri', required=True)
    parser.add_argument('--output', default='directional_svs.tsv')
    args = parser.parse_args()

    # Загружаем аннотацию с брейкпоинтами (уже содержит перцентили)
    ann = pd.read_csv(args.annotation, sep='\t')
    
    # Загружаем оригинальный SYRI
    syri = pd.read_csv(args.syri, sep='\t')

    results = []
    
    for idx, row in ann.iterrows():
        base = {
            'sv_id': idx,
            'type': row['type'],
            'ref_chrom': row['ref_chrom'],
            'ref_start': row['ref_start'],
            'ref_end': row['ref_end'],
            'qry_chrom': row['qry_chrom'],
            'qry_start': row['qry_start'],
            'qry_end': row['qry_end']
        }
        
        # Копируем все поля для каждого брейкпоинта
        for bp in ['ref_start', 'ref_end', 'qry_start', 'qry_end']:
            for genome in ['genome1', 'genome2']:
                for field in ['coverage', 'region', 'gene', 'percentile', 'significant']:
                    col = f'{bp}_{genome}_{field}'
                    if col in row:
                        base[col] = row[col]
        
        # Определяем направление на основе значимости
        sig1_count = 0
        sig2_count = 0
        
        for bp in ['ref_start', 'ref_end', 'qry_start', 'qry_end']:
            if row.get(f'{bp}_genome1_significant', False):
                sig1_count += 1
            if row.get(f'{bp}_genome2_significant', False):
                sig2_count += 1
        
        if sig1_count > 0 and sig2_count == 0:
            direction = 'ancestral_in_genome2'
        elif sig1_count == 0 and sig2_count > 0:
            direction = 'ancestral_in_genome1'
        else:
            direction = 'ambiguous'
        
        base['direction'] = direction
        results.append(base)
    
    # Сохраняем результат
    pd.DataFrame(results).to_csv(args.output, sep='\t', index=False)
    print(f"Saved to {args.output}")

if __name__ == '__main__':
    main()

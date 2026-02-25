import pandas as pd
import numpy as np
from collections import defaultdict
import sys
import argparse

def main():
    # Парсинг аргументов командной строки
    parser = argparse.ArgumentParser(description='Сравнение результатов SyRI и EagleC2')
    parser.add_argument('--syri', required=True, help='Файл с результатами SyRI')
    parser.add_argument('--eaglec2', required=True, help='Файл с результатами EagleC2')
    parser.add_argument('--output', default='matches.tsv', help='Выходной файл с совпадениями')
    parser.add_argument('--pad', type=int, default=100000, help='Допустимое отклонение (в bp)')
    parser.add_argument('--sample', default='sample', help='Имя образца для отчета')
    
    args = parser.parse_args()
    
    PAD = args.pad
    
    SV_TYPES = {
        'INV', 'TRANS', 'DUP', 'INVTR', 'INVDUP'
    }
    
    # 1. Чтение EagleC2 результатов
    eaglec_data = []
    try:
        with open(args.eaglec2) as f:
            next(f)  # пропускаем заголовок, если есть
            for idx, line in enumerate(f, 1):
                parts = line.strip().split()
                if len(parts) < 4:
                    continue
                chrom1, pos1, chrom2, pos2 = parts[0], int(parts[1]), parts[2], int(parts[3])
                x2, y2 = sorted([pos1, pos2])
                eaglec_data.append((chrom1, x2, y2, f'eaglec_{idx}'))
    except FileNotFoundError:
        print(f"Ошибка: файл EagleC2 не найден: {args.eaglec2}")
        sys.exit(1)
    
    # 2. Группировка EagleC2 по хромосомам
    eaglec_by_chrom = defaultdict(list)
    for chrom, x, y, eid in eaglec_data:
        eaglec_by_chrom[chrom].append((x, y, eid))
    
    # 3. Чтение SyRI результатов и поиск совпадений
    matches = []
    syri_count = 0

    try:
        # Read mapped_svs.tsv (tab-separated with header)
        df = pd.read_csv(args.syri, sep='\t')
    
        for idx, row in df.iterrows():
            # Get columns from mapped file
            sv_type = str(row['type'])  # 'SYN1', 'HDR3', etc.
            if sv_type not in SV_TYPES:
                continue
            
            syri_count += 1
            # chrom = str(row['ref_chrom'])  # Original chromosome name
            chrom = str(row['qry_chrom'])
            x1, y1 = sorted([int(row['ref_start']), int(row['ref_end'])])
        
            if chrom not in eaglec_by_chrom:
                continue
            
            for x2, y2, eid in eaglec_by_chrom[chrom]:
                dist_x = abs(x1 - x2)
                dist_y = abs(y1 - y2)
            
                if dist_x < PAD and dist_y < PAD:
                    matches.append([
                        chrom, x1, y1, x2, y2,
                        dist_x, dist_y,
                        f'syri_{idx+1}', eid, sv_type
                    ])
                
    except FileNotFoundError:
        print(f"Ошибка: файл SyRI не найден: {args.syri}")
        sys.exit(1)
    except Exception as e:
        print(f"Ошибка чтения файла SyRI: {e}")
        sys.exit(1)
    
    # 4. Сохранение результатов
    if matches:
        df = pd.DataFrame(matches, columns=[
            'chrom', 'syri_x1', 'syri_y1', 'eaglec_x2', 'eaglec_y2',
            'distance_x', 'distance_y', 'syri_id', 'eaglec_id', 'sv_type'
        ])
        df.to_csv(args.output, sep='\t', index=False)
        print(f"Совпадения сохранены в: {args.output}")
    else:
        df = pd.DataFrame()
        print("Совпадений не найдено")
    
    # 5. Вывод статистики
    print(f"\n=== РЕЗУЛЬТАТЫ ДЛЯ {args.sample.upper()} ===")
    print(f"Всего SyRI SV: {syri_count}")
    print(f"Всего EagleC2 SV: {len(eaglec_data)}")
    print(f"Найдено совпадений: {len(matches)}")
    
    if len(matches) > 0:
        print(f"Уникальных SyRI с совпадениями: {df['syri_id'].nunique()}")
        print(f"Уникальных EagleC2 с совпадениями: {df['eaglec_id'].nunique()}")
    
        print("\nТоп-10 совпадений по точности:")
        df_sorted = df.sort_values(['distance_x', 'distance_y']).head(10)
        for _, row in df_sorted.iterrows():
            print(f"{row['syri_id']} ↔ {row['eaglec_id']}: "
                  f"Δx={row['distance_x']:,} bp, Δy={row['distance_y']:,} bp")
    else:
        with open(args.output, 'w') as f:
            f.write('chrom\tsyri_x1\tsyri_y1\teaglec_x2\teaglec_y2\tdistance_x\tdistance_y\tsyri_id\teaglec_id\tsv_type\n')
    # Возвращаем код выхода для Nextflow
    return 0

if __name__ == "__main__":
    sys.exit(main())

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
    parser.add_argument('--output-syri', default='syri_matched.tsv', help='Выходной файл SyRI только с совпавшими SV')
    parser.add_argument('--output-matches', default='matches.tsv', help='Выходной файл с деталями совпадений')
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
    matched_syri_indices = set()  # множество индексов syri, которые совпали

    try:
        # Читаем mapped_svs.tsv (tab-separated with header)
        df_syri = pd.read_csv(args.syri, sep='\t')

        for idx, row in df_syri.iterrows():
            sv_type = str(row['type'])
            if sv_type not in SV_TYPES:
                continue

            syri_count += 1
            chrom = str(row['ref_chrom'])
            x1, y1 = sorted([int(row['ref_start']), int(row['ref_end'])])

            if chrom not in eaglec_by_chrom:
                continue

            for x2, y2, eid in eaglec_by_chrom[chrom]:
                dist_x = abs(x1 - x2)
                dist_y = abs(y1 - y2)

                if dist_x < PAD and dist_y < PAD:
                    # Запоминаем, что этот индекс syri совпал
                    matched_syri_indices.add(idx)
                    matches.append([
                        chrom, x1, y1, x2, y2,
                        dist_x, dist_y,
                        f'syri_{idx+1}', eid, sv_type
                    ])
                    # Можно break, если нужно только первое совпадение, 
                    # но оставим все возможные совпадения

    except FileNotFoundError:
        print(f"Ошибка: файл SyRI не найден: {args.syri}")
        sys.exit(1)
    except Exception as e:
        print(f"Ошибка чтения файла SyRI: {e}")
        sys.exit(1)

    # 4. Сохраняем отфильтрованный SyRI файл
    if matched_syri_indices:
        df_syri_matched = df_syri.iloc[list(matched_syri_indices)]
        df_syri_matched.to_csv(args.output_syri, sep='\t', index=False)
        print(f"Отфильтрованный SyRI файл сохранен в: {args.output_syri}")
        print(f"Сохранено {len(df_syri_matched)} строк из {syri_count} ({len(df_syri_matched)/syri_count*100:.1f}%)")
    else:
        # Создаем пустой файл с заголовком
        if syri_count > 0:
            # Если были syri строки, но совпадений нет, сохраняем пустой файл с заголовком
            pd.DataFrame(columns=df_syri.columns).to_csv(args.output_syri, sep='\t', index=False)
            print(f"Совпадений не найдено, создан пустой файл: {args.output_syri}")
        else:
            print("Нет SyRI строк для обработки")

    # 5. Сохранение деталей совпадений
    if matches:
        df_matches = pd.DataFrame(matches, columns=[
            'chrom', 'syri_x1', 'syri_y1', 'eaglec_x2', 'eaglec_y2',
            'distance_x', 'distance_y', 'syri_id', 'eaglec_id', 'sv_type'
        ])
        df_matches.to_csv(args.output_matches, sep='\t', index=False)
        print(f"Детали совпадений сохранены в: {args.output_matches}")
    else:
        # Создаем пустой файл с заголовком
        pd.DataFrame(columns=[
            'chrom', 'syri_x1', 'syri_y1', 'eaglec_x2', 'eaglec_y2',
            'distance_x', 'distance_y', 'syri_id', 'eaglec_id', 'sv_type'
        ]).to_csv(args.output_matches, sep='\t', index=False)
        print("Совпадений не найдено, создан пустой файл matches")

    # 6. Вывод статистики
    print(f"\n=== РЕЗУЛЬТАТЫ ДЛЯ {args.sample.upper()} ===")
    print(f"Всего SyRI SV: {syri_count}")
    print(f"Всего EagleC2 SV: {len(eaglec_data)}")
    print(f"Найдено совпадений: {len(matches)}")
    print(f"Уникальных SyRI с совпадениями: {len(matched_syri_indices)}")

if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
import argparse
import os
from collections import defaultdict
from collections import deque
from copy import deepcopy

import matplotlib
import matplotlib.font_manager
import matplotlib.patches as patches
import numpy as np
from matplotlib.patches import Rectangle
from matplotlib.path import Path
from matplotlib.pyplot import get_cmap
from pandas import DataFrame

VARS = ['INV', 'TRANS', 'INVTR', 'DUP', 'INVDP']

TRACK_HEIGHT = 0.6
TOP_TRACK_POS = 7.7
TOP_CHR_POS = 6.7
BOTTOM_CHR_POS = 4.3
BOTTOM_TRACK_POS = 3.3
LEGEND_SPACE = 0.3
CHR_LINE_HEIGHT = 0.1
TRACK_LABEL_SPACE = 0.7
TRACK_SPACING = 0.1

FONT_SIZE = 6
DPI = 300
BBOX = [0, 1.01, 0.5, 0.3]
BBOXMAR = 0.15
BBOX_GENOMES = BBOX.copy()
BBOX_ANNOTATIONS = [BBOX[0] + BBOXMAR, BBOX[1], BBOX[2], BBOX[3]]
BBOX_TRACKS = [BBOX[0] + 2*BBOXMAR, BBOX[1], BBOX[2], BBOX[3]]
BBOX_WINDOW = [BBOX[0] + 3*BBOXMAR, 0.905]

COLOR_INVERSION = '#FFA500'
COLOR_TRANSLOCATION = '#9ACD32'
COLOR_DUPLICATION = '#00BBFF'
COLOR_ENRICHED = 'red'
COLOR_REPEAT_LINE = '#9932CC'
COLOR_REPEAT_BACKGROUND = '#DDA0DD'
COLOR_GENE_LINE = '#2E8B57'
COLOR_GENE_BACKGROUND = '#90EE90'
REFERENCE_GENOME_COLOR = matplotlib.colors.to_hex(get_cmap('tab10')(0))
QUERY_GENOME_COLOR = matplotlib.colors.to_hex(get_cmap('tab10')(1))

def merge_ranges(ranges):
    if len(ranges) < 2:
        return ranges
    for i in ranges:
        if i[0] > i[1]:
            i[1], i[0] = i[0], i[1]
    ranges = ranges[ranges[:, 0].argsort()]
    min_value = ranges[0, 0]
    max_value = ranges[0, 1]
    out_range = deque()
    for i in ranges[1:]:
        if i[0] > max_value:
            out_range.append([min_value, max_value])
            min_value = i[0]
            max_value = i[1]
        elif i[1] > max_value:
            max_value = i[1]
    out_range.append([min_value, max_value])
    return np.array(out_range)

def read_fasta(file_path):
    out = {}
    current_chrom = ''
    current_seq = deque()
    with open(file_path, 'r') as fin:
        for line in fin:
            if '>' in line:
                if current_chrom:
                    out[current_chrom] = ''.join(current_seq)
                    current_seq = deque()
                current_chrom = line.strip().split('>')[1].split(' ')[0]
            else:
                current_seq.append(line.strip())

    if current_chrom:
        out[current_chrom] = ''.join(current_seq)
    return out

def read_sv_tsv(file_path):
    sv_data = []
    with open(file_path, 'r') as fin:
        header = fin.readline().strip().split('\t')

        for line in fin:
            fields = line.strip().split('\t')
            if len(fields) < len(header):
                continue

            fields[2] = fields[2].strip()
            fields[5] = fields[5].strip()
            sv_data.append(fields)

    df = DataFrame(sv_data)
    df.columns = header

    df = df[['ref_chrom', 'ref_start', 'ref_end', 'qry_chrom',
             'qry_start', 'qry_end', 'type',
             'ancestral_start', 'ancestral_end']]

    str_columns = ['ref_chrom', 'qry_chrom', 'type', 'ancestral_start', 'ancestral_end']
    int_columns = ['ref_start', 'ref_end', 'qry_start', 'qry_end']

    df[str_columns] = df[str_columns].astype(str)
    df[int_columns] = df[int_columns].astype(int)

    return df

class BaseTrack:
    def __init__(self, file_path, genome_name, track_name):
        self.file_path = file_path
        self.name = track_name
        self.genome = genome_name
        self.line_color = 'black'
        self.background_color = 'white'
        self.name_color = 'black'
        self.name_size = matplotlib.rcParams['font.size']
        self.name_font = 'Arial'
        self.line_width = 1
        self.background_alpha = 0.7
        self.track_alpha = 1
        self.bin_counts = {}

    def read_data(self, chromosome_lengths):
        raise NotImplementedError

class RepeatTrack(BaseTrack):
    def __init__(self, file_path, genome_name, bin_width=100000):
        super().__init__(file_path, genome_name, "repeats")
        self.line_color = COLOR_REPEAT_LINE
        self.background_color = COLOR_REPEAT_BACKGROUND
        self.bin_width = bin_width

    def read_data(self, chromosome_lengths):
        print(f"Reading repeats: {self.file_path}")
        bin_width = int(self.bin_width)

        genome_idx = 0 if self.genome == "reference" else 1
        available_chroms = set(chromosome_lengths[genome_idx][1].keys())

        bin_counts = defaultdict(deque)
        current_chrom = ''
        positions = deque()
        added_chroms = []

        with open(self.file_path, 'r') as fin:
            for line in fin:
                fields = line.strip().split()
                chrom = fields[0]

                if chrom not in available_chroms:
                    continue

                if not current_chrom:
                    current_chrom = chrom
                    positions.append([int(fields[1]), int(fields[2])])
                elif current_chrom == chrom:
                    positions.append([int(fields[1]), int(fields[2])])
                else:
                    self._process_chromosome(
                        chromosome_lengths, genome_idx, bin_width,
                        bin_counts, added_chroms,
                        current_chrom, positions
                    )
                    current_chrom = chrom
                    positions = deque([[int(fields[1]), int(fields[2])]])

            self._process_chromosome(
                chromosome_lengths, genome_idx, bin_width,
                bin_counts, added_chroms,
                current_chrom, positions
            )

        self.bin_counts = bin_counts

    def _process_chromosome(self, chromosome_lengths, genome_idx, bin_width,
                            bin_counts, added_chroms,
                            chrom, positions):
        ranges = merge_ranges(np.array(positions))
        all_positions = np.array(list(set(
            i for rng in ranges for i in range(rng[0], rng[1])
        )))

        step = bin_width // 2
        bins = np.arange(
            0, chromosome_lengths[genome_idx][1][chrom] + bin_width, step
        )

        bin_values = np.histogram(all_positions, bins)[0]
        bin_counts[chrom] = deque([
            ((bins[i] + bins[i + 1]) / 2, bin_values[i] / bin_width)
            for i in range(len(bin_values))
        ])
        added_chroms.append(chrom)

class GeneTrack(BaseTrack):

    def __init__(self, file_path, genome_name):
        super().__init__(file_path, genome_name, "genes")
        self.line_color = COLOR_GENE_LINE
        self.background_color = COLOR_GENE_BACKGROUND

    def read_data(self, chromosome_lengths):
        print(f"Reading genes: {self.file_path}")

        genome_idx = 0 if self.genome == "reference" else 1
        available_chroms = set(chromosome_lengths[genome_idx][1].keys())

        bin_counts = defaultdict(list)

        with open(self.file_path, 'r') as fin:
            for line in fin:
                fields = line.strip().split()
                chrom = fields[0]

                if chrom not in available_chroms:
                    continue

                start = int(fields[1])
                end = int(fields[2])
                bin_counts[chrom].append((start, end))

        self.bin_counts = bin_counts

class Genome:

    def __init__(self, file_path, name, color):
        self.file_path = file_path
        self.name = name
        self.line_color = color
        self.line_width = 1
        self.lengths = None

    def read_data(self):
        sequences = read_fasta(self.file_path)
        self.lengths = {chrom: len(seq) for chrom, seq in sequences.items()}

class CompositeCoordinateSystem:

    def __init__(self, ref_chromosomes, query_chromosomes, chromosome_lengths):
        self.ref_chromosomes = ref_chromosomes
        self.query_chromosomes = query_chromosomes

        self.ref_offsets, self.ref_total_length = self._compute_offsets(
            ref_chromosomes, chromosome_lengths[0][1]
        )

        self.query_offsets, self.query_total_length = self._compute_offsets(
            query_chromosomes, chromosome_lengths[1][1]
        )

        self.max_length = max(self.ref_total_length, self.query_total_length)

    @staticmethod
    def _compute_offsets(chromosomes, lengths):
        offsets = {}
        offset = 0
        for chrom in chromosomes:
            offsets[chrom] = offset
            offset += lengths[chrom]
        return offsets, offset

    def convert_ref_coord(self, chrom, position):
        return self.ref_offsets.get(chrom, 0) + position

    def convert_query_coord(self, chrom, position):
        return self.query_offsets.get(chrom, 0) + position

def load_tracks(args, chromosome_lengths):
    tracks = []

    if args.ref_repeats:
        track = RepeatTrack(args.ref_repeats, "reference")
        track.read_data(chromosome_lengths)
        tracks.append(track)

    if args.ref_genes:
        track = GeneTrack(args.ref_genes, "reference")
        track.read_data(chromosome_lengths)
        tracks.append(track)

    if args.qry_repeats:
        track = RepeatTrack(args.qry_repeats, "query")
        track.read_data(chromosome_lengths)
        tracks.append(track)

    if args.qry_genes:
        track = GeneTrack(args.qry_genes, "query")
        track.read_data(chromosome_lengths)
        tracks.append(track)

    return tracks

def validate_alignment_to_fasta(alignments, genomes_file):
    output = deque()
    genomes = deque()
    with open(genomes_file, 'r') as fin:
        for i, line in enumerate(fin):
            path, name = line.strip().split("\t")
            color = REFERENCE_GENOME_COLOR if i == 0 else QUERY_GENOME_COLOR
            genome = Genome(path, name, color)
            genome.read_data()
            if i == 0:
                chromosomes = set(np.unique(alignments['ref_chrom']))
            else:
                chromosomes = set(np.unique(alignments['qry_chrom']))
            filtered_lengths = {
                chrom: genome.lengths[chrom]
                for chrom in chromosomes
                if chrom in genome.lengths
            }
            output.append((name, filtered_lengths))
            genomes.append(genome)
    return output, genomes

def filter_input_data(df):
    df = df.loc[((df['ref_end'] - df['ref_start']) >= 10000) |
                ((df['qry_end'] - df['qry_start']) >= 10000)].copy()
    df.sort_values(['qry_chrom', 'qry_start', 'qry_end'], inplace=True)
    df.sort_values(['ref_chrom', 'ref_start', 'ref_end'], inplace=True)
    return df

def fix_inversion_coordinates(df):
    inversion_mask = ['INV' in typ for typ in df['type']]

    df.loc[inversion_mask, 'qry_start'] = df.loc[inversion_mask, 'qry_start'] + df.loc[inversion_mask, 'qry_end']
    df.loc[inversion_mask, 'qry_end'] = df.loc[inversion_mask, 'qry_start'] - df.loc[inversion_mask, 'qry_end']
    df.loc[inversion_mask, 'qry_start'] = df.loc[inversion_mask, 'qry_start'] - df.loc[inversion_mask, 'qry_end']

    return df

def draw_axes(ax, composite):
    bottom_limit = BOTTOM_TRACK_POS - TRACK_HEIGHT - 0.5
    upper_limit = TOP_TRACK_POS + TRACK_HEIGHT + 0.5

    ax.set_ylim(bottom_limit, upper_limit)
    ax.yaxis.set_visible(False)

    ax.set_xlim(0, composite.max_length)
    ax.xaxis.grid(True, which='both', linestyle='--')
    ax.ticklabel_format(axis='x', useOffset=False, style='plain')
    ax.set_axisbelow(True)

    xticks = ax.get_xticks()
    if composite.max_length >= 1_000_000_000:
        scale_factor = 1_000_000_000
        xlabel = 'Chromosome position (in Gbp)'
    elif composite.max_length >= 1_000_000:
        scale_factor = 1_000_000
        xlabel = 'Chromosome position (in Mbp)'
    elif composite.max_length >= 1_000:
        scale_factor = 1_000
        xlabel = 'Chromosome position (in Kbp)'
    else:
        scale_factor = 1
        xlabel = 'Chromosome position'

    xticks_scaled = xticks / scale_factor
    ax.set_xticks(xticks[:-1])
    ax.set_xticklabels(xticks_scaled[:-1])
    ax.set_xlabel(xlabel)

    return ax


def draw_chromosomes(ax, composite, genomes, chromosome_lengths):
    chromosome_lines = []

    ref_genome = genomes[0]
    for i, ref_chrom in enumerate(composite.ref_chromosomes):
        y_pos = TOP_CHR_POS
        start = composite.ref_offsets[ref_chrom]
        end = start + chromosome_lengths[0][1][ref_chrom]

        line = ax.hlines(y_pos, start, end,
                         color=ref_genome.line_color,
                         linewidth=ref_genome.line_width,
                         label=ref_genome.name if i == 0 else "",
                         zorder=2)
        ax.vlines(start, y_pos, y_pos + 0.1, color=ref_genome.line_color, linewidth=1)

        if i == 0:
            chromosome_lines.append(line)

        ax.text(start + 0.01 * chromosome_lengths[0][1][ref_chrom], y_pos + 0.02,
                ref_chrom, fontsize=6, ha='left', va='bottom')

    query_genome = genomes[1]
    for i, qry_chrom in enumerate(composite.query_chromosomes):
        y_pos = BOTTOM_CHR_POS
        start = composite.query_offsets[qry_chrom]
        end = start + chromosome_lengths[1][1][qry_chrom]

        line = ax.hlines(y_pos, start, end,
                         color=query_genome.line_color,
                         linewidth=query_genome.line_width,
                         label=query_genome.name if i == 0 else "",
                         zorder=2)
        ax.vlines(start, y_pos - 0.1, y_pos, color=query_genome.line_color, linewidth=1)

        if i == 0:
            chromosome_lines.append(line)

        ax.text(start + 0.01 * chromosome_lengths[1][1][qry_chrom], y_pos - 0.05,
                qry_chrom, fontsize=6, ha='left', va='top')

    indent_positions = [TOP_CHR_POS, BOTTOM_CHR_POS]
    return ax, indent_positions, chromosome_lines


def get_enriched_vertices(point0, point1, point2, point3):
    x0, y0 = point0
    x1, y1 = point1
    x2, y2 = point2
    x3, y3 = point3

    Q0 = point0
    Q1 = ((x0 + x1) / 2, (y0 + y1) / 2)
    Q2 = ((x0 + 2 * x1 + x2) / 4, (y0 + 2 * y1 + y2) / 4)
    Q3 = ((x0 + 3 * x1 + 3 * x2 + x3) / 8, (y0 + 3 * y1 + 3 * y2 + y3) / 8)

    return [Q0, Q1, Q2, Q3]

def create_bezier_patch(ref_start, ref_end, qry_start, qry_end,
                        ref_y, qry_y, color, alpha,
                        ancestral_start, ancestral_end,
                        label='', line_width=0, zorder=0):
    start_mid = (qry_start - ref_start) / 2
    end_mid = (qry_end - ref_end) / 2
    height_mid = (qry_y - ref_y) / 2

    vertices = [
        (ref_start, ref_y),
        (ref_start, ref_y + height_mid),
        (ref_start + 2 * start_mid, ref_y + height_mid),
        (ref_start + 2 * start_mid, ref_y + 2 * height_mid),
        (qry_end, qry_y),
        (qry_end, qry_y - height_mid),
        (qry_end - 2 * end_mid, qry_y - height_mid),
        (qry_end - 2 * end_mid, qry_y - 2 * height_mid),
        (ref_start, ref_y),
    ]

    codes = [
        Path.MOVETO,
        Path.CURVE4,
        Path.CURVE4,
        Path.CURVE4,
        Path.LINETO,
        Path.CURVE4,
        Path.CURVE4,
        Path.CURVE4,
        Path.CLOSEPOLY,
    ]

    path = Path(vertices, codes)
    main_patch = patches.PathPatch(
        path,
        facecolor=color,
        linewidth=line_width,
        alpha=alpha,
        label=label,
        edgecolor=color,
        zorder=zorder
    )

    patches_list = []
    enriched_codes = [
        Path.MOVETO,
        Path.CURVE4,
        Path.CURVE4,
        Path.CURVE4,
    ]
    enriched_line_width = line_width + 0.1

    enriched_segments = [
        ("genome1", ancestral_start,
         [(ref_start, ref_y),
          (ref_start, ref_y + height_mid),
          (ref_start + 2 * start_mid, ref_y + height_mid),
          (ref_start + 2 * start_mid, ref_y + 2 * height_mid)]),

        ("genome2", ancestral_start,
         [(ref_start + 2 * start_mid, ref_y + 2 * height_mid),
          (ref_start + 2 * start_mid, ref_y + height_mid),
          (ref_start, ref_y + height_mid),
          (ref_start, ref_y)]),

        ("genome2", ancestral_end,
         [(qry_end, qry_y),
          (qry_end, qry_y - height_mid),
          (qry_end - 2 * end_mid, qry_y - height_mid),
          (qry_end - 2 * end_mid, qry_y - 2 * height_mid)]),

        ("genome1", ancestral_end,
         [(qry_end - 2 * end_mid, qry_y - 2 * height_mid),
          (qry_end - 2 * end_mid, qry_y - height_mid),
          (qry_end, qry_y - height_mid),
          (qry_end, qry_y)]),
    ]

    for expected, actual, points in enriched_segments:
        if actual == expected:
            enriched_vertices = get_enriched_vertices(*points)
            enriched_path = Path(enriched_vertices, enriched_codes)
            enriched_patch = patches.PathPatch(
                enriched_path,
                facecolor='none',
                edgecolor=COLOR_ENRICHED,
                linewidth=enriched_line_width,
                alpha=1,
                zorder=zorder + 1
            )
            patches_list.append(enriched_patch)

    patches_list.append(main_patch)
    return patches_list

def draw_structural_variants(ax, alignments, indent_positions):
    alpha = 0.8
    legend_added = {'INV': False, 'TRANS': False, 'DUP': False}
    sv_handles = dict()

    df = deepcopy(alignments)

    df.loc[df['type'] == 'INVTR', 'type'] = 'TRANS'
    df.loc[df['type'] == 'INVDP', 'type'] = 'DUP'

    color_map = {
        'INV': COLOR_INVERSION,
        'TRANS': COLOR_TRANSLOCATION,
        'DUP': COLOR_DUPLICATION
    }
    label_map = {
        'INV': 'Inversion',
        'TRANS': 'Translocation',
        'DUP': 'Duplication'
    }
    df['label'] = [label_map[t] for t in df['type']]
    df.loc[df.duplicated(['label']), 'label'] = ''

    df['color'] = [color_map[t] for t in df['type']]
    df['line_width'] = 0.1
    df['zorder'] = 1
    df['ref_y'] = indent_positions[0]
    df['qry_y'] = indent_positions[1]

    for row in df.itertuples(index=False):
        patches_list = create_bezier_patch(
            row.ref_start, row.ref_end,
            row.qry_start, row.qry_end,
            row.ref_y, row.qry_y,
            row.color, alpha,
            row.ancestral_start, row.ancestral_end,
            label=row.label, line_width=row.line_width, zorder=row.zorder)

        for patch in patches_list:
            ax.add_patch(patch)

        if row.label and not legend_added[row.type]:
            sv_handles[row.type] = patches_list[-1]
            legend_added[row.type] = True

    legend_items = [sv_handles[t] for t in ['INV', 'TRANS', 'DUP'] if t in sv_handles]
    return ax, legend_items

def draw_track_data(ax, track, chromosome, y_base, offset, margin):
    if chromosome not in track.bin_counts:
        return

    if track.name == 'genes':
        for start, end in track.bin_counts[chromosome]:
            x_start = start + offset
            x_end = end + offset
            height = TRACK_HEIGHT * 0.5
            y_center = y_base + (TRACK_HEIGHT - height) / 2
            rect = Rectangle(
                (x_start, y_center),
                x_end - x_start,
                height,
                facecolor=track.line_color,
                alpha=0.7,
                linewidth=0,
                zorder=2
            )
            ax.add_patch(rect)
    else:
        positions = [k[0] for k in track.bin_counts[chromosome]]
        densities = [k[1] for k in track.bin_counts[chromosome]]

        if not positions:
            return

        max_density = max(densities)
        style = {
            'color': track.line_color,
            'linewidth': track.line_width,
            'zorder': 2,
            'alpha': track.track_alpha
        }

        shifted_positions = [pos + offset for pos in positions]
        y_positions = [(d * TRACK_HEIGHT / max_density) + y_base for d in densities]

        ax.fill_between(shifted_positions, y_positions, y_base, **style)

        ax.text(-margin * 2, y_base, '0%', fontsize=3.5, ha='right', va='bottom', color='black')
        ax.text(-margin * 2, y_base + TRACK_HEIGHT, '100%', fontsize=3.5, ha='right', va='top', color='black')

def draw_tracks(ax, tracks, composite):
    margin = composite.max_length / 300

    tracks_up = 0
    tracks_down = 0
    legend_added = {'repeats': False, 'genes': False}
    track_handles = {}

    for track in tracks:
        if track.genome == "reference":
            y_base = TOP_TRACK_POS - TRACK_HEIGHT + tracks_up * (TRACK_SPACING + TRACK_HEIGHT)
            tracks_up += 1
            ax.add_patch(Rectangle(
                (0, y_base),
                composite.ref_total_length,
                TRACK_HEIGHT,
                linewidth=0,
                facecolor=track.background_color,
                alpha=track.background_alpha,
                zorder=1
            ))

            for ref_chrom in composite.ref_chromosomes:
                if ref_chrom in track.bin_counts:
                    offset = composite.ref_offsets[ref_chrom]
                    draw_track_data(ax, track, ref_chrom, y_base, offset, margin)
        else:
            y_base = BOTTOM_TRACK_POS - tracks_down * (TRACK_SPACING + TRACK_HEIGHT)
            tracks_down += 1
            ax.add_patch(Rectangle(
                (0, y_base),
                composite.query_total_length,
                TRACK_HEIGHT,
                linewidth=0,
                facecolor=track.background_color,
                alpha=track.background_alpha,
                zorder=1
            ))

            for qry_chrom in composite.query_chromosomes:
                if qry_chrom in track.bin_counts:
                    offset = composite.query_offsets[qry_chrom]
                    draw_track_data(ax, track, qry_chrom, y_base, offset, margin)

        track_type = track.name
        if not legend_added[track_type]:
            if track_type == 'repeats':
                window_size_kbp = tracks[0].bin_width // 1000
                step_size_kbp = tracks[0].bin_width // 2000
                label = f'repeats (window size = {window_size_kbp}Kbp, step = {step_size_kbp}Kbp)'
            else:
                label = track_type

            handle = Rectangle(
                (0, 0),1, 1,
                facecolor=track.line_color,
                alpha=track.track_alpha,
                label=label
            )
            track_handles[track_type] = handle
            legend_added[track_type] = True

    legend_items = [track_handles[t] for t in ['repeats', 'genes'] if t in track_handles]
    return ax, legend_items

def read_mapping_file(file_path):
    pairs = []
    with open(file_path, 'r') as f:
        for line in f:
            parts = line.strip().split(':')
            ref_chroms = [c for c in parts[0].split(',') if c]
            query_chroms = [c for c in parts[1].split(',') if c]
            pairs.append((ref_chroms, query_chroms))

    print(f"Loaded {len(pairs)} chromosome groups from config file")
    return pairs

def filter_alignments_for_group(alignments, ref_chroms, query_chroms, composite):
    mask = (alignments['ref_chrom'].isin(ref_chroms)) & (alignments['qry_chrom'].isin(query_chroms))
    group_alignments = alignments[mask].copy()

    if group_alignments.empty:
        return group_alignments

    for ref_chrom in ref_chroms:
        ref_mask = group_alignments['ref_chrom'] == ref_chrom
        group_alignments.loc[ref_mask, 'ref_start'] = group_alignments.loc[ref_mask, 'ref_start'].apply(
            lambda x: composite.convert_ref_coord(ref_chrom, x)
        )
        group_alignments.loc[ref_mask, 'ref_end'] = group_alignments.loc[ref_mask, 'ref_end'].apply(
            lambda x: composite.convert_ref_coord(ref_chrom, x)
        )

    for qry_chrom in query_chroms:
        query_mask = group_alignments['qry_chrom'] == qry_chrom
        group_alignments.loc[query_mask, 'qry_start'] = group_alignments.loc[query_mask, 'qry_start'].apply(
            lambda x: composite.convert_query_coord(qry_chrom, x)
        )
        group_alignments.loc[query_mask, 'qry_end'] = group_alignments.loc[query_mask, 'qry_end'].apply(
            lambda x: composite.convert_query_coord(qry_chrom, x)
        )

    return group_alignments

def visualizer(args):
    print('Start visualizing')
    matplotlib.use('agg')

    alignments = read_sv_tsv(args.sr.name)

    mapping_pairs = read_mapping_file(args.mapping.name)

    chromosome_lengths, genomes = validate_alignment_to_fasta(
        alignments, args.genomes.name
    )

    alignments = filter_input_data(alignments)
    alignments = fix_inversion_coordinates(alignments)

    plt = matplotlib.pyplot
    plt.rcParams['font.size'] = FONT_SIZE

    tracks = load_tracks(args, chromosome_lengths)

    for group_idx, (ref_chroms, query_chroms) in enumerate(mapping_pairs, 1):
        print(f"Processing group: {','.join(ref_chroms)} -> {','.join(query_chroms)}")

        composite = CompositeCoordinateSystem(
            ref_chroms, query_chroms, chromosome_lengths
        )

        group_alignments = filter_alignments_for_group(
            alignments, ref_chroms, query_chroms, composite
        )

        fig = plt.figure(figsize=[10, 3])
        ax = fig.add_subplot(111, frameon=False)

        ax = draw_axes(ax, composite)

        ax, indent_positions, chromosome_handles = draw_chromosomes(
            ax, composite, genomes, chromosome_lengths
        )

        legend1 = plt.legend(
            handles=chromosome_handles,
            loc='upper left',
            bbox_to_anchor=BBOX_GENOMES,
            ncol=1,
            borderaxespad=0.,
            frameon=False,
            title='Genomes'
        )
        legend1._legend_box.align = "left"
        plt.gca().add_artist(legend1)

        ax, sv_handles = draw_structural_variants(
            ax, group_alignments, indent_positions
        )

        legend2 = plt.legend(
            handles=sv_handles,
            loc='upper left',
            bbox_to_anchor=BBOX_ANNOTATIONS,
            ncol=1,
            mode='expand',
            borderaxespad=0.,
            frameon=False,
            title='Annotations'
        )
        legend2._legend_box.align = "left"
        plt.gca().add_artist(legend2)

        if tracks:
            ax, track_handles = draw_tracks(ax, tracks, composite)
            legend3 = plt.legend(
                handles=track_handles,
                loc='upper left',
                bbox_to_anchor=BBOX_TRACKS,
                ncol=1,
                borderaxespad=0.,
                frameon=False,
                title='Features'
            )
            legend3._legend_box.align = "left"

        output_file = f"group{group_idx}"
        fig.savefig(
            output_file,
            dpi=DPI,
            bbox_inches='tight',
            pad_inches=0.2
        )
        print(f"Plot {output_file} generated")

    print('Finish visualizing')

def main():
    parser = argparse.ArgumentParser(
        "Plotting structural rearrangements between genomes",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    input_group = parser.add_argument_group("Input/Output files")
    input_group.add_argument('--sr', type=argparse.FileType('r'), required=True)
    input_group.add_argument('--genomes', type=argparse.FileType('r'), required=True)
    input_group.add_argument('--mapping', type=argparse.FileType('r'), required=True)

    track_group = parser.add_argument_group("Track files")
    track_group.add_argument('--ref-repeats', type=str)
    track_group.add_argument('--ref-genes', type=str)
    track_group.add_argument('--qry-repeats', type=str)
    track_group.add_argument('--qry-genes', type=str)

    args = parser.parse_args()
    visualizer(args)


if __name__ == '__main__':
    main()

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
import re

from .vcf import infer_source_set_from_info, parse_info_field


UPLOADED_VARIATION_IDX = 0
LOCATION_IDX = 1
FEATURE_IDX = 4
FEATURE_TYPE_IDX = 5
CONSEQUENCE_IDX = 6
PROTEIN_POSITION_IDX = 9
AMINO_ACIDS_IDX = 10
EXTRA_IDX = 13

CONSEQUENCE_PRIORITY = {
    'frameshift_variant': 0,
    'inframe_insertion': 1,
    'inframe_deletion': 1,
    'missense_variant': 2,
    'synonymous_variant': 3,
}


@dataclass(frozen=True)
class VEPStrandnessStats:
    total_vep_rows_read: int
    total_unique_uploaded_variations: int
    rna_edit_uploaded_variations_considered: int
    rna_edit_variations_with_strand_conflict: int
    vep_rows_dropped_by_strandness_control: int


@dataclass(frozen=True)
class ParsedVEPRow:
    raw_line: str
    uploaded_variation: str
    location: str
    feature: str
    feature_type: str
    consequence_terms: tuple[str, ...]
    protein_position: str
    amino_acids: str
    extra_fields: dict
    strand: int | None


def _normalize_vep_fields(raw_line):
    fields = raw_line.rstrip('\n').split('\t')
    if len(fields) <= EXTRA_IDX:
        fields.extend([''] * (EXTRA_IDX + 1 - len(fields)))
    return fields


def _normalize_chromosome(chromosome):
    chrom = re.sub(r'^chr', '', chromosome.strip())
    if re.match(r'^0[0-9]+$', chrom):
        chrom = chrom.lstrip('0')
    return chrom if chrom else '0'


def _location_key(location):
    if ':' not in location:
        return None
    chrom, pos = location.split(':', 1)
    pos = pos.split('-')[0].strip()
    if not pos:
        return None
    return _normalize_chromosome(chrom), pos


def _has_informative_value(value):
    return value not in ('', '-', '.', '?', './.')


def _has_protein_annotation(row):
    return _has_informative_value(row.protein_position) and _has_informative_value(row.amino_acids)


def _aa_change_key(row):
    return row.amino_acids if _has_informative_value(row.amino_acids) else None


def _consequence_rank(row):
    if not row.consequence_terms:
        return len(CONSEQUENCE_PRIORITY) + 1
    return min(CONSEQUENCE_PRIORITY.get(term, len(CONSEQUENCE_PRIORITY) + 1) for term in row.consequence_terms)


def _row_score(row, aa_change_counts):
    aa_change = _aa_change_key(row)
    return (
        _consequence_rank(row),
        0 if _has_protein_annotation(row) else 1,
        1 if 'NMD_transcript_variant' in row.consequence_terms else 0,
        1 if 'FLAGS' in row.extra_fields else 0,
        -aa_change_counts.get(aa_change, 0),
        row.feature,
        row.raw_line,
    )


def parse_vep_row(raw_line):
    fields = _normalize_vep_fields(raw_line)
    extra_fields = parse_info_field(fields[EXTRA_IDX]) if fields[EXTRA_IDX] else {}
    strand = extra_fields.get('STRAND', None)
    try:
        strand = int(strand) if strand is not None else None
    except ValueError:
        strand = None
    return ParsedVEPRow(
        raw_line=raw_line if raw_line.endswith('\n') else raw_line + '\n',
        uploaded_variation=fields[UPLOADED_VARIATION_IDX].strip(),
        location=fields[LOCATION_IDX].strip(),
        feature=fields[FEATURE_IDX].strip(),
        feature_type=fields[FEATURE_TYPE_IDX].strip(),
        consequence_terms=tuple(
            term.strip() for term in fields[CONSEQUENCE_IDX].split(',') if term.strip()
        ),
        protein_position=fields[PROTEIN_POSITION_IDX].strip(),
        amino_acids=fields[AMINO_ACIDS_IDX].strip(),
        extra_fields=extra_fields,
        strand=strand,
    )


def read_vep_file(vep_path):
    header_lines = []
    vep_rows = []
    with open(vep_path) as handle:
        for line in handle:
            if line.startswith('#'):
                header_lines.append(line)
            else:
                vep_rows.append(line)
    return header_lines, vep_rows


def collect_rna_edit_positions(vcf_path):
    rna_edit_positions = set()
    with open(vcf_path) as handle:
        for line in handle:
            if line.startswith('#'):
                continue
            columns = line.rstrip('\n').split('\t')
            if len(columns) < 8:
                continue
            source_set = infer_source_set_from_info(columns[7])
            if source_set == 'RNA_EDIT':
                rna_edit_positions.add((_normalize_chromosome(columns[0]), columns[1].strip()))
    return rna_edit_positions


def choose_keep_strand(rows):
    aa_change_counts = Counter(_aa_change_key(row) for row in rows if _aa_change_key(row) is not None)
    rows_by_strand = defaultdict(list)
    for row in rows:
        rows_by_strand[row.strand].append(row)

    def strand_score(strand):
        strand_rows = rows_by_strand[strand]
        best_row = min(strand_rows, key=lambda row: _row_score(row, aa_change_counts))
        return _row_score(best_row, aa_change_counts) + (-len(strand_rows), strand)

    return min(rows_by_strand, key=strand_score)


def filter_vep_rows_by_rna_edit_strandness(vep_rows, rna_edit_positions):
    parsed_rows = []
    rows_by_uploaded_variation = defaultdict(list)
    uploaded_variation_order = []

    for raw_line in vep_rows:
        if not raw_line.strip():
            continue
        row = parse_vep_row(raw_line)
        parsed_rows.append(row)
        if row.uploaded_variation not in rows_by_uploaded_variation:
            uploaded_variation_order.append(row.uploaded_variation)
        rows_by_uploaded_variation[row.uploaded_variation].append(row)

    keep_strand_by_uploaded_variation = {}
    rna_edit_uploaded_variations_considered = 0
    rna_edit_variations_with_strand_conflict = 0

    for uploaded_variation in uploaded_variation_order:
        rows = rows_by_uploaded_variation[uploaded_variation]
        location_keys = {_location_key(row.location) for row in rows}
        is_rna_edit = any(location_key in rna_edit_positions for location_key in location_keys if location_key is not None)
        if not is_rna_edit:
            continue

        rna_edit_uploaded_variations_considered += 1
        strand_rows = [
            row for row in rows
            if row.feature_type == 'Transcript' and row.strand in (-1, 1)
        ]
        strands = {row.strand for row in strand_rows}
        if len(strands) <= 1:
            continue

        rna_edit_variations_with_strand_conflict += 1
        keep_strand_by_uploaded_variation[uploaded_variation] = choose_keep_strand(strand_rows)

    kept_rows = []
    dropped_rows = 0
    for row in parsed_rows:
        keep_strand = keep_strand_by_uploaded_variation.get(row.uploaded_variation, None)
        if keep_strand is not None and row.feature_type == 'Transcript' and row.strand in (-1, 1) and row.strand != keep_strand:
            dropped_rows += 1
            continue
        kept_rows.append(row.raw_line)

    stats = VEPStrandnessStats(
        total_vep_rows_read=len(parsed_rows),
        total_unique_uploaded_variations=len(uploaded_variation_order),
        rna_edit_uploaded_variations_considered=rna_edit_uploaded_variations_considered,
        rna_edit_variations_with_strand_conflict=rna_edit_variations_with_strand_conflict,
        vep_rows_dropped_by_strandness_control=dropped_rows,
    )
    return kept_rows, stats


def filter_vep_file_by_rna_edit_strandness(input_path, output_path, vcf_path):
    header_lines, vep_rows = read_vep_file(input_path)
    rna_edit_positions = collect_rna_edit_positions(vcf_path)
    kept_rows, stats = filter_vep_rows_by_rna_edit_strandness(vep_rows, rna_edit_positions)
    with open(output_path, 'w') as out:
        for line in header_lines:
            out.write(line)
        for line in kept_rows:
            out.write(line)
    return stats

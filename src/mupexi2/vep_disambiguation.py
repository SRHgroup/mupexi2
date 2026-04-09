from __future__ import annotations

from collections import Counter
from dataclasses import dataclass

from .vcf import parse_info_field


UPLOADED_VARIATION_IDX = 0
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
class TranscriptDisambiguationStats:
    total_vep_rows_read: int
    total_unique_uploaded_variations: int
    mutations_with_multiple_transcript_rows: int
    transcript_rows_dropped_by_disambiguation: int


@dataclass(frozen=True)
class ParsedVEPRow:
    raw_line: str
    uploaded_variation: str
    feature: str
    feature_type: str
    consequence_terms: tuple[str, ...]
    protein_position: str
    amino_acids: str
    extra_fields: dict


def _normalize_vep_fields(raw_line):
    fields = raw_line.rstrip('\n').split('\t')
    if len(fields) <= EXTRA_IDX:
        fields.extend([''] * (EXTRA_IDX + 1 - len(fields)))
    return fields


def parse_vep_row(raw_line):
    fields = _normalize_vep_fields(raw_line)
    consequence_terms = tuple(
        term.strip() for term in fields[CONSEQUENCE_IDX].split(',') if term.strip()
    )
    extra_fields = parse_info_field(fields[EXTRA_IDX]) if fields[EXTRA_IDX] else {}
    return ParsedVEPRow(
        raw_line=raw_line if raw_line.endswith('\n') else raw_line + '\n',
        uploaded_variation=fields[UPLOADED_VARIATION_IDX].strip(),
        feature=fields[FEATURE_IDX].strip(),
        feature_type=fields[FEATURE_TYPE_IDX].strip(),
        consequence_terms=consequence_terms,
        protein_position=fields[PROTEIN_POSITION_IDX].strip(),
        amino_acids=fields[AMINO_ACIDS_IDX].strip(),
        extra_fields=extra_fields,
    )


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


def _selection_key(row, aa_change_counts):
    aa_change = _aa_change_key(row)
    return (
        0 if row.feature_type == 'Transcript' else 1,
        _consequence_rank(row),
        0 if _has_protein_annotation(row) else 1,
        1 if 'NMD_transcript_variant' in row.consequence_terms else 0,
        1 if 'FLAGS' in row.extra_fields else 0,
        -aa_change_counts.get(aa_change, 0),
        row.feature,
        row.raw_line,
    )


def select_best_vep_rows(vep_rows):
    groups = {}
    uploaded_variation_order = []
    total_rows_read = 0

    for raw_line in vep_rows:
        if not raw_line.strip():
            continue
        row = parse_vep_row(raw_line)
        total_rows_read += 1
        if row.uploaded_variation not in groups:
            groups[row.uploaded_variation] = []
            uploaded_variation_order.append(row.uploaded_variation)
        groups[row.uploaded_variation].append(row)

    kept_rows = []
    mutations_with_multiple_transcript_rows = 0
    for uploaded_variation in uploaded_variation_order:
        rows = groups[uploaded_variation]
        transcript_rows = [row for row in rows if row.feature_type == 'Transcript']
        if len(transcript_rows) > 1:
            mutations_with_multiple_transcript_rows += 1

        aa_change_counts = Counter(_aa_change_key(row) for row in rows if _aa_change_key(row) is not None)
        best_row = min(rows, key=lambda row: _selection_key(row, aa_change_counts))
        kept_rows.append(best_row.raw_line)

    stats = TranscriptDisambiguationStats(
        total_vep_rows_read=total_rows_read,
        total_unique_uploaded_variations=len(uploaded_variation_order),
        mutations_with_multiple_transcript_rows=mutations_with_multiple_transcript_rows,
        transcript_rows_dropped_by_disambiguation=total_rows_read - len(kept_rows),
    )
    return kept_rows, stats


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


def filter_vep_file(input_path, output_path):
    header_lines, vep_rows = read_vep_file(input_path)
    kept_rows, stats = select_best_vep_rows(vep_rows)
    with open(output_path, 'w') as out:
        for line in header_lines:
            out.write(line)
        for line in kept_rows:
            out.write(line)
    return stats

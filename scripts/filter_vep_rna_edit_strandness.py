#!/usr/bin/env python3

from __future__ import annotations

import argparse
import os
import sys


def build_parser():
    parser = argparse.ArgumentParser(
        description='Filter a VEP tabular file with RNA_EDIT strandness control.'
    )
    parser.add_argument('input_vep', help='Input VEP file')
    parser.add_argument('--vcf', required=True, help='VEP-compatible VCF used to identify RNA_EDIT variants')
    parser.add_argument('-o', '--output', required=True, help='Output filtered VEP file')
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)

    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    src_path = os.path.join(repo_root, 'src')
    if src_path not in sys.path:
        sys.path.insert(0, src_path)

    from mupexi2.vep_strandness_control import filter_vep_file_by_rna_edit_strandness

    stats = filter_vep_file_by_rna_edit_strandness(args.input_vep, args.output, args.vcf)
    print(
        'RNA_EDIT strandness control: read {} rows across {} Uploaded_variation values; {} RNA_EDIT variants considered; {} had transcript rows on both strands; {} rows dropped'.format(
            stats.total_vep_rows_read,
            stats.total_unique_uploaded_variations,
            stats.rna_edit_uploaded_variations_considered,
            stats.rna_edit_variations_with_strand_conflict,
            stats.vep_rows_dropped_by_strandness_control,
        )
    )
    return 0


if __name__ == '__main__':
    raise SystemExit(main())

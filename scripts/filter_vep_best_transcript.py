#!/usr/bin/env python3

from __future__ import annotations

import argparse
import os
import sys


def build_parser():
    parser = argparse.ArgumentParser(
        description='Filter a VEP tabular file down to one best transcript row per Uploaded_variation.'
    )
    parser.add_argument('input_vep', help='Input VEP file')
    parser.add_argument('-o', '--output', required=True, help='Output VEP file')
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)

    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    src_path = os.path.join(repo_root, 'src')
    if src_path not in sys.path:
        sys.path.insert(0, src_path)

    from mupexi2.vep_disambiguation import filter_vep_file

    stats = filter_vep_file(args.input_vep, args.output)
    print(
        'Filtered VEP rows: read {} rows across {} Uploaded_variation values; {} mutations had >1 transcript row; {} rows dropped'.format(
            stats.total_vep_rows_read,
            stats.total_unique_uploaded_variations,
            stats.mutations_with_multiple_transcript_rows,
            stats.transcript_rows_dropped_by_disambiguation,
        )
    )
    return 0


if __name__ == '__main__':
    raise SystemExit(main())

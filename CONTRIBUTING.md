# Contributing to MuPeXi2

Thanks for helping improve MuPeXi2. This project is used in an academic HPC setting, so changes must be safe, reproducible, and avoid altering outputs unless explicitly intended and tested.

## Goals
- Improve, update, and refactor MuPeXi2 while keeping output stable.
- Keep multi-source neoantigen workflows reproducible on HPC.

## Current Stage
MuPeXi2 currently supports:
- Multi-source merged VCF workflows with `SOURCE_SET` (`SOMATIC`, `GERMLINE`, `RNA_EDIT`).
- Somatic + germline context-aware peptide generation with phased-context logic.
- RNA-editing-aware runs with known-site filtering (`KNOWN_RNAEDIT_DB` by default, configurable).
- Superpeptide extraction with context provenance columns (`n_mutations`, context mutation IDs/types, `Mutation_Origin`, `edit_sig`).
- Fusion peptide extraction from Arriba output (chimeric proteins / junction antigens).
- Reproducibility/reporting improvements (run-level JSON report, documented output schema).
- Optional per-k parallel extraction (`--parallel-k`) with explicit CPU sufficiency checks.

## Current Context
This fork already supports:
- Fusion peptide extraction from Arriba output (gene-gene junction antigens).
- Somatic + germline aware peptide generation using a phased VCF (germline + somatic).
- RNA-edit integration in merged VCF mode. 

- Runs on HPC:
  - Software is provided via `module load ...`.
  - Prefer **no new dependencies** unless absolutely necessary.
- Installations must work **offline** (`--no-deps`) or via `PYTHONPATH`.

## How We Run on Computerome
Preferred (no install):
- `export PYTHONPATH="$REPO/src:$PYTHONPATH"`
- `python -m mupexi2.cli ...`

If installed locally:
- `pip install --user -e . --no-deps`
- Run as `mupexi2 ...`

Wrapper script: `scripts/run_hpc_example.sh` (your smoke test).

## Definition of Done
- `python -m mupexi2.cli -h` works.
- Smoke test `scripts/run_hpc_example.sh` and better `scripts/run_supertest_hpc_template.sh` passed. The latter ensures a combination of multiple parameters still produce the output.
- Output formats (e.g., `.mupexi`) do **not** change unless the change is documented in `docs/CHANGELOG.md`.

## Style / Workflow
- Keep commits small and reviewable.
- Avoid large refactors without tests.
- Do not change algorithmic logic unless a test demonstrates the need or documents the change.
- Prefer clear error messages over silent failures.

## Roadmap (Near-Term)
Planned next extensions:
- Germline inclusion into fusion-derived peptide extraction.
- Abrupt splicing neoantigen peptide extraction (D. Kwok et al 2025, Nature).
- Evaluation and potential support for circular RNA-derived peptides. (L. Koch RNA: translated circular RNAs 2017, Nature genetics)
- Evaluation and potential support for additional antigens from cryptic promoter activation. (N. Shah et al 2023, Nature genetics)
- Evidence-based feature suggestions with literature support are welcome.

## Scientific/Technical Constraint
RNA-based neoantigen discovery faces capture overlap limitations (e.g., WES kit vs poly-A RNA-seq), which complicates:
- excluding normal variation for RNA-editing calls,
- phasing RNA-derived variants with germline reliably.

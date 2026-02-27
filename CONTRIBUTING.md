# Contributing to MuPeXi2

Thanks for helping improve MuPeXi2. This project is used in an academic HPC setting, so changes must be safe, reproducible, and avoid altering outputs unless explicitly intended and tested.

## Goals
- Improve, update, and refactor MuPeXi2 while keeping output stable.
- Add RNA editing support and richer peptide enumeration over time.

## Current Context
This fork already supports:
- Fusion peptide extraction from Arriba output (chimeric proteins / junction antigens).
- Somatic + germline aware peptide generation using a phased VCF (germline + somatic).

## Hard Constraints
- Patient/sensitive data (VCFs etc.) **must not** be committed to git.
- Runs on HPC (Computerome):
  - Software is provided via `module load ...`.
  - Internet may be limited/unavailable.
  - Prefer **no new dependencies** unless absolutely necessary.
- Installations must work **offline** (`--no-deps`) or via `PYTHONPATH`.

## How We Run on Computerome
Preferred (no install):
- `export PYTHONPATH="$REPO/src:$PYTHONPATH"`
- `python -m mupexi2.cli ...`

If installed locally:
- `pip install --user -e . --no-deps`
- Run as `mupexi2 ...`

Wrapper script: `scripts/run_cprome_example.sh` (keep it working).

## Definition of Done
- `python -m mupexi2.cli -h` works.
- Smoke/regression tests pass.
- Output formats (e.g., `.mupexi`) do **not** change unless:
  - a regression test is updated, **and**
  - the change is documented in `docs/CHANGELOG.md`.

## Style / Workflow
- Keep commits small and reviewable.
- Avoid large refactors without tests.
- Do not change algorithmic logic unless a test demonstrates the need or documents the change.
- Prefer clear error messages over silent failures.

## Roadmap (Near-Term)
Add RNA editing support and richer peptide enumeration:
- Introduce an RNA-editing variant input (VCF-like or dedicated format).
- Extend parsing + peptide generation to:
  - extract multiple peptide versions when multiple mutation sources overlap/stack in the same peptide window,
  - label each peptide/row with mutation origin (`somatic`, `germline`, `rna_editing`, and combinations).
- Make changes incrementally and protect existing behavior with tests.

## Scientific/Technical Constraint
RNA-based neoantigen discovery faces capture overlap limitations (e.g., WES kit vs poly-A RNA-seq), which complicates:
- excluding normal variation for RNA-editing calls,
- phasing RNA-derived variants with germline reliably.

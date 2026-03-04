# MuPeXi2 Output Schema

This document defines the main columns written to `*.mupexi` output files.

## SNV Output Columns

Core columns:
- `HLA_allele`
- `Mut_peptide`
- `Norm_peptide`
- `Mismatches`
- `Cancer_Driver_Gene`
- `Expression_Level`
- `Expression_score`
- `Chr`
- `Gene_ID`
- `Gene_Symbol`
- `Transcript_ID`
- `Mutation_Consequence`
- `Mutation_Origin` (`SOMATIC` or `RNA_EDIT`)
- `superpeptide` (`Yes`/`No`)
- `context_mutation_types` (`NA` or colon-joined source letters such as `G:S`)
- `context_mutation_ids` (`NA` or colon-joined VEP-like mutation IDs)
- `n_mutations` (format: `S=<int>:G=<int>:R=<int>`)
- `edit_sig` (RNA edit signature label when present; otherwise `NA`)
- `has_germline` (`Yes`/`No`)
- `Germline_positions` (list of context germline protein positions)

NetMHC EL columns:
- `Mut_MHCrank_EL`
- `Mut_MHCscore_EL`
- `Mutant_affinity_score`
- `Norm_MHCrank_EL`
- `Norm_MHCscore_EL`
- `Normal_affinity_score`

NetMHC BA columns (when BA mode is enabled):
- `Mut_MHCrank_BA`
- `Mut_MHCscore_BA`
- `Mut_MHCaffinity`
- `Norm_MHCrank_BA`
- `Norm_MHCscore_BA`
- `Norm_MHCaffinity`

SNV-specific columns:
- `mutation_id_vep`
- `Amino_Acid_Change`
- `peptide_position` (1-based mutation position within peptide)
- `unique_peptide_id` (`<mutation_id_vep>|<peptide_position>`)
- `tumor_vaf`
- `normal_vaf`
- `t_depth`
- `n_depth`
- `t_alt_count`
- `Genomic_Position`
- `Protein_position`
- `priority_Score`

## Fusion Output Columns

Fusion runs reuse shared core and NetMHC columns, and add:
- `fusion_id`
- `Chr2`
- `breakpoints`
- `Gene_ID2`
- `Transcript_ID2`
- `Gene_Symbol2`
- `junction_coverage`
- `spanning_coverage`
- `left_breakpoint_distance`
- `sites`

## Run Reproducibility Report

Each run also writes a JSON report:
- `<prefix>_run_report.json`

The report captures:
- runtime timestamps and duration
- command line
- Python/platform/package versions
- selected run parameters
- phasing/context counters
- VEP and peptide summary counters (when available)

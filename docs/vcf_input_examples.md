# VCF Input Examples (Synthetic)

These files are **synthetic examples** (not real patient variants) showing the expected VCF structure for MuPeXi2.

## Files

- `docs/vcf_examples/merged_phased_rna_dna_example.vcf`
- `docs/vcf_examples/dna_only_mutect2_example.vcf`

## 1) Merged phased VCF (SOMATIC + GERMLINE + RNA_EDIT)

Use this style when running with options like:
- `--vcf-type merged`
- `--germlines true`
- `--rna-edit true`

Expected sample tracks in this project are:
- `DNA_NORMAL`
- `TUMOR`

`RNA_EDIT` records are distinguished by `INFO/SOURCE_SET=RNA_EDIT`; they do not require a separate `RNA_TUMOR` sample column.

Required/important fields:
- `INFO/SOURCE_SET` is required when `--germlines=true` or `--rna-edit=true`.
  - Allowed values expected by MuPeXi2: `SOMATIC`, `GERMLINE`, `RNA_EDIT`.
- `INFO/KNOWN_RNAEDIT_DB` is required by default for RNA-edit records to pass filtering.
  - Default behavior is equivalent to `--rnaedit-known-only`.
  - If this key is absent on RNA-edit records, use `--rnaedit-allow-novel` to include them.
  - The key name is configurable with `--rnaedit-known-key` (for example if your pipeline uses `KNOWN_DB` instead of `KNOWN_RNAEDIT_DB`).
  - Common values can include `ASAOKA`, `RADAR`, `REDI_PORTAL`, and `APOBEC3_MOTIF`; MuPeXi2 treats all of them as known when this key is present.
- For phasing-aware context application, tumor sample genotype/phasing should use:
  - `GT` with `|` (phased genotype), and
  - `PS` (phase set).

Notes:
- By default MuPeXi2 keeps:
  - `SOMATIC` records with `FILTER=PASS`
  - `GERMLINE` records with `FILTER=PASS` (when `--germlines=true`)
  - `RNA_EDIT` records with `KNOWN_RNAEDIT_DB` present (when `--rna-edit=true` and known-only mode), including `APOBEC3_MOTIF`

## 2) DNA-only (somatic) VCF

Use this style for somatic-only runs, for example:
- `--vcf-type mutect2`
- `--germlines false`
- `--rna-edit false`

In this mode, `SOURCE_SET` is not required.

## Minimal run examples

Merged/phased:

```bash
python3 -m mupexi2.cli \
  -v docs/vcf_examples/merged_phased_rna_dna_example.vcf \
  --vcf-type merged \
  --germlines true \
  --rna-edit true \
  --tumor-sample TUMOR \
  --normal-sample DNA_NORMAL \
  -l 9-11 \
  -a HLA-A26:01,HLA-A24:02,HLA-B27:02,HLA-B15:09,HLA-C02:02,HLA-C07:04
```

DNA-only:

```bash
python3 -m mupexi2.cli \
  -v docs/vcf_examples/dna_only_mutect2_example.vcf \
  --vcf-type mutect2 \
  --germlines false \
  --rna-edit false \
  --tumor-sample TUMOR \
  --normal-sample DNA_NORMAL \
  -l 9-11 \
  -a HLA-A26:01,HLA-A24:02,HLA-B27:02,HLA-B15:09,HLA-C02:02,HLA-C07:04
```

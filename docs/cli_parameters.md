# MuPeXi2 CLI Parameters and Input Files

This page is the canonical reference for MuPeXi2 CLI usage.

## Input files

- `--vcf-file` (`-v`): SNV/indel VCF. Can be plain `.vcf` or gzipped `.vcf.gz`.
- `--fusion-file` (`-z`): fusion TSV (Arriba-like format) for junction/chimeric peptide workflow.
- `--expression-file` (`-e`): tab-separated expression table with ID and mean expression.
  - With `--expression-type transcript`, IDs are transcript-like (for example `ENST...`).
  - With `--expression-type gene`, IDs are gene-like (for example `ENSG...`).
- `--config-file` (`-c`): ini with local paths to references/tools (VEP, NetMHCpan, references).

At least one of `--vcf-file` or `--fusion-file` must be provided.

## Parameter reference

### Core options

- `-v, --vcf-file <file>`: SNV/indel VCF input.
- `-z, --fusion-file <file>`: fusion TSV input.
- `-a, --alleles <list>`: comma-separated HLA alleles.
  - Default: `HLA-A02:01`
- `-l, --length <spec>`: peptide lengths (`9`, `9-11`, `9,10,11`).
  - Default: `9`
- `-e, --expression-file <file>`: expression table path.
  - Default: none
- `-E, --expression-type <type>`: expression ID namespace (`transcript` or `gene`).
  - Default: `transcript`

### Output options

- `-o, --output-file <name>`: output table name.
  - Default: `<prefix>.mupexi`
- `-d, --out-dir <dir>`: output directory.
  - Default: current directory
- `-p, --prefix <name>`: prefix for output artifacts.
  - Default: VCF basename (if VCF is used)
- `-L, --log-file <name>`: run log filename.
  - Default: `<prefix>.log`
- `-f, --make-fasta`: write long-peptide FASTA output.
  - Default: disabled
- `-t, --keep-temp`: keep temporary files.
  - Default: disabled
- `-M, --mismatch-only`: print mismatch-only normal peptide representation.
  - Default: disabled
- `-n, --netmhc-full-anal`: run NetMHCpan in EL+BA mode.
  - Default: EL-only mode

### Runtime/reference options

- `-c, --config-file <file>`: config ini path.
  - Default: package `config.ini`
- `-A, --assembly <name>`: assembly passed to VEP.
  - Default: species default
- `-s, --species <name>`: species (`human`, `mouse`, `mouse_black6`, `mouse_balbc`).
  - Default: `human`
- `-F, --fork <int>`: VEP fork count; must be > 1.
  - Default: `2`
- `-g, --liftover`: liftover hg19 VCF to GRCh38 (requires configured tools).
  - Default: disabled
- `-w, --webface`: webserver mode.
  - Default: disabled

### Multi-source VCF options

- `--vcf-type <type>`: VCF style (`mutect2` or `merged`).
  - Default: `mutect2`
- `--vcf_type <type>`: alias for `--vcf-type`.
- `--germlines <bool>`: include `GERMLINE` as context variants.
  - Default: `false`
- `--rna-edit <bool>`: include `RNA_EDIT` as primary variants.
  - Default: `false`
- `--rna_edit <bool>`: alias for `--rna-edit`.
- `--tumor-sample <name>`: explicit tumor sample column name.
  - Default: auto-detect from VCF header
- `--normal-sample <name>`: explicit normal sample column name.
  - Default: auto-detect from VCF header
- `--rnaedit-known-only`: keep only known RNA-edit records.
  - Default behavior: enabled
- `--rnaedit-allow-novel`: include novel RNA-edit records.
  - Overrides `--rnaedit-known-only`
- `--rnaedit-known-key <key>`: INFO key used to identify known RNA-edit sites.
  - Default: `KNOWN_RNAEDIT_DB`
  - Any value under that key is treated as known, including labels such as `ASAOKA`, `RADAR`, `REDI_PORTAL`, and `APOBEC3_MOTIF`
- `--phasing-mode <mode>`: phasing policy for germline context (`auto`, `strict`).
  - Default: `auto`

### Context and superpeptide options

- `--superpeptides <bool>`: allow phased nearby non-germline context (`SOMATIC`/`RNA_EDIT`) in peptide construction.
  - Default: `false`
- `--parallel-k <bool>`: process each k-length in parallel.
  - Default: `false`
- `--parallel_k <bool>`: alias for `--parallel-k`.

Context QC filters for context variants only:

- `--context-rna-min-depth <n>`
  - Purpose: minimum tumor depth for `RNA_EDIT` context variant to be used.
  - Default: `0`
- `--context-rna-min-alt-count <n>`
  - Purpose: minimum tumor alt count for `RNA_EDIT` context variant to be used.
  - Default: `0`
- `--context-rna-min-vaf <f>`
  - Purpose: minimum tumor VAF for `RNA_EDIT` context variant to be used.
  - Default: `0`
- `--context-somatic-min-depth <n>`
  - Purpose: minimum tumor depth for `SOMATIC` context variant to be used.
  - Default: `0`
- `--context-somatic-min-alt-count <n>`
  - Purpose: minimum tumor alt count for `SOMATIC` context variant to be used.
  - Default: `0`
- `--context-somatic-min-vaf <f>`
  - Purpose: minimum tumor VAF for `SOMATIC` context variant to be used.
  - Default: `0`

Notes:
- Germline context variants are not filtered by these context QC thresholds.
- Thresholds apply only when that variant is considered as context.

## Run artifacts

In addition to `*.mupexi` output and optional FASTA/log/temp files, each run writes:

- `<prefix>_run_report.json`: reproducibility summary (parameters, versions, counters, timing)

Output column definitions are documented in:

- `docs/output_schema.md`

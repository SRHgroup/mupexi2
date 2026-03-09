# mupexi2 (lab repo)

This repository packages the current **MuPeXI2** workflow for safer iteration and reproducible runs.

## Install

```bash
conda create -n mupexi2 python=3.11 -y
conda activate mupexi2
pip install -e ".[dev]"
```

For offline/HPC runs you can also use:

```bash
export PYTHONPATH="$REPO/src:${PYTHONPATH:-}"
python3 -m mupexi2.cli -h
```

## Documentation

- Full CLI/input reference: `docs/cli_parameters.md`
- Synthetic VCF input examples: `docs/vcf_input_examples.md`
- Output column schema: `docs/output_schema.md`
- Contributor workflow: `CONTRIBUTING.md`

## Quick usage

Somatic-only style:

```bash
python3 -m mupexi2.cli \
  -v /path/to/sample.vcf.gz \
  -l 9-11 \
  -a HLA-A26:01,HLA-A24:02,HLA-B27:02,HLA-B15:09,HLA-C02:02,HLA-C07:04 \
  -c /path/to/config.ini \
  -p Sample_01 \
  -d /path/to/outdir \
  -e /path/to/expression.tsv
```

Merged multi-source VCF (SOMATIC + GERMLINE + RNA_EDIT):

```bash
python3 -m mupexi2.cli \
  -v /path/to/merged_phased.vcf.gz \
  --vcf-type merged \
  --germlines true \
  --rna-edit true \
  --tumor-sample TUMOR \
  --normal-sample DNA_NORMAL \
  --phasing-mode auto \
  --superpeptides true \
  --context-rna-min-depth 10 \
  --context-rna-min-alt-count 3 \
  --context-rna-min-vaf 0.05 \
  -l 9-11 \
  -a HLA-A26:01,HLA-A24:02,HLA-B27:02,HLA-B15:09,HLA-C02:02,HLA-C07:04 \
  -c /path/to/config.ini \
  -p Sample_01 \
  -d /path/to/outdir \
  -e /path/to/expression.tsv
```

## Tests

Current tests are lightweight and avoid large reference bundles:
- CLI help exits cleanly
- peptide-length parsing works
- default config path resolves inside the package

```bash
pytest -q
```

## Example artifacts

`examples/` contains:
- `your_server.ini` (example cluster config)
- `output_example.mupexi` (example output for regression/reference)

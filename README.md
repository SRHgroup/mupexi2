# mupexi2 

This repository packages the current **MuPeXI2** workflow for neoantigen peptide extraction and immunogenicity annotation from cancer mutation calling files.

## Install

Clone the repo
```bash
git clone https://github.com/SRHgroup/mupexi2.git
cd mupexi2
```

We reccomend creating dedicated conda enviornment to avoid dependencies issues
```
conda create -n mupexi2 python=3.11 -y
conda activate mupexi2
```

Install
```
pip install -e ".[dev]"
```

For offline/HPC runs you can also use:

```bash
cd /path/to/mupexi2
export PYTHONPATH="$REPO/src:${PYTHONPATH:-}"
python3 -m mupexi2.cli -h
```

## Documentation

- Full input parameters: `docs/cli_parameters.md`
- VCF input example schema: `docs/vcf_input_examples.md`
- VCF example: `docs/vcf_examples`
- Output file schema: `docs/output_schema.md`
- Output file example: `docs/output_example.mupexi`
- CONFIG file example: `docs/your_server.ini`
- Contributor workflow: `CONTRIBUTING.md`

## Quick usage

Somatic-only mode:

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

Merged multi-source VCF mode (SOMATIC + GERMLINE + RNA_EDIT):

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


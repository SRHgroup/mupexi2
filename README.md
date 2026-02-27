# mupexi2 (lab repo)

This repository turns the current **MuPeXI2** script into a small, testable Python package.
The goal is **safe iteration**: every behavioural change should be backed by a regression test.

## Install (recommended: conda)

```bash
conda create -n mupexi2 python=3.11 -y
conda activate mupexi2
pip install -e ".[dev]"
```

## Run

The CLI keeps the legacy MuPeXI options (getopt). Example (adapt paths):

```bash
mupexi2 \
  -v /path/to/sample.vcf \
  -l 9-11 \
  -a HLA-A26:01,HLA-A24:02,HLA-B27:02,HLA-B15:09,HLA-C02:02,HLA-C07:04 \
  -t -f -n \
  -c /path/to/your_server.ini \
  -p Patient_P3 \
  -d /path/to/outdir \
  -e /path/to/expression/*.Rstat.txt
```

## Tests

Current tests are **lightweight** and run without large reference files:
- CLI help exits cleanly
- peptide-length parsing works
- default config path resolves inside the package

Add a full end-to-end regression test once you can share **small** example inputs
(VCF + minimal reference set or a path-based HPC test harness).

```bash
pytest -q
```

## Included example artifacts

`examples/` contains:
- `your_server.ini` (example config from your cluster)
- `output_example.mupexi` (example output for reference / future regression testing)

#!/usr/bin/env bash
set -euo pipefail

# Example HPC runner (adapt modules/paths to your cluster)
module load tools ngs anaconda3/2023.03 netmhcpan/4.0a perl ensembl-tools/87

# After installing mupexi2 (pip install -e .), you can run:
mupexi2 \
  -v /path/to/Patient_P3_3.2.Filtered.som.AF.vcf \
  -l 9-11 \
  -a HLA-A26:01,HLA-A24:02,HLA-B27:02,HLA-B15:09,HLA-C02:02,HLA-C07:04 \
  -t -f -n \
  -c examples/your_server.ini \
  -p Patient_P3 \
  -d /path/to/outdir \
  -e "/path/to/kallisto/Patient_P3_T*_1.4.RunStatBootstrapMean.Rstat.txt"

#!/usr/bin/env bash
set -euo pipefail

# Example HPC runner (adapt modules/paths to your cluster)
module load tools ngs anaconda3/2025.06-1 netmhcpan/4.0a perl/5.36.1 ensembl-tools/90

REPO_DIR="/home/projects/SRHgroup/apps/mupexi2"
export PYTHONPATH="$REPO_DIR/src:${PYTHONPATH:-}"

# After installing mupexi2 (pip install -e .), you can run:
python3 -m mupexi2.cli \
        -v /home/projects/SRHgroup/projects/MuPeXI_germline/data/dna/vcf/Patient_11_3.2.Filtered.som.AF.vcf \
        -l 9-11 \
        -a HLA-A26:01,HLA-A24:02,HLA-B27:02,HLA-B15:09,HLA-C02:02,HLA-C07:04 \
        -t -f -n \
        -c /home/projects/SRHgroup/apps/netMHCpan/4.0a/config.computerome.ini \
        -p Patient_11 \
        -d /home/projects/SRHgroup/projects/MuPeXI_germline/data/mupexi/new_tests \
        -e /home/projects/SRHgroup/projects/MuPeXI_germline/data/rna/kallisto/Patient_11_T*_1.4.RunStatBootstrapMean.Rstat.txt

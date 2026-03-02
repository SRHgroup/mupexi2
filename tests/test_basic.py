import os
import subprocess
import sys

import pytest

from mupexi2.core import extract_peptide_length, read_options


def test_extract_peptide_length_range():
    assert extract_peptide_length("9-11") == [9, 10, 11]


def test_extract_peptide_length_csv():
    assert extract_peptide_length("8,10,12") == [8, 10, 12]


def test_cli_help_exits_zero():
    # Use module invocation so this works in editable installs without relying on console scripts
    p = subprocess.run([sys.executable, "-m", "mupexi2.cli", "-h"], capture_output=True, text=True)
    assert p.returncode == 0
    # usage() prints help to stdout
    assert "MuPeXI" in (p.stdout + p.stderr)


def test_default_config_path_is_packaged(tmp_path):
    # read_options does not touch filesystem. We only check it picks a non-root default.
    opts = read_options(["-v", "dummy.vcf"])
    assert opts.config.endswith("config.ini")
    # Important: should not resolve to /config.ini when running as a console script.
    assert not os.path.isabs(opts.config) or "site-packages" in opts.config or "src" in opts.config


def test_new_vcf_mode_defaults():
    opts = read_options(["-v", "dummy.vcf"])
    assert opts.vcf_type == "mutect2"
    assert opts.germlines is False
    assert opts.rna_edit is False
    assert opts.rnaedit_known_only is True
    assert opts.rnaedit_known_key == "KNOWN_RNAEDIT_DB"


def test_new_vcf_mode_parsing():
    opts = read_options([
        "-v", "dummy.vcf",
        "--vcf-type", "merged",
        "--germlines", "true",
        "--rna-edit", "true",
        "--tumor-sample", "TUMOR",
        "--normal-sample", "DNA_NORMAL",
        "--rnaedit-allow-novel",
        "--rnaedit-known-key", "MY_KEY",
    ])
    assert opts.vcf_type == "merged"
    assert opts.germlines is True
    assert opts.rna_edit is True
    assert opts.tumor_sample == "TUMOR"
    assert opts.normal_sample == "DNA_NORMAL"
    assert opts.rnaedit_known_only is False
    assert opts.rnaedit_known_key == "MY_KEY"

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

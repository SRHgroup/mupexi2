import os
import subprocess
import sys

import pytest

from mupexi2.core import (
    parse_genotype_state,
    should_apply_context_variant,
    create_vep_compatible_vcf,
    extract_peptide_length,
    infer_source_set_from_info,
    read_options,
    resolve_sample_indices,
    vcf_has_sourceset,
)


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


def test_infer_source_set_from_info():
    assert infer_source_set_from_info("AC=1;SOURCE_SET=RNA_EDIT;DP=10") == "RNA_EDIT"
    assert infer_source_set_from_info("AC=1;SOMATIC;DP=10") == "SOMATIC"
    assert infer_source_set_from_info("AC=1;GERMLINE;DP=10") == "GERMLINE"
    assert infer_source_set_from_info("AC=1;DP=10") == "UNKNOWN"


def test_vcf_has_sourceset(tmp_path):
    p = tmp_path / "with_ss.vcf"
    p.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR\tDNA_NORMAL\n"
        "1\t10\t.\tA\tG\t.\tPASS\tSOURCE_SET=RNA_EDIT\tGT\t0/1\t0/0\n"
    )
    assert vcf_has_sourceset(str(p), None) is True

    p2 = tmp_path / "without_ss.vcf"
    p2.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR\tDNA_NORMAL\n"
        "1\t10\t.\tA\tG\t.\tPASS\tDP=3\tGT\t0/1\t0/0\n"
    )
    assert vcf_has_sourceset(str(p2), None) is False


def test_create_vep_compatible_vcf_filters_rna_edit_known_only(tmp_path):
    vcf = tmp_path / "input.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR\tDNA_NORMAL\n"
        "chr1\t10\t.\tA\tG\t.\tPASS\tSOURCE_SET=SOMATIC\tGT:AD:AF\t0/1:10,5:0.33\t0/0:12,0:0.0\n"
        "chr1\t11\t.\tA\tT\t.\tq10\tSOURCE_SET=SOMATIC\tGT:AD:AF\t0/1:10,5:0.33\t0/0:12,0:0.0\n"
        "chr1\t12\t.\tA\tC\t.\tstrand_bias\tSOURCE_SET=RNA_EDIT;KNOWN_RNAEDIT_DB=REDIportal\tGT:AD:AF\t0/1:10,5:0.33\t0/0:12,0:0.0\n"
        "chr1\t13\t.\tA\tC\t.\tstrand_bias\tSOURCE_SET=RNA_EDIT\tGT:AD:AF\t0/1:10,5:0.33\t0/0:12,0:0.0\n"
        "chr1\t14\t.\tA\tC\t.\tPASS\tSOURCE_SET=GERMLINE\tGT:AD:AF\t0/1:10,5:0.33\t0/0:12,0:0.0\n"
    )
    out = create_vep_compatible_vcf(
        str(vcf), None, None, str(tmp_path), "x", str(tmp_path), None,
        germlines=True, rna_edit=True, rnaedit_known_only=True, rnaedit_known_key="KNOWN_RNAEDIT_DB"
    )
    with open(out.name) as fh:
        lines = [ln.strip() for ln in fh if ln.strip() and not ln.startswith("#")]
    assert len(lines) == 3
    assert lines[0].startswith("1\t10")
    assert any("\t12\t" in ln for ln in lines)
    assert any("\t14\t" in ln for ln in lines)


def test_resolve_sample_indices_new_names():
    header = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "DNA_NORMAL", "TUMOR"]
    t_idx, n_idx, t_name, n_name = resolve_sample_indices(header)
    assert t_idx == 10
    assert n_idx == 9
    assert t_name == "TUMOR"
    assert n_name == "DNA_NORMAL"


def test_parse_genotype_state_and_context_rule():
    primary = parse_genotype_state("0|1", "17005")
    same_hap = parse_genotype_state("0|1", "17005")
    other_hap = parse_genotype_state("1|0", "17005")
    hom_alt = parse_genotype_state("1/1", "")
    assert should_apply_context_variant(primary, same_hap) is True
    assert should_apply_context_variant(primary, other_hap) is False
    assert should_apply_context_variant(primary, hom_alt) is True

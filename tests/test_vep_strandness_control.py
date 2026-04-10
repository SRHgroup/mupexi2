from types import SimpleNamespace

from mupexi2.core import build_vep_info
from mupexi2.vep_strandness_control import (
    choose_keep_strand,
    collect_rna_edit_positions,
    filter_vep_file_by_rna_edit_strandness,
    filter_vep_rows_by_rna_edit_strandness,
    parse_vep_row,
)


def make_vep_row(
    uploaded_variation,
    feature,
    consequence,
    strand,
    amino_acids="I/V",
    protein_position="34",
    feature_type="Transcript",
    extra_suffix="",
    location="1:100",
):
    extra = "SYMBOL=GENE1;STRAND={}".format(strand)
    if extra_suffix:
        extra = "{};{}".format(extra, extra_suffix)
    fields = [
        uploaded_variation,
        location,
        "G",
        "ENSG000001",
        feature,
        feature_type,
        consequence,
        "100",
        "100",
        protein_position,
        amino_acids,
        "Att/Gtt",
        "-",
        extra,
    ]
    return "\t".join(fields) + "\n"


def make_vcf_text(source_set="RNA_EDIT"):
    return (
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tDNA_NORMAL\tTUMOR\n"
        "1\t100\t.\tA\tG\t.\tPASS\tSOURCE_SET={}\tGT:AD:AF:PS\t0/0:12,0:0.0:.\t0|1:10,5:0.33:17005\n".format(source_set)
    )


def test_rna_edit_strandness_control_prefers_non_nmd_strand():
    rows = [
        make_vep_row("var1", "ENST000002", "missense_variant,NMD_transcript_variant", -1, amino_acids="Y/C"),
        make_vep_row("var1", "ENST000001", "missense_variant", 1, amino_acids="I/V"),
    ]
    kept_rows, stats = filter_vep_rows_by_rna_edit_strandness(rows, {("1", "100")})
    assert kept_rows == [make_vep_row("var1", "ENST000001", "missense_variant", 1, amino_acids="I/V")]
    assert stats.rna_edit_variations_with_strand_conflict == 1
    assert stats.vep_rows_dropped_by_strandness_control == 1


def test_rna_edit_strandness_control_prefers_more_frequent_amino_acid_change():
    rows = [
        make_vep_row("var1", "ENST000001", "missense_variant", 1, amino_acids="I/V"),
        make_vep_row("var1", "ENST000002", "missense_variant", 1, amino_acids="I/V"),
        make_vep_row("var1", "ENST000003", "missense_variant", -1, amino_acids="Y/C"),
    ]
    kept_rows, _ = filter_vep_rows_by_rna_edit_strandness(rows, {("1", "100")})
    assert kept_rows == rows[:2]


def test_rna_edit_strandness_control_leaves_single_strand_unchanged():
    rows = [
        make_vep_row("var1", "ENST000001", "missense_variant", 1),
        make_vep_row("var1", "ENST000002", "synonymous_variant", 1, amino_acids="-", protein_position="-"),
    ]
    kept_rows, stats = filter_vep_rows_by_rna_edit_strandness(rows, {("1", "100")})
    assert kept_rows == rows
    assert stats.vep_rows_dropped_by_strandness_control == 0


def test_strandness_control_does_not_touch_non_rna_edit_variants():
    rows = [
        make_vep_row("var1", "ENST000001", "missense_variant", 1),
        make_vep_row("var1", "ENST000002", "missense_variant", -1, amino_acids="Y/C"),
    ]
    kept_rows, stats = filter_vep_rows_by_rna_edit_strandness(rows, set())
    assert kept_rows == rows
    assert stats.rna_edit_uploaded_variations_considered == 0


def test_rna_edit_strandness_control_breaks_ties_deterministically():
    rows = [
        make_vep_row("var1", "ENST000020", "missense_variant", 1, amino_acids="I/V"),
        make_vep_row("var1", "ENST000010", "missense_variant", -1, amino_acids="I/V"),
    ]
    parsed_rows = [parse_vep_row(row) for row in rows]
    assert choose_keep_strand(parsed_rows) == -1
    kept_rows, _ = filter_vep_rows_by_rna_edit_strandness(rows, {("1", "100")})
    assert kept_rows == [rows[1]]


def test_filter_vep_file_by_rna_edit_strandness_and_build_vep_info(tmp_path):
    vep_path = tmp_path / "input.vep"
    vep_path.write_text(
        "# VEP synthetic header\n"
        + make_vep_row("var1", "ENST000002", "missense_variant,NMD_transcript_variant", -1, amino_acids="Y/C")
        + make_vep_row("var1", "ENST000001", "missense_variant", 1, amino_acids="I/V")
    )
    vcf_path = tmp_path / "input.vcf"
    vcf_path.write_text(make_vcf_text())
    filtered_path = tmp_path / "filtered.vep"

    stats = filter_vep_file_by_rna_edit_strandness(str(vep_path), str(filtered_path), str(vcf_path))
    filtered_rows = [line for line in filtered_path.read_text().splitlines() if line and not line.startswith("#")]

    assert len(filtered_rows) == 1
    assert "\tENST000001\t" in filtered_rows[0]
    assert stats.vep_rows_dropped_by_strandness_control == 1

    vep_info, vep_counters, transcript_info, protein_positions = build_vep_info(
        SimpleNamespace(name=str(filtered_path)),
        None,
        SimpleNamespace(name=str(vcf_path)),
        [9],
        tumor_sample="TUMOR",
        normal_sample="DNA_NORMAL",
        vep_strandness_stats=stats,
    )

    assert len(vep_info) == 1
    assert vep_info[0].trans_id == "ENST000001"
    assert vep_counters.rna_edit_variations_with_strand_conflict == 1
    assert vep_counters.vep_rows_dropped_by_strandness_control == 1
    assert list(transcript_info.values())[0] == {"ENSG000001": ["ENST000001"]}
    assert list(protein_positions.values())[0] == {"ENSG000001": {"ENST000001": "34"}}


def test_collect_rna_edit_positions_from_vcf(tmp_path):
    vcf_path = tmp_path / "input.vcf"
    vcf_path.write_text(make_vcf_text())
    assert collect_rna_edit_positions(str(vcf_path)) == {("1", "100")}

from types import SimpleNamespace

from mupexi2.core import build_vep_info
from mupexi2.vep_disambiguation import select_best_vep_rows


def make_vep_row(
    uploaded_variation,
    feature,
    consequence,
    amino_acids="I/V",
    protein_position="34",
    feature_type="Transcript",
    extra="SYMBOL=GENE1;STRAND=1",
    location="1:100",
):
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


def test_select_best_vep_rows_prefers_non_nmd_row():
    kept_rows, stats = select_best_vep_rows([
        make_vep_row("var1", "ENST000002", "missense_variant,NMD_transcript_variant"),
        make_vep_row("var1", "ENST000001", "missense_variant"),
    ])
    assert len(kept_rows) == 1
    assert "\tENST000001\t" in kept_rows[0]
    assert stats.mutations_with_multiple_transcript_rows == 1
    assert stats.transcript_rows_dropped_by_disambiguation == 1


def test_select_best_vep_rows_prefers_majority_amino_acid_change():
    kept_rows, _ = select_best_vep_rows([
        make_vep_row("var1", "ENST000001", "missense_variant", amino_acids="I/V"),
        make_vep_row("var1", "ENST000002", "missense_variant", amino_acids="I/V"),
        make_vep_row("var1", "ENST000003", "missense_variant", amino_acids="Y/C"),
    ])
    assert kept_rows == [make_vep_row("var1", "ENST000001", "missense_variant", amino_acids="I/V")]


def test_select_best_vep_rows_leaves_single_annotation_unchanged():
    row = make_vep_row("var1", "ENST000001", "missense_variant")
    kept_rows, stats = select_best_vep_rows([row])
    assert kept_rows == [row]
    assert stats.transcript_rows_dropped_by_disambiguation == 0


def test_select_best_vep_rows_does_not_filter_by_edit_signature():
    kept_rows, _ = select_best_vep_rows([
        make_vep_row("var1", "ENST000001", "missense_variant", extra="SYMBOL=GENE1;STRAND=-1;EDIT_SIG=NON_ADAR"),
    ])
    assert len(kept_rows) == 1
    assert "EDIT_SIG=NON_ADAR" in kept_rows[0]


def test_select_best_vep_rows_breaks_ties_deterministically():
    kept_rows, _ = select_best_vep_rows([
        make_vep_row("var1", "ENST000020", "missense_variant", amino_acids="I/V"),
        make_vep_row("var1", "ENST000010", "missense_variant", amino_acids="I/V"),
    ])
    assert kept_rows == [make_vep_row("var1", "ENST000010", "missense_variant", amino_acids="I/V")]


def test_build_vep_info_uses_disambiguated_vep_rows(tmp_path):
    vep_path = tmp_path / "input.vep"
    vep_path.write_text(
        "# VEP synthetic header\n"
        + make_vep_row("var1", "ENST000002", "missense_variant,NMD_transcript_variant")
        + make_vep_row("var1", "ENST000001", "missense_variant")
    )

    vcf_path = tmp_path / "input.vcf"
    vcf_path.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tDNA_NORMAL\tTUMOR\n"
        "1\t100\t.\tA\tG\t.\tPASS\tSOURCE_SET=RNA_EDIT\tGT:AD:AF:PS\t0/0:12,0:0.0:.\t0|1:10,5:0.33:17005\n"
    )

    vep_info, vep_counters, transcript_info, protein_positions = build_vep_info(
        SimpleNamespace(name=str(vep_path)),
        None,
        SimpleNamespace(name=str(vcf_path)),
        [9],
        tumor_sample="TUMOR",
        normal_sample="DNA_NORMAL",
    )

    assert len(vep_info) == 1
    assert vep_info[0].trans_id == "ENST000001"
    assert vep_counters.total_vep_rows_read == 2
    assert vep_counters.total_unique_uploaded_variations == 1
    assert vep_counters.mutations_with_multiple_transcript_rows == 1
    assert vep_counters.transcript_rows_dropped_by_disambiguation == 1
    assert list(transcript_info.values())[0] == {"ENSG000001": ["ENST000001"]}
    assert list(protein_positions.values())[0] == {"ENSG000001": {"ENST000001": "34"}}

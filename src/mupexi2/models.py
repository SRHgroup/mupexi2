from dataclasses import dataclass


@dataclass(frozen=True)
class RunConfig:
    vcf_file: str
    fusion_file: str
    peptide_length: str
    output: str
    logfile: str
    HLA_alleles: str
    config: str
    expression_file: str
    fasta_file_name: str
    webserver: str
    outdir: str
    keep_temp: str
    prefix: str
    print_mismatch: str
    liftover: str
    expression_type: str
    num_mismatches: str
    assembly: str
    fork: str
    species: str
    netmhc_anal: str
    vcf_type: str
    germlines: bool
    rna_edit: bool
    tumor_sample: str
    normal_sample: str
    rnaedit_known_only: bool
    rnaedit_known_key: str
    phasing_mode: str
    superpeptides: bool
    parallel_k: bool


@dataclass(frozen=True)
class ContextMutation:
    consequence: str
    protein_pos: str
    aa_change: str
    source: str
    mutation_id: str


@dataclass(frozen=True)
class VariantRecord:
    gene_id: str
    trans_id: str
    mutation_consequence: str
    chr: str
    pos: str
    cdna_pos: str
    prot_pos: int
    prot_pos_to: str
    aa_normal: str
    aa_mut: str
    codon_normal: str
    codon_mut: str
    alt_allele: str
    symbol: str
    variant_type: str
    edit_sig: str
    nearby_germlines: dict


@dataclass(frozen=True)
class PeptideSequenceInfo:
    chop_normal_sequence: str
    mutation_sequence: str
    normal_sequence: str
    mutation_position: object
    consequence: str
    relative_germline_positions: object
    relative_context_sources: object


@dataclass(frozen=True)
class PepMatchInfo:
    normal_peptide: str
    mismatch: int
    mismatch_peptide: str

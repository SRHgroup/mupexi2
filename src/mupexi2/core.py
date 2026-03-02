#!/usr/bin/env python3

"""
MuPeXI - Mutant peptide extractor and Informer
Date: 2026-02-01
By: Grigorii Nos, Enrique Matrinez
MuPeXI.py extracts user defined peptides lengths around missense variant mutations, indels and frameshifts.
Information from each mutation is annotated together with the mutant and normal peptides in the file output.
"""

# Import modules
from __future__ import absolute_import
from __future__ import print_function
from collections import defaultdict, namedtuple
from datetime import datetime
from six.moves import map
from six.moves import range

from Bio.Seq import Seq
from Bio import BiopythonWarning
from six.moves.configparser import ConfigParser
from tempfile import NamedTemporaryFile
import psutil
import sys, re, getopt, itertools, warnings, string, subprocess, os, os.path, math, tempfile, shutil, numpy, pandas
import six
from six.moves import zip
import warnings


def main(args):
    start_time = datetime.now()

    # State input in variables
    input_ = read_options(args)
    version = '2.0.0'

    # Redirect std error when run on webserver
    webserver_err_redirection(input_.webserver)

    # Check input file paths
    species = define_species(input_.species)
    peptide_length = extract_peptide_length(input_.peptide_length)

    paths = check_input_paths(input_, peptide_length, species)
    tmp_dir = create_tmp_dir()
    www_tmp_dir = create_webserver_tmp_dir(input_.webserver)

    # Read in data
    print_ifnot_webserver('\nReading in data', input_.webserver)
    expression = build_expression(input_.expression_file, input_.webserver, input_.expression_type, species)
    proteome_reference, sequence_count = build_proteome_reference(paths.proteome_ref_file, input_.webserver, species)
    genome_reference = build_genome_reference(paths.genome_ref_file, input_.webserver, species)
    cancer_genes = build_cancer_genes(paths.cosmic_file, input_.webserver)

    # Extract reference peptides
    reference_peptides, reference_peptide_counters, \
    reference_peptide_file_names = reference_peptide_extraction(proteome_reference, peptide_length, tmp_dir,
                                                                input_.webserver, input_.keep_temp, input_.prefix,
                                                                input_.outdir, input_.config)

    """
    VEP: Ensembls Variant effect predictor 
    Detecting vcf file input, extracting allele frequencies and running VEP
    """
    start_time_vep = datetime.now()

    if input_.vcf_file is not None:

        print_ifnot_webserver('\nVEP: Starting process for running the Ensembl Variant Effect Predictor',
                              input_.webserver)

        variant_caller = detect_variant_caller(input_.vcf_file, input_.webserver)
        vcf_file = liftover_hg19(input_.liftover, input_.webserver, input_.vcf_file, input_.keep_temp, input_.outdir,
                                 input_.prefix, tmp_dir, input_.config)
        if (input_.germlines or input_.rna_edit) and not vcf_has_sourceset(vcf_file, input_.liftover):
            sys.exit('ERROR: SOURCE_SET is required when --germlines=true or --rna-edit=true')
        vcf_sorted_file = create_vep_compatible_vcf(vcf_file, input_.webserver, input_.keep_temp, input_.outdir,
                                                    input_.prefix, tmp_dir, input_.liftover, input_.germlines,
                                                    input_.rna_edit, input_.rnaedit_known_only,
                                                    input_.rnaedit_known_key)

        snv_qc = extract_snv_qc(vcf_sorted_file, input_.webserver, variant_caller, input_.tumor_sample,
                                input_.normal_sample)
        vep_file = run_vep(vcf_sorted_file, input_.webserver, tmp_dir, paths.vep_path, paths.vep_dir, input_.keep_temp,
                           input_.prefix, input_.outdir, input_.assembly, input_.fork, species)
        vep_info, vep_counters, transcript_info, protein_positions = build_vep_info(vep_file, input_.webserver, vcf_sorted_file, peptide_length)

        end_time_vep = datetime.now()

        """
        MuPeX: Mutant peptide extraction
        Extracting user defined peptides lengths around missense variant mutations, indels and frameshifts
        """
        start_time_mupex = datetime.now()
        print_ifnot_webserver('\nMuPeX: Starting mutant peptide extraction for SNV mutations', input_.webserver)

        # extract mutant peptides
        peptide_info, peptide_counters, fasta_printout, \
        pepmatch_file_names = peptide_extraction(peptide_length, vep_info, proteome_reference, genome_reference,
                                                 reference_peptides, reference_peptide_file_names,
                                                 input_.fasta_file_name,
                                                 paths.peptide_match, tmp_dir, input_.webserver, input_.print_mismatch,
                                                 input_.keep_temp, input_.prefix, input_.outdir, input_.num_mismatches)

        """
        MuPeI: Mutant peptide Informer 
        Information from each mutation is annotated together with the mutant and normal peptides in the file output
        """
        start_time_mupei = datetime.now()
        print_ifnot_webserver('\nMuPeI: Starting mutant peptide informer', input_.webserver)

        # run netMHCpan
        unique_mutant_peptide_count, peptide_file = write_peptide_file(peptide_info, tmp_dir, input_.webserver,
                                                                       input_.keep_temp, input_.prefix, input_.outdir)
        netMHCpan_runtime, unique_alleles, \
        netMHC_EL_file, netMHC_BA_file = run_netMHCpan(input_.HLA_alleles, paths.netMHC, peptide_file, tmp_dir,
                                                       input_.webserver, input_.keep_temp, input_.prefix,
                                                       input_.outdir, input_.netmhc_anal)
        net_mhc_BA = build_netMHC(netMHC_BA_file, input_.webserver, 'YES') if not netMHC_BA_file is None else None
        net_mhc_EL = build_netMHC(netMHC_EL_file, input_.webserver, 'NO')

        # write files
        output_file = write_output_file(peptide_info=peptide_info, expression=expression, net_mhc_BA=net_mhc_BA,
                                        net_mhc_EL=net_mhc_EL, unique_alleles=unique_alleles, cancer_genes=cancer_genes,
                                        tmp_dir=tmp_dir, webserver=input_.webserver,
                                        print_mismatch=input_.print_mismatch,
                                        snv_qc=snv_qc, expression_file_type=input_.expression_type,
                                        transcript_info=transcript_info, reference_peptides=reference_peptides,
                                        proteome_reference=proteome_reference, protein_positions=protein_positions,
                                        version=version, mode='SNV')

        if 'fus_fasta_printout' not in locals(): fus_fasta_printout = {}

        # print fasta file
    if input_.fusion_file is not None:

        print_ifnot_webserver('\nMuPeX: Starting mutant peptide extraction for fusion mutations', input_.webserver)

        fus_info, fus_transcript_info, fus_counter = build_fusion_info(fusion_file=input_.fusion_file,
                                                                       discarded_fusion_file='{}_discarded.tsv'.format(
                                                                           input_.fusion_file.replace('.tsv', '')),
                                                                       webserver=None, junction_filter=3,
                                                                       spanning_filter=1)

        fus_peptide_info, fus_peptide_counters, fus_fasta_printout, fus_pepmatch = fusion_peptide_extraction(
            peptide_lengths=peptide_length, fus_info=fus_info,
            reference_peptide_file_names=reference_peptide_file_names,
            fasta_file_name=f'{input_.fasta_file_name}_fus_peptides',
            tmp_dir=tmp_dir, keep_tmp=True, file_prefix=input_.prefix,
            webserver=input_.webserver, num_mismatches=input_.num_mismatches,
            peptide_match=paths.peptide_match, print_mismatch=input_.print_mismatch,
            outdir=input_.outdir)

        # peptide_info.update(fus_peptide_info)

        unique_fus_mutant_peptide_count, fus_peptide_file = write_peptide_file(peptide_info=fus_peptide_info,
                                                                               tmp_dir=tmp_dir,
                                                                               webserver=input_.webserver,
                                                                               keep_tmp=input_.keep_temp,
                                                                               file_prefix='{}_fusions'.format(
                                                                                   input_.prefix), outdir=input_.outdir)

        fus_netMHCpan_runtime, unique_alleles, \
        fus_netMHC_EL_file, fus_netMHC_BA_file = run_netMHCpan(HLA_alleles=input_.HLA_alleles,
                                                               netMHCpan_path=paths.netMHC,
                                                               peptide_file=fus_peptide_file, tmp_dir=tmp_dir,
                                                               webserver=input_.webserver,
                                                               keep_tmp=input_.keep_temp,
                                                               file_prefix='{}_fusions'.format(input_.prefix),
                                                               outdir=input_.outdir, netmhc_anal=input_.netmhc_anal)

        fus_net_mhc_BA = build_netMHC(fus_netMHC_BA_file, input_.webserver,
                                      'YES') if not fus_netMHC_BA_file is None else None
        fus_net_mhc_EL = build_netMHC(fus_netMHC_EL_file, input_.webserver, 'NO')

        fus_output_file = write_output_file(peptide_info=fus_peptide_info, expression=expression,
                                            net_mhc_BA=fus_net_mhc_BA,
                                            net_mhc_EL=fus_net_mhc_EL, unique_alleles=unique_alleles,
                                            cancer_genes=cancer_genes,
                                            tmp_dir=tmp_dir, webserver=input_.webserver, print_mismatch=None,
                                            snv_qc=None,
                                            expression_file_type=input_.expression_type,
                                            transcript_info=fus_transcript_info,
                                            reference_peptides=None, proteome_reference=None, protein_positions=None,
                                            version=version, mode='FUS')

        if 'fasta_printout' not in locals(): fasta_printout = {}

    fasta_printout.update(fus_fasta_printout)

    fasta_file = write_fasta(tmp_dir, fasta_printout, input_.webserver)

    end_time_mupex = datetime.now()

    try:
        log_file = write_log_file(sys.argv, peptide_length, sequence_count, reference_peptide_counters, vep_counters,
                                  # gotta make this possible for fusion onlys
                                  peptide_counters, start_time_mupex, start_time_mupei, start_time, end_time_mupex,
                                  input_.HLA_alleles, netMHCpan_runtime, unique_mutant_peptide_count,
                                  unique_alleles, tmp_dir, input_.webserver, version)
    except:
        log_file = None

    output_files = {None if input_.vcf_file is None else output_file: '{}_snv.mupexi'.format(input_.prefix),
                    None if input_.fusion_file is None else fus_output_file: '{}_fus.mupexi'.format(input_.prefix)}

    move_output_files(outdir=input_.outdir, log_file=log_file, logfile_name=input_.logfile,
                      fasta_file=fasta_file, fasta_file_name=input_.fasta_file_name, output_files=output_files,
                      webserver=input_.webserver, www_tmp_dir=www_tmp_dir)

    # clean up
    clean_up(tmp_dir)

    webserver_print_output(input_.webserver, www_tmp_dir, input_.output, input_.logfile, input_.fasta_file_name)


"""
CHECK FILE PATHS
"""


def read_dataset(file, sep='\t', header=0, index_col=0, comment=None):
    """
    See pandas.read_csv documentation - this is just a wrapper
    :param file:
    :param sep:
    :param header:
    :param index_col:
    :param comment:
    :return:
    """
    return pandas.read_csv(file, sep=sep, header=header, index_col=index_col,
                           na_values=['Na', 'NA', 'NAN', '.'], comment=comment)


def check_input_paths(input_, peptide_lengths, species):
    file_names = defaultdict(dict)  # empty dictionary
    category_list, id_list = create_lits_for_path_checkup(peptide_lengths)

    # parse files stated in config file
    for category, ID in zip(category_list, id_list):
        file_path = config_parse(input_.config, category, ID)
        check_path(file_path)
        file_names[ID] = file_path

    # Exit program if genome or proteome reference is not stated
    if file_names['cDNA'] is None:
        usage()
        sys.exit('cDNA reference file is missing!')
    if file_names['pep'] is None:
        usage()
        sys.exit('Proteome reference file is missing!')

    if input_.expression_file is not None:
        check_path(input_.expression_file)
        check_file_size(input_.webserver, input_.expression_file, 'expression file')

    if input_.vcf_file is not None:  # gotta create check_fusion_file too
        check_vcf_file(input_.vcf_file, input_.liftover, species, input_.webserver)

    check_netMHC_path(file_names['MHC'])
    check_path(input_.outdir)

    # Create and fill named tuple
    Paths = namedtuple('files', ['gene_symbol_file', 'cosmic_file', 'genome_ref_file', 'proteome_ref_file',
                                 'netMHC', 'peptide_match', 'vep_path', 'vep_dir'])
    paths = Paths(file_names['symbol'], file_names['cosmic'], file_names['cDNA'], file_names['pep'],
                  file_names['MHC'], file_names['PM'], file_names['VEP'], file_names['VEPdir'])

    return paths


def create_lits_for_path_checkup(peptide_lengths):
    for length in peptide_lengths:
        category_list = ['netMHC', 'PeptideMatch', 'EnsemblVEP', 'EnsemblVEP', 'References', 'References', 'References']
        id_list = ['MHC', 'PM', 'VEP', 'VEPdir', 'cosmic', 'cDNA', 'pep']
        category_list.append('References')
        id_list.append('pep' + str(length))
    return category_list, id_list


def config_parse(config_file, category, ID):
    # check if config.ini file exists
    if not os.path.exists(config_file):
        usage()
        sys.exit(
            'ERROR: Path to OR config.ini does not exist!\nERROR: Use -c option to specify path/to/config.ini file\n')

    # parse config file
    config = ConfigParser()
    config.read(config_file)
    path = os.path.expandvars(config.get(category, ID)) if config.has_option(category, ID) else None

    return path


def check_path(path):
    if path is not None:
        if not os.path.exists(path):
            usage()
            sys.exit('ERROR: {} path or file does not exist\n'.format(path))


def check_vcf_file(vcf_file, liftover, species, webserver):
    check_path(vcf_file)

    # Exit program if file size are exceeded
    check_file_size(webserver, vcf_file, 'VCF file')

    with open(vcf_file) as f:
        first_line = f.readline()
        if not first_line.startswith('##fileformat=VCF'):
            usage()
            sys.exit('ERROR: {} file is not a VCF file\n'.format(vcf_file))
        for line in f.readlines():
            if not webserver is None:
                if species == 'human':
                    if '##reference' in line:
                        if 'GRCh37' in line or 'hg19' in line or 'HG19' in line:
                            if liftover is None and species.assembly == "GRCh38":
                                usage()
                                sys.exit(
                                    'ERROR: The VCF file is aligned to HG19 / GRCh37\nINFO: {}\nINFO: Please run NGS '
                                    'analysis aligning to GRCh38, use the hg19 liftover option (-g/--hg19), or the '
                                    'GRCh37 prediction option (-a/ --assembly GRCh37)\n'.format(line.strip()))
                            else:
                                continue
            if '#CHROM' in line:
                break
            elif not line.startswith('##'):
                usage()
                sys.exit('ERROR: {} file is not a VCF file, #CHROM header is missing\n'.format(vcf_file))


def check_netMHC_path(netMHC_path):
    if 'netH2pan' in netMHC_path:
        print('\tMouse specific MHC binding predictor netH2pan used')
    elif 'netMHCpan' in netMHC_path:
        if not '4.0' in netMHC_path:
            usage()
            sys.exit(
                'ERROR:\tnetMHCpan version 4.0 not stated in path {}\n\tOnly this version is supported'.format(
                    netMHC_path))
    else:
        usage()
        sys.exit(
            'ERROR:\tnetMHCpan / netH2pan not stated in path {} \n\t\tCheck the correct path to the netXpan binding '
            'predictor is annotated in the config.ini'.format(netMHC_path))


def check_file_size(webserver, input_file_path, file_tag):
    if webserver is not None:
        if os.path.getsize(input_file_path) > 20000000:
            sys.exit('STOP: The input {} size exceeds the limit of 20M\n\tPlease check your file\n\t'
                     'If correct use the command-line tool instead or contact us: '
                     'ambj@cbs.dtu.dk or eklund@cbs.dtu.dk'.format(file_tag))


def create_tmp_dir():
    tmp_dir = tempfile.mkdtemp(prefix='MuPeXI_')
    os.system('chmod -R a+rwx {}'.format(tmp_dir))
    return tmp_dir


def print_ifnot_webserver(string, webserver):
    if webserver is None:
        print(string)


def create_webserver_tmp_dir(webserver):
    if webserver is not None:
        www_tmp_dir = tempfile.mkdtemp(prefix='MuPeXI_', dir='/usr/opt/www/pub/CBS/services/MuPeXI-1.1/tmp')
    else:
        www_tmp_dir = None
    return www_tmp_dir


def webserver_err_redirection(webserver):
    if webserver is not None:
        sys.stderr = sys.stdout


def print_mem_usage():
    # total memory usage
    process = psutil.Process(os.getpid())
    process_mb = float(process.memory_info().rss / 1000000)
    print('Total mem usage: {} MB'.format(process_mb))


"""
READ IN DATA
"""


def build_proteome_reference(proteome_ref_file, webserver, species):
    print_ifnot_webserver('\tCreating proteome reference dictionary', webserver)
    proteome_reference = defaultdict(dict)  # empty dictionary
    sequence_count = 0  # sequence count

    with open(proteome_ref_file) as f:
        for line in f.readlines():
            if line.startswith('>'):  # fasta header (>)
                # test species compatibility
                if not species.gene_id_prefix in line:
                    usage()
                    sys.exit(
                        'ERROR: species prefix {} for {} not found in proteome reference file\nINFO: '
                        'use --species option if another species than {} is used'.format(
                            species.gene_id_prefix, species.species, species.species))
                # Save gene and transcript Ensembl ID (re = regular expression)
                geneID = re.search(r'gene:({}\d+)'.format(species.gene_id_prefix), line).group(1).strip()
                transID = re.search(r'transcript:({}\d+)'.format(species.trans_id_prefix), line).group(1).strip()
                # Insert gene and transcript ID in directory, assign empty value
                proteome_reference[geneID][transID] = ""
                sequence_count += 1
            else:
                # Assign amino acid sequence as value in directory (the following lines until next header)
                proteome_reference[geneID][transID] += line.strip()
    return proteome_reference, sequence_count


def build_genome_reference(genome_ref_file, webserver, species):
    print_ifnot_webserver('\tCreating genome reference dictionary', webserver)
    genome_reference = defaultdict(dict)  # empty dictionary

    with open(genome_ref_file) as f:
        for line in f.readlines():
            if line.startswith('>'):  # fasta header (>)
                # test species compatibility
                if not species.gene_id_prefix in line:
                    usage()
                    sys.exit(
                        'ERROR: species prefix {} for {} not found in cDNA reference file\nINFO: '
                        'use --species option if another species than {} is used'.format(
                            species.gene_id_prefix, species.species, species.species))
                # Save gene and transcript Ensembl ID (re = regular expression)
                geneID = re.search(r'gene:({}\d+)'.format(species.gene_id_prefix), line).group(1).strip()
                transID = re.search(r'>({}\d+)'.format(species.trans_id_prefix), line).group(1).strip()
                # Insert gene and transcript ID in directory, assign empty value
                genome_reference[geneID][transID] = ""
            else:
                # Assign amino acid sequence as value in directory (the following lines until next header)
                genome_reference[geneID][transID] += line.strip()
    return genome_reference


def build_expression(expression_file, webserver, expression_type, species):
    if expression_file is not None:
        print_ifnot_webserver('\tCreating expression file dictionary', webserver)
        expression = defaultdict(dict)  # empty dictionary

        with open(expression_file) as f:
            for line in f.readlines():
                line = line.split()
                if not line[0] == 'target_id':
                    check_expression_file_type(expression_type, line, species)
                    # save line information
                    if '.' in line[0]:  # example: ENST00000415118.3
                        ensembl_id = line[0].split('.')[0]
                    else:  # example: ENST00000415118
                        ensembl_id = line[0]
                    mean = line[1]
                    # fill dictionary
                    expression[ensembl_id] = mean
    else:
        expression = None

    return expression


def define_species(species):
    if species == 'human':
        trans_id_prefix = 'ENST'
        gene_id_prefix = 'ENSG'
        assembly = 'GRCh38'
        species = 'homo_sapiens'
    elif species == 'mouse':
        trans_id_prefix = 'ENSMUST'
        gene_id_prefix = 'ENSMUSG'
        assembly = 'GRCm38'
        species = 'mus_musculus'
    elif species == 'mouse_balbc':
        trans_id_prefix = 'MGP_BALBcJ_T'
        gene_id_prefix = 'MGP_BALBcJ_G'
        assembly = 'BALB_cJ_v1'
        species = 'mus_musculus_balbcj'
    elif species == 'mouse_black6':
        trans_id_prefix = 'MGP_C57BL6NJ_T'
        gene_id_prefix = 'MGP_C57BL6NJ_G'
        assembly = 'C57BL_6NJ_v1'
        species = 'mus_musculus_c57bl6nj'
    else:
        usage()
        sys.exit('ERROR:\tSpecies {} not recognized \n'.format(species))

    Species = namedtuple('species', ['gene_id_prefix', 'trans_id_prefix', 'assembly', 'species'])
    species = Species(gene_id_prefix, trans_id_prefix, assembly, species)

    return species


def check_expression_file_type(expression_type, line, species):
    if expression_type == 'gene':
        if species.trans_id_prefix in line[0]:
            usage()
            sys.exit(
                'ERROR:\tEnsembl transcript id detected: {}\n\tWhile expression type option "gene" was used\n'.format(
                    line[0]))
    if expression_type == 'transcript':
        if species.gene_id_prefix in line[0]:
            usage()
            sys.exit(
                'ERROR:\tEnsembl gene id detected: {}\n\tWhile expression type option "transcript" was used\n'.format(
                    line[0]))


def build_cancer_genes(cosmic_file, webserver):
    if cosmic_file is not None:
        print_ifnot_webserver('\tCreating cancer genes list', webserver)
        cancer_genes = set()

        with open(cosmic_file) as f:
            for line in f.readlines():
                if not line.startswith('Gene Symbol'):
                    gene_symbol = line.split()[0].strip()
                    cancer_genes.add(gene_symbol)
    else:
        cancer_genes = None

    return cancer_genes


def extract_peptide_length(peptide_length):
    peptide_length_list = list()
    if isinstance(peptide_length, int):
        peptide_length_list.append(peptide_length)
    else:
        if any(i.isdigit() for i in peptide_length):
            if '-' in peptide_length:
                length_range = list(map(int, peptide_length.split('-')))
                assert length_range[1] > length_range[
                    0], 'peptide length range should be stated from lowest to highest. your input: {}'.format(
                    peptide_length)
                peptide_length_list = list(range(length_range[0], length_range[1] + 1))
            elif ',' in peptide_length:
                peptide_length_list = list(map(int, peptide_length.split(',')))
            else:
                peptide_length_list.append(int(peptide_length))
        else:
            sys.exit('ERROR: peptide_length argument must contain integers')

    if sum(n < 0 for n in peptide_length_list) > 0:
        sys.exit('ERROR: peptide_length argument must not have negative numbers')

    return peptide_length_list


"""
VEP
"""


def detect_variant_caller(vcf_file, webserver):
    print_ifnot_webserver('\tDetecting variant caller', webserver)
    with open(vcf_file) as f:
        for line in f.readlines():
            if any(ids in line for ids in ['ID=MuTect2,', 'ID=Mutect2,']):
                variant_caller = 'MuTect2'
                print_ifnot_webserver('\t\tMuTect2', webserver)
                break
            if 'ID=MuTect,' in line:
                variant_caller = 'MuTect'
                print_ifnot_webserver('\t\tMuTect', webserver)
                break
            elif not line.startswith('##'):
                variant_caller = None
                print_ifnot_webserver(
                    '\t\tVariant caller not detected in VCF file.\n\t\tNOTE:\tGenomic allele frequency is only '
                    'taken into account\n\t\t\twith variant calls from MuTect or MuTect2!',
                    webserver)
                break
    return variant_caller


def parse_info_field(info_field):
    info = {}
    for token in info_field.strip().split(';'):
        if not token:
            continue
        if '=' in token:
            key, value = token.split('=', 1)
            info[key] = value
        else:
            info[token] = True
    return info


def infer_source_set_from_info(info_field):
    info = parse_info_field(info_field)
    source_set = info.get('SOURCE_SET', None)
    if source_set is not None:
        source_set = source_set.strip().upper()
        if source_set in ('SOMATIC', 'GERMLINE', 'RNA_EDIT'):
            return source_set
    if 'SOMATIC' in info:
        return 'SOMATIC'
    if 'GERMLINE' in info:
        return 'GERMLINE'
    return 'UNKNOWN'


def vcf_has_sourceset(vcf_file, liftover):
    vcf_file_name = vcf_file if liftover is None else vcf_file.name
    with open(vcf_file_name) as f:
        for line in f:
            if line.startswith('#'):
                continue
            columns = line.strip().split('\t')
            if len(columns) > 7 and 'SOURCE_SET=' in columns[7]:
                return True
    return False


def liftover_hg19(liftover, webserver, vcf_file, keep_tmp, outdir, file_prefix, tmp_dir, config_file):
    if liftover is not None:
        print_ifnot_webserver('\tPerforming Liftover', webserver)
        chain_file = config_parse(config_file, 'LiftOver', 'chain')
        fasta_file = config_parse(config_file, 'LiftOver', 'fasta')
        picard_path = config_parse(config_file, 'LiftOver', 'picard')
        java8 = config_parse(config_file, 'LiftOver', 'java8')
        check_path(chain_file)
        check_path(fasta_file)

        vcf_liftover_file = NamedTemporaryFile(delete=False, dir=tmp_dir, suffix='.vcf')
        rejected_records_file = NamedTemporaryFile(delete=False, dir=tmp_dir, suffix='.vcf')

        # Run picard tools lift-over - only on web-server.
        p1 = subprocess.Popen([java8, '-jar', '{}/picard.jar'.format(picard_path), 'LiftoverVcf',
                               'I={}'.format(vcf_file),
                               'O={}'.format(vcf_liftover_file.name),
                               'CHAIN={}'.format(chain_file),
                               'REJECT={}'.format(rejected_records_file.name),
                               'REFERENCE_SEQUENCE={}'.format(fasta_file)],
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE)
        output, error = p1.communicate()
        vcf_liftover_file.close()

        # Test if VEP file is empty
        if os.stat(vcf_liftover_file.name).st_size <= 23000:
            sys.exit('ERROR: Liftover non fuctioning \nLiftOver {}'.format(error))

        keep_temp_file(keep_tmp, 'vcf', vcf_liftover_file.name, file_prefix, outdir, None, 'vcf_liftover')
        keep_temp_file(keep_tmp, 'vcf', rejected_records_file.name, file_prefix, outdir, None,
                       'rejected_records_liftover')

        vcf_final_file = vcf_liftover_file
    else:
        vcf_final_file = vcf_file

    return vcf_final_file


def create_vep_compatible_vcf(vcf_file, webserver, keep_tmp, outdir, file_prefix, tmp_dir, liftover, germlines=False,
                              rna_edit=False, rnaedit_known_only=True, rnaedit_known_key='KNOWN_RNAEDIT_DB'):
    print_ifnot_webserver('\tChange VCF to the VEP compatible', webserver)
    vcf_file_name = vcf_file if liftover is None else vcf_file.name
    vcf_sorted_file = NamedTemporaryFile(delete=False, dir=tmp_dir)
    with open(vcf_file_name) as source, open(vcf_sorted_file.name, 'w') as out:
        for line in source:
            if line.startswith('#'):
                out.write(line)
                continue

            columns = line.rstrip('\n').split('\t')
            if len(columns) < 8:
                continue

            chrom = re.sub(r'^chr', '', columns[0])
            if re.match(r'^0[0-9]+$', chrom):
                chrom = chrom.lstrip('0')
            columns[0] = chrom if chrom else '0'

            source_set = infer_source_set_from_info(columns[7])
            filter_value = columns[6].strip()
            info_fields = parse_info_field(columns[7])
            keep_record = False

            if source_set == 'SOMATIC':
                keep_record = filter_value == 'PASS'
            elif source_set == 'GERMLINE':
                keep_record = germlines and filter_value == 'PASS'
            elif source_set == 'RNA_EDIT':
                if rna_edit:
                    keep_record = (rnaedit_known_key in info_fields) if rnaedit_known_only else True
            else:
                # Legacy mutect2 style: keep PASS.
                keep_record = filter_value == 'PASS'

            if keep_record:
                out.write('\t'.join(columns) + '\n')
    vcf_sorted_file.close()

    keep_temp_file(keep_tmp, 'vcf', vcf_sorted_file.name, file_prefix, outdir, None, 'vep_compatible_vcf')

    return vcf_sorted_file


def resolve_sample_indices(header, tumor_sample=None, normal_sample=None):
    sample_names = header[9:]
    sample_index = {name: i + 9 for i, name in enumerate(sample_names)}

    def find_unique(candidates):
        hits = [name for name in sample_names if candidates(name)]
        return hits[0] if len(hits) == 1 else None

    if tumor_sample is not None and tumor_sample not in sample_index:
        sys.exit('ERROR: tumor sample "{}" is not present in VCF header'.format(tumor_sample))
    if normal_sample is not None and normal_sample not in sample_index:
        sys.exit('ERROR: normal sample "{}" is not present in VCF header'.format(normal_sample))

    tumor_name = tumor_sample
    normal_name = normal_sample

    if tumor_name is None:
        tumor_name = find_unique(lambda n: n == 'TUMOR') or \
                     find_unique(lambda n: n.endswith('_T') or n.endswith('-T') or n.endswith('.T')) or \
                     find_unique(lambda n: 'TUMOR' in n.upper())
    if normal_name is None:
        normal_name = find_unique(lambda n: n in ('NORMAL', 'DNA_NORMAL')) or \
                      find_unique(lambda n: n.endswith('_N') or n.endswith('-N') or n.endswith('.N')) or \
                      find_unique(lambda n: 'NORMAL' in n.upper())

    if tumor_name is None or normal_name is None or tumor_name == normal_name:
        if len(sample_names) < 2:
            sys.exit('ERROR: VCF must contain at least two sample columns for SNV QC extraction')
        tumor_name = sample_names[0]
        normal_name = sample_names[1]

    return sample_index[tumor_name], sample_index[normal_name], tumor_name, normal_name


def extract_tumor_vaf(vcf_sorted_file, webserver, variant_caller, tumor_prefix='_T'):
    """
    :param vcf_sorted_file: vcf file as the output of create_vep_compatible_vcf() function
    :param webserver: bool, whether you perform analysis on the webserver
    :param variant_caller: str, what program was used for the variant colling, Mutec and Mutec2 preffered
    :param tumor_prefix: str, prefix in the VCF file that identifies column with tumor mutational metrics
    :return: defaultdict, dict with mutation_id : tumor_vaf values
    """

    print_ifnot_webserver('\tExtracting allele frequencies', webserver)
    allele_fractions = defaultdict(dict)

    if variant_caller not in ('MuTect', 'MuTect2'):
        allele_fractions = None  # no variant caller detected
    else:
        with open(vcf_sorted_file.name) as f:
            for line in f.readlines():
                if line.startswith('##'):
                    continue
                elif line.startswith('#'):
                    header = line.split('\t')
                    header = [re.sub(r'\n|#', '', file) for file in header]

                    chr_i = header.index('CHROM')
                    pos_i = header.index('POS')
                    ref_i = header.index('REF')
                    alt_i = header.index('ALT')
                    format_i = header.index('FORMAT')

                    tumor_meta_i = [i for i, item in enumerate(header) if item.endswith(tumor_prefix)][0]
                else:
                    columns = line.split('\t')
                    chromosome = columns[chr_i].strip()
                    genomic_position = columns[pos_i].strip()
                    reference_allele = columns[ref_i].strip()
                    altered_allele = columns[alt_i].strip()
                    format_fields = columns[format_i].strip().split(':')
                    if not len(reference_allele) == len(altered_allele):
                        altered_allele = altered_allele[1:] if len(reference_allele) < len(
                            altered_allele) else altered_allele
                        altered_allele = '-' if len(reference_allele) > len(altered_allele) else altered_allele
                    if variant_caller == 'MuTect':
                        allele_fraction = columns[tumor_meta_i].split(':')[
                            format_fields.index('FA')].strip()  # GT:AD:BQ:DP:FA
                    elif variant_caller == 'MuTect2':
                        allele_fraction = columns[tumor_meta_i].split(':')[
                            format_fields.index('AF')].strip()  # GT:AD:AF:ALT_F1R2:ALT_F2R1:FOXOG:QSS:REF_F1R2:REF_F2R1
                    ID = '{chromosome}_{genomic_position}_{altered_allele}'.format(
                        chromosome=chromosome,
                        genomic_position=genomic_position if not altered_allele == '-' else int(genomic_position) + 1,
                        altered_allele=altered_allele)
                    allele_fractions[ID] = allele_fraction

    return allele_fractions


def extract_snv_qc(vcf_sorted_file, webserver, variant_caller, tumor_sample=None, normal_sample=None):
    """
    :param vcf_sorted_file: vcf file as the output of create_vep_compatible_vcf() function
    :param webserver: bool, whether you perform analysis on the webserver
    :param variant_caller: str, what program was used for the variant colling, MuTect and MuTect2 preffered

    :return: defaultdict, dict of dicts
    {mutation_id : {tumor_vaf, normal_vaf, tumor_depth, normal_depth, alt_read_counts}}
    """

    print_ifnot_webserver('\tExtracting SNV QC metrics', webserver)
    snv_qc = defaultdict(dict)
    with open(vcf_sorted_file.name) as f:
        for line in f.readlines():
            if line.startswith('##'):
                continue
            elif line.startswith('#'):
                header = line.split('\t')
                header = [re.sub(r'\n|#', '', file) for file in header]

                chr_i = header.index('CHROM')
                pos_i = header.index('POS')
                ref_i = header.index('REF')
                alt_i = header.index('ALT')
                info_i = header.index('INFO')
                format_i = header.index('FORMAT')

                tumor_meta_i, normal_meta_i, detected_tumor, detected_normal = resolve_sample_indices(
                    header, tumor_sample, normal_sample)
                if tumor_sample is None or normal_sample is None:
                    print_ifnot_webserver('\t\tUsing samples tumor="{}" normal="{}" for QC'.format(
                        detected_tumor, detected_normal), webserver)
            else:
                columns = line.split('\t')
                chromosome = columns[chr_i].strip()
                genomic_position = columns[pos_i].strip()
                reference_allele = columns[ref_i].strip()
                altered_allele = columns[alt_i].strip()
                source_set = infer_source_set_from_info(columns[info_i].strip())
                if source_set == 'GERMLINE':
                    continue

                if not len(reference_allele) == len(altered_allele):
                    altered_allele = altered_allele[1:] if len(reference_allele) < len(altered_allele) else altered_allele
                    altered_allele = '-' if len(reference_allele) > len(altered_allele) else altered_allele

                format_fields = columns[format_i].strip().split(':')
                tumor_values = columns[tumor_meta_i].strip().split(':')
                normal_values = columns[normal_meta_i].strip().split(':')
                tumor_map = dict(zip(format_fields, tumor_values))
                normal_map = dict(zip(format_fields, normal_values))

                preferred_vaf_key = 'FA' if variant_caller == 'MuTect' else 'AF'
                vaf_column = preferred_vaf_key if preferred_vaf_key in format_fields else (
                    'AF' if 'AF' in format_fields else ('FA' if 'FA' in format_fields else None))

                def parse_ad(ad_value):
                    if ad_value is None or ad_value in ('.', './.'):
                        return None
                    parts = ad_value.split(',')
                    if len(parts) < 2:
                        return None
                    try:
                        ref = int(parts[0])
                        alt = int(parts[1])
                    except ValueError:
                        return None
                    return ref, alt

                tumor_ad = parse_ad(tumor_map.get('AD'))
                normal_ad = parse_ad(normal_map.get('AD'))

                if vaf_column is not None:
                    tumor_vaf = tumor_map.get(vaf_column, 'NA')
                    normal_vaf = normal_map.get(vaf_column, 'NA')
                else:
                    tumor_vaf = 'NA'
                    normal_vaf = 'NA'

                if tumor_vaf in ('NA', '.', './.') and tumor_ad is not None:
                    t_depth_calc = tumor_ad[0] + tumor_ad[1]
                    tumor_vaf = '{:.6g}'.format(float(tumor_ad[1]) / t_depth_calc) if t_depth_calc > 0 else 'NA'
                if normal_vaf in ('NA', '.', './.') and normal_ad is not None:
                    n_depth_calc = normal_ad[0] + normal_ad[1]
                    normal_vaf = '{:.6g}'.format(float(normal_ad[1]) / n_depth_calc) if n_depth_calc > 0 else 'NA'

                t_alt_count = tumor_ad[1] if tumor_ad is not None else 'NA'
                t_depth = (tumor_ad[0] + tumor_ad[1]) if tumor_ad is not None else 'NA'
                n_depth = (normal_ad[0] + normal_ad[1]) if normal_ad is not None else 'NA'

                ID = '{chromosome}_{genomic_position}_{altered_allele}'.format(
                    chromosome=chromosome,
                    genomic_position=genomic_position if not altered_allele == '-' else int(genomic_position) + 1,
                    altered_allele=altered_allele)
                snv_qc[ID] = {'tumor_vaf': tumor_vaf, 'normal_vaf': normal_vaf,
                              't_depth': t_depth, 'n_depth': n_depth, 't_alt_count': t_alt_count}

    return snv_qc


def run_vep(vcf_sorted_file, webserver, tmp_dir, vep_path, vep_dir, keep_tmp, file_prefix, outdir, input_assembly, fork,
            species):
    print_ifnot_webserver('\tRunning VEP', webserver)
    vep_file = NamedTemporaryFile(delete=False, dir=tmp_dir)

    assembly = species.assembly if input_assembly is None else input_assembly

    popen_args = [
        vep_path,
        '-fork', str(fork),
        '--offline',
        '--quiet',
        '--assembly', assembly,
        '--species', species.species,
        '--dir', vep_dir,
        '-i', vcf_sorted_file.name,
        '--force_overwrite',
        '--symbol',  # print gene symbol in output
        '-o', vep_file.name]

    p1 = subprocess.Popen(popen_args,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE)
    output, error = p1.communicate()
    vep_file.close()

    # Test if VEP file is empty
    if os.stat(vep_file.name).st_size == 0:
        print(('\nERROR:\tVEP output file empty\nVEP:\t{}'.format(error)))
        if "use an undefined value as a symbol reference at" in error:
            print("NOTE:\tCheck the VCF file have been generated with The GRCh38 assembly - or use the liftover option")
        sys.exit()

    keep_temp_file(keep_tmp, 'vep', vep_file.name, file_prefix, outdir, None, 'vep')

    return vep_file


def build_vep_info(vep_file, webserver, vep_compatible_vcf, peptide_length):
    print_ifnot_webserver('\tCreating mutation information dictionary', webserver)
    vep_info = []  # empty list
    previous_mutation_id = previous_mutation_id_vep = ''  # empty string
    # Creating named tuple
    Mutation_Info = namedtuple('mutation_info',
                               ['gene_id', 'trans_id', 'mutation_consequence', 'chr', 'pos', 'cdna_pos', 'prot_pos',
                                'prot_pos_to', 'aa_normal', 'aa_mut', 'codon_normal', 'codon_mut', 'alt_allele',
                                'symbol','variant_type','nearby_germlines'])
    transcript_info = defaultdict(dict)
    protein_positions = defaultdict(lambda: defaultdict(dict))
    germline_info = {}

    non_used_mutation_count, misssense_variant_count, inframe_insertion_count, \
    inframe_deletion_count, frameshift_variant_count = 0, 0, 0, 0, 0

    with open(vep_file.name) as f, open(vep_compatible_vcf.name) as vcf:

        variant_types = {}
        germline_positions = {}

        for line in vcf.readlines():
            if line.startswith('#'):
                continue
            columns = line.strip().split("\t")
            chrom_pos = (columns[0], columns[1])
            variant_types[chrom_pos] = infer_source_set_from_info(columns[7])

        for line in f.readlines():
            if line.startswith('#'):
                continue
            consequence = line.strip().split('\t')[6]

            if consequence == 'missense_variant' or consequence == 'inframe_insertion':
                chrom_pos = (line.strip().split('\t')[1].split(':')[0], line.strip().split('\t')[1].split(':')[1])
                variant_type = variant_types.get((chrom_pos[0], chrom_pos[1]))
            elif consequence == 'inframe_deletion':
                chrom_pos = (line.strip().split('\t')[1].split(':')[0], line.strip().split('\t')[1].split(':')[1])
                vep_to_vcf_pos = str(int(chrom_pos[1].split('-')[0])-1)
                variant_type = variant_types.get((chrom_pos[0], vep_to_vcf_pos))
            else:
                continue

            if variant_type == 'GERMLINE':
                variant_info = [consequence, line.strip().split('\t')[9], line.strip().split('\t')[10]]
                germline_positions[chrom_pos] = variant_info

    with open(vep_file.name) as f:

        for i, line in enumerate(f.readlines()):
            if line.startswith('#'):  # skip lines starting with #
                continue
            # Go to the next line if it is not a missense variant mutations
            mutation_consequence = ['missense_variant', 'inframe_insertion', 'inframe_deletion', 'frameshift_variant']
            mutation_id = line.split('\t')[0].strip()
            if not any(wanted_consequence in line for wanted_consequence in mutation_consequence):
                if not mutation_id == previous_mutation_id:
                    non_used_mutation_count += 1
                    # save previous mutation ID
                    previous_mutation_id = mutation_id
                continue
            if 'stop' in line:
                continue
            line = line.split('\t')
            # save relevant information from line
            mutation_consequence = line[6].split(',')[0].strip()
            chr_, genome_pos = line[1].split(':')
            alt_allele = line[2].strip()
            geneID = line[3].strip()
            transID = line[4].strip()
            prot_pos = line[9].strip()
            cdna_pos = line[7].strip()
            symbol = re.search(r'SYMBOL=(\d*\w*\d*\w*\d*\w*)',
                               line[13]).group(1).strip() if not re.search(r'SYMBOL', line[13]) is None else '-'
            aa_normal, aa_mutation = line[10].split('/')
            codon_normal, codon_mut = line[11].split('/')
            if '-' in line[9]:
                prot_pos, prot_pos_to = line[9].split('-')
            else:
                prot_pos, prot_pos_to = line[9].strip(), None
            variant_type = variant_types.get((chr_, genome_pos), 'UNKNOWN')
            if variant_type == 'SOMATIC':
                nearby_germlines = get_nearby_germlines(chr_, genome_pos, germline_positions, peptide_length)
            elif variant_type == 'GERMLINE':
                nearby_germlines = None
            else:
                nearby_germlines = None
            mutation_id_vep = '{}_{}_{}/{}'.format(chr_, genome_pos, aa_normal, aa_mutation)
            # Generate dict of dicts (dependent on both mutation ID and gene ID)
            # then set the default value of the key to be a list and append the transcript id
            transcript_info[mutation_id_vep].setdefault(geneID, []).append(transID)
            # ad protein position
            protein_positions[mutation_id_vep][geneID][transID] = prot_pos
            # append information from the line to the list of named tuples - fill tuple
            vep_info.append(
                Mutation_Info(geneID, transID, mutation_consequence, chr_, genome_pos, cdna_pos, int(prot_pos),
                              prot_pos_to, aa_normal, aa_mutation, codon_normal, codon_mut, alt_allele, symbol, variant_type, nearby_germlines))

            # count independent mutation mutation consequences
            if (not mutation_id_vep == previous_mutation_id_vep) and mutation_consequence == 'missense_variant':
                misssense_variant_count += 1
            if (not mutation_id_vep == previous_mutation_id_vep) and mutation_consequence == 'inframe_insertion':
                inframe_insertion_count += 1
            if (not mutation_id_vep == previous_mutation_id_vep) and mutation_consequence == 'inframe_deletion':
                inframe_deletion_count += 1
            if (not mutation_id_vep == previous_mutation_id_vep) and mutation_consequence == 'frameshift_variant':
                frameshift_variant_count += 1

            # save previous mutation ID
            previous_mutation_id_vep = mutation_id_vep

    if not vep_info:
        sys.exit(
            '\nNO RELEVANT MUTATIONS (missense variant, inframe insertion / deletion or frameshift variant) '
            'found in VEP file.\nMuPeXI run stopped\nPrint temporary files to check if this is correct (option -t)\n')

    # Count unique gene and transcript IDs
    gene_IDs, transcript_IDs = set(), set()
    for mutation in vep_info:
        gene_IDs.add(mutation.gene_id)
        transcript_IDs.add(mutation.trans_id)
    gene_count, transcript_Count = len(gene_IDs), len(transcript_IDs)

    # Create and fill counter named tuple
    VEPCounters = namedtuple('vep_counters',
                             ['non_used_mutation_count', 'misssense_variant_count', 'gene_count', 'transcript_Count',
                              'inframe_insertion_count', 'inframe_deletion_count', 'frameshift_variant_count'])
    vep_counters = VEPCounters(non_used_mutation_count, misssense_variant_count, gene_count, transcript_Count,
                               inframe_insertion_count, inframe_deletion_count, frameshift_variant_count)

    return vep_info, vep_counters, transcript_info, protein_positions

#### function to get germlines surrounding a somatic with a given peptide length. returns a dict with position:AAchange

def get_nearby_germlines(chrom, pos, germline_positions, peptide_length):
    nearby_germlines = {}
    max_distance = peptide_length[-1]*3
    for (g_chrom, g_pos), info in germline_positions.items():
        if germline_positions[(g_chrom, g_pos)][0] == 'missense_variant':
            if '-' in g_pos or '-' in pos:
                continue
            if g_chrom == chrom and abs(int(g_pos) - int(pos)) <= max_distance:
                nearby_germlines[(g_chrom, g_pos)] = info
        elif germline_positions[(g_chrom, g_pos)][0] == 'inframe_insertion':
            start, stop = g_pos.split('-')
            if '-' in g_pos or '-' in pos:
                continue
            if g_chrom == chrom and (abs(int(start)-int(pos)) <= max_distance or abs(int(stop)-int(pos)) <= max_distance):
                nearby_germlines[(g_chrom, g_pos)] = info
        elif germline_positions[(g_chrom, g_pos)][0] == 'inframe_deletion':
            start, stop = g_pos.split('-')
            if '-' in g_pos or '-' in pos:
                continue
            if g_chrom == chrom and abs(int(start)-int(pos)) <= max_distance:
                nearby_germlines[(g_chrom, g_pos)] = info
    nearby_germlines = dict(sorted(nearby_germlines.items()))
    if not nearby_germlines:
        nearby_germlines = None
    return nearby_germlines


######
# FUS
######


def build_fusion_info(fusion_file, webserver, discarded_fusion_file=None, junction_filter=3, spanning_filter=1):
    if discarded_fusion_file is not None:  # should be if else with discarded read by default, and second option for star-fusion outs
        fus_df = pandas.concat([read_dataset(fusion_file, index_col=None),
                                read_dataset(discarded_fusion_file, index_col=None)], axis=0)
    else:
        fus_df = read_dataset(fusion_file, index_col=None)

    fus_df = fus_df.query('peptide_sequence.notna()', engine='python')  ### warning ###
    fus_df['split_reads'] = fus_df['split_reads1'] + fus_df['split_reads2']

    fus_df = fus_df[((fus_df['discordant_mates'] > spanning_filter) & (fus_df['split_reads'] > junction_filter)) |
                    (fus_df['confidence'].isin(['medium', 'high']))]

    fus_info = []  # empty list
    # Creating named tuple
    Fusion_Info = namedtuple('Fusion_info',
                             ['fusion_id', 'gene_id1', 'gene_id2', 'trans_id1', 'trans_id2',
                              'junction_coverage', 'spanning_coverage', 'reading_frame',
                              'sites', 'breakpoints', 'protein_sequence', 'left_breakpoint_distance', 'symbols',
                              'normal_peptide', 'mismatches'])

    fus_transcript_info = defaultdict(dict)
    protein_positions = defaultdict(lambda: defaultdict(dict))

    for index, row in fus_df.iterrows():

        protein_sequence = row['peptide_sequence']

        if '|' not in protein_sequence:
            continue

        geneID1 = row['gene_id1']
        geneID2 = row['gene_id1']
        transID1 = row['transcript_id1']
        transID2 = row['transcript_id1']
        junction_coverage = row['discordant_mates']
        spanning_coverage = row['split_reads']
        reading_frame = row['reading_frame']
        sites = '{}|{}'.format(row['site1'], row['site2'])
        breakpoints = '{}|{}'.format(row['breakpoint1'], row['breakpoint2'])
        normal_peptide = None
        mismatches = None

        left_breakpoint_distance = ''

        symbols = '{}|{}'.format(row['#gene1'], row['gene2'])

        fusion_id = '{}_{}'.format(symbols, breakpoints)

        fus_info.append(Fusion_Info(fusion_id, geneID1, geneID2, transID1, transID2,
                                    junction_coverage, spanning_coverage, reading_frame,
                                    sites, breakpoints, protein_sequence, left_breakpoint_distance, symbols,
                                    normal_peptide, mismatches))

        fus_transcript_info[fusion_id].setdefault('{}|{}'.format(geneID1, geneID2), []).append([transID1, transID1])

    if not fus_info:
        print_ifnot_webserver('\nNO RELEVANT FUSIONS found', webserver=webserver)

    fus_counter = fus_df['reading_frame'].value_counts().to_dict()

    return fus_info, fus_transcript_info, fus_counter


"""
MuPeX
"""


def reference_peptide_extraction(proteome_reference, peptide_length, tmp_dir, webserver, keep_tmp, file_prefix, outdir,
                                 config_file):
    print_ifnot_webserver('\tExtracting all possible peptides from reference', webserver)
    reference_peptides = set()
    reference_peptide_count = 0
    reference_peptide_file_names = defaultdict(dict)  # empty dictionary

    # loop trough the proteome reference and chop up in all peptides
    for length in peptide_length:
        # check if peptide references are stated in config file
        reference_peptide_file_path = config_parse(config_file, 'References', 'pep' + str(length))
        # Create peptide reference file if not stated in the config file
        if not reference_peptide_file_path is None:
            print_ifnot_webserver('\t\tPeptide reference file of {} aa are stated in config file'.format(length),
                                  webserver)
            reference_peptide_file_names[length] = reference_peptide_file_path
            with open(reference_peptide_file_path) as f:
                for line in f.readlines():
                    reference_peptides.add(line.strip())
                    reference_peptide_count += 1
        else:
            reference_length_peptides = set()
            print_ifnot_webserver('\t\tPeptides of {} aa are being extracted'.format(length), webserver)
            reference_peptides_file = NamedTemporaryFile(delete=False, dir=tmp_dir)
            for geneID in proteome_reference:
                for transID in proteome_reference[geneID]:
                    aa_sequence = proteome_reference[geneID][transID]
                    peptide_list = chopchop(aa_sequence, length)
                    for pep in peptide_list:
                        if not len(pep) == length:
                            continue
                        if '*' in pep:
                            continue
                        if 'X' in pep:
                            continue
                        if 'U' in pep:
                            continue
                        reference_peptide_count += 1
                        reference_length_peptides.add(pep)
                        reference_peptides.add(pep)
            reference_peptides_file.write(str.encode('{}\n'.format('\n'.join(reference_length_peptides))))
            reference_peptides_file.close()
            reference_peptide_file_names[length] = reference_peptides_file

    ReferencePetideCounters = namedtuple('reference_peptide_counters', ['total_peptide_count', 'unique_peptide_count'])
    reference_peptide_counters = ReferencePetideCounters(reference_peptide_count, len(reference_peptides))

    keep_temp_file(keep_tmp, 'txt', reference_peptide_file_names, file_prefix, outdir, peptide_length,
                   'reference_peptide')

    return reference_peptides, reference_peptide_counters, reference_peptide_file_names


def chopchop(aaSeq, peptide_length, mut_type='SNV', reading_frame='in-frame'):
    if mut_type == 'SNV':
        peptides = []
        for i in range(len(aaSeq)):
            pep = aaSeq[i:i + peptide_length]
            if len(pep) < peptide_length:
                break
            peptides.append(pep)

    elif mut_type == 'FUS':

        if reading_frame == 'in-frame':
            pep_range = peptide_length - 1
        elif reading_frame == 'out-of-frame':
            pep_range = len(aaSeq) - 1
        else:
            sys.exit('reading_frame argument must be either in-frame or out-of-frame')

        if '|' in aaSeq:
            peptides = []
            breakpoint_pos = aaSeq.find('|')
            for i in range(pep_range):

                pep_start = (breakpoint_pos - (peptide_length - i - 1))

                if i < peptide_length - 1:
                    pep_end = breakpoint_pos + i + 2
                else:
                    pep_end = breakpoint_pos + i + 1

                fus_seq_i = aaSeq[pep_start:pep_end]

                if len(fus_seq_i.replace('|', '')) != peptide_length:
                    continue

                peptides.append(fus_seq_i)
        else:
            warnings.warn('Fusion protein sequence must have | as a breakpoint symbol')
    else:
        sys.exit('mut_type argument should be either SNV or FUS')

    return peptides


def peptide_extraction(peptide_lengths, vep_info, proteome_reference, genome_reference, reference_peptides,
                       reference_peptide_file_names, fasta_file_name, peptide_match, tmp_dir, webserver, print_mismatch,
                       keep_tmp, file_prefix, outdir, num_mismatches):
    print_ifnot_webserver('\tPeptide extraction begun', webserver)
    peptide_count, normal_match_count, removal_count = 0, 0, 0
    peptide_info = defaultdict(dict)  # empty dictionary
    fasta_printout = defaultdict(dict) if not fasta_file_name is None else None
    pepmatch_file_names = defaultdict(dict)  # empty dictionary

    for p_length in peptide_lengths:
        mutated_peptides_missing_normal = set()
        for mutation_info in vep_info:
            # create mutation counters
            intermediate_peptide_counters = {'mutation_peptide_count': 0, 'mutation_normal_match_count': 0,
                                             'peptide_removal_count': 0}
            # extract sequence
            peptide_sequence_info = mutation_sequence_creation(mutation_info, proteome_reference, genome_reference,
                                                               p_length)
            if peptide_sequence_info is not None:
                if fasta_printout is not None:
                    fasta_printout = long_peptide_fasta_creation(peptide_sequence_info, mutation_info, fasta_printout)
                normpeps, mutpeps = chopchop(peptide_sequence_info.chop_normal_sequence, p_length), chopchop(
                    peptide_sequence_info.mutation_sequence, p_length)
                peptide_mutation_position = peptide_mutation_position_annotation(mutpeps,
                                                                                 peptide_sequence_info.mutation_position,
                                                                                 p_length)
                peptide_info, intermediate_peptide_counters = peptide_selection(normpeps, mutpeps,
                                                                                peptide_mutation_position,
                                                                                intermediate_peptide_counters,
                                                                                peptide_sequence_info, peptide_info,
                                                                                mutation_info, p_length,
                                                                                reference_peptides)
                mutated_peptides_missing_normal = normal_peptide_identification(peptide_info,
                                                                                mutated_peptides_missing_normal,
                                                                                mutpeps, mutation_info)

            # Accumulate counters
            peptide_count += intermediate_peptide_counters['mutation_peptide_count']
            normal_match_count += intermediate_peptide_counters['mutation_normal_match_count']
            removal_count += intermediate_peptide_counters['peptide_removal_count']

        peptide_info, pepmatch_file_names = normal_peptide_correction(mutated_peptides_missing_normal, mutation_info,
                                                                      p_length, reference_peptide_file_names,
                                                                      peptide_info, peptide_match, tmp_dir,
                                                                      pepmatch_file_names, webserver, print_mismatch,
                                                                      num_mismatches)

    # Create and fill counter named-tuple
    PeptideCounters = namedtuple('peptide_counters', ['peptide_count', 'normal_match_count', 'removal_count'])
    peptide_counters = PeptideCounters(peptide_count, normal_match_count, removal_count)

    keep_temp_file(keep_tmp, 'txt', pepmatch_file_names, file_prefix, outdir, peptide_lengths, 'pepmatch')

    return peptide_info, peptide_counters, fasta_printout, pepmatch_file_names


##### VERY BAD DESCISION, BETTER REDO DICT OF DICTS STRUCTURE AND STORE NORM PEPTIDE WITHIN A NAMED TUPLE
def unlist_dict_of_dicts(dicti):
    for key, value in dicti.items():
        for key2, value2 in value.items():

            if isinstance(value2[0], list):
                dicti[key][key2][0] = dicti[key][key2][0][0]
                print('w', key2, value2[0])
    return dicti


def fusion_peptide_extraction(peptide_lengths, fus_info, fasta_file_name, reference_peptide_file_names,
                              tmp_dir, keep_tmp, file_prefix, webserver, num_mismatches,
                              peptide_match, print_mismatch, outdir):
    print_ifnot_webserver('\tFusion Peptide extraction begun', webserver)
    peptide_info = defaultdict(dict)  # empty dictionary
    fasta_printout = defaultdict(dict) if not fasta_file_name is None else None
    pepmatch_file_names = defaultdict(dict)  # empty dictionary

    for p_length in peptide_lengths:
        print_ifnot_webserver('\t\tFusion Peptides of {} aa are being extracted'.format(p_length), webserver)

        for fusion_info_i in fus_info:
            protein_sequence = fusion_info_i.protein_sequence.replace('*', '').replace('?', '').upper()

            peptides = chopchop(aaSeq=protein_sequence, peptide_length=p_length,
                                mut_type='FUS', reading_frame=fusion_info_i.reading_frame)

            for pep in peptides:
                left_breakpoint_dis = pep.find('|')
                pep = pep.replace('|', '')

                fusion_info_i = fusion_info_i._replace(left_breakpoint_distance=left_breakpoint_dis)

                peptide_info[pep] = {pep: fusion_info_i}

                fasta_printout[
                    f'>DTU_FUSION_{fusion_info_i.fusion_id}_length={p_length}_LBD={left_breakpoint_dis}'] = pep

        peptide_info, pepmatch_file_names = normal_peptide_correction(
            mutated_peptides_missing_normal=set(peptide_info.keys()),
            mutation_info=fus_info,
            peptide_length=p_length,
            reference_peptide_file_names=reference_peptide_file_names,
            peptide_info=peptide_info,
            peptide_match=peptide_match, tmp_dir=tmp_dir,
            pepmatch_file_names=pepmatch_file_names,
            webserver=webserver, print_mismatch=print_mismatch,
            num_mismatches=num_mismatches,
            mode='FUS')

    peptide_info = unlist_dict_of_dicts(peptide_info)

    keep_temp_file(keep_tmp, 'txt', pepmatch_file_names, file_prefix, outdir, peptide_lengths, 'pepmatch')
    peptide_counters = len(peptide_info)

    return peptide_info, peptide_counters, fasta_printout, pepmatch_file_names


# peptide_extraction
def mutation_sequence_creation(mutation_info, proteome_reference, genome_reference, peptide_length):
    # Create empty named tuple
    PeptideSequenceInfo = namedtuple('peptide_sequence_info',
                                     ['chop_normal_sequence', 'mutation_sequence', 'normal_sequence',
                                      'mutation_position', 'consequence', 'relative_germline_positions'])

    if mutation_info.mutation_consequence == 'missense_variant':
        peptide_sequence_info = missense_variant_peptide(proteome_reference, mutation_info, peptide_length,
                                                         PeptideSequenceInfo)
    elif mutation_info.mutation_consequence == 'inframe_insertion':
        peptide_sequence_info = insertion_peptide(proteome_reference, mutation_info, peptide_length,
                                                  PeptideSequenceInfo)
    elif mutation_info.mutation_consequence == 'inframe_deletion':
        peptide_sequence_info = deletion_peptide(proteome_reference, mutation_info, peptide_length, PeptideSequenceInfo)
    elif mutation_info.mutation_consequence == 'frameshift_variant':
        peptide_sequence_info = frame_shift_peptide(genome_reference, proteome_reference, mutation_info, peptide_length,
                                                    PeptideSequenceInfo)
    else:
        peptide_sequence_info = None

    return peptide_sequence_info


# peptide_extraction > mutation_sequence_creation
def missense_variant_peptide(proteome_reference, mutation_info, peptide_length, PeptideSequenceInfo):
    asserted_proteome = reference_assertion(proteome_reference, mutation_info, reference_type='proteome')
    aaseq = asserted_proteome[0]
    index = index_creator(mutation_info.prot_pos, peptide_length, aaseq, cdna_pos_start=None, index_type='amino_acid',
                          codon=None, frame_type=None)
    normal_sequence = aaseq[index.lower_index:index.higher_index]
    mutation_sequence = normal_sequence[:index.mutation_peptide_position - 1] + mutation_info.aa_mut + normal_sequence[
                                                                                                       index.mutation_peptide_position:]
    consequence = 'M'
    relative_germline_positions = None ####
    
    if mutation_info.nearby_germlines:
        relative_germline_positions = []
        somatic_pos = int(mutation_info.prot_pos)
        for chrom_pos, info in mutation_info.nearby_germlines.items():
            germline_consequence = info[0]
            germline_dist_to_somatic = abs(somatic_pos - int(info[1]))
            if germline_consequence == 'missense_variant' and germline_dist_to_somatic < peptide_length:
                germline_position = info[1]
                AA_change = info[2].split('/')[1]
                germline_index = int(germline_position) - index.lower_index
                mutation_sequence = mutation_sequence[:germline_index - 1] + AA_change + mutation_sequence[germline_index:]
                relative_germline_positions.append(germline_index)

            elif germline_consequence == 'inframe_insertion':  #### to do
                continue

            elif germline_consequence == 'inframe_deletion': #### to do
                continue

    # Return long peptide (created) and information
    return PeptideSequenceInfo(normal_sequence, mutation_sequence, normal_sequence, index.mutation_peptide_position,
                               consequence, relative_germline_positions)


# peptide_extraction > mutation_sequence_creation
def insertion_peptide(proteome_reference, mutation_info, peptide_length, PeptideSequenceInfo):
    asserted_proteome = reference_assertion(proteome_reference, mutation_info, reference_type='proteome')
    aaseq = asserted_proteome[0]
    index = index_creator(mutation_info.prot_pos, peptide_length, aaseq, cdna_pos_start=None, index_type='amino_acid',
                          codon=None, frame_type=None)
    normal_sequence = aaseq[index.lower_index:index.higher_index]
    relative_germline_positions = None

    if mutation_info.aa_normal == '-':
        mutation_sequence = normal_sequence[:index.mutation_peptide_position] + \
                            mutation_info.aa_mut + normal_sequence[index.mutation_peptide_position:]
    else:
        mutation_sequence = normal_sequence[
                            :index.mutation_peptide_position - 1] + \
                            mutation_info.aa_mut + normal_sequence[index.mutation_peptide_position:]

    chop_normal_sequence = '-' * len(mutation_sequence)
    insertion_range = '{}:{}'.format(index.mutation_peptide_position,
                                     index.mutation_peptide_position + len(mutation_info.aa_mut))
    consequence = 'I'

    # Return long peptide (created) and information
    return PeptideSequenceInfo(chop_normal_sequence, mutation_sequence, normal_sequence, insertion_range, consequence, relative_germline_positions)


# peptide_extraction > mutation_sequence_creation
def deletion_peptide(proteome_reference, mutation_info, peptide_length, PeptideSequenceInfo):
    asserted_proteome = reference_assertion(proteome_reference, mutation_info, reference_type='proteome')
    aaseq = asserted_proteome[0]
    index = index_creator(mutation_info.prot_pos, peptide_length, aaseq, cdna_pos_start=None, index_type='amino_acid',
                          codon=None, frame_type=None)
    normal_sequence = aaseq[index.lower_index:index.higher_index]
    relative_germline_positions = None

    if mutation_info.aa_mut == '-':
        mutation_sequence = normal_sequence[:index.mutation_peptide_position - 1] + normal_sequence[
                                                                                    index.mutation_peptide_position:]
    else:
        mutation_sequence = normal_sequence[
                            :index.mutation_peptide_position - 1] + \
                            mutation_info.aa_mut + normal_sequence[index.mutation_peptide_position +
                                                                   len(mutation_info.aa_mut):]

    chop_normal_sequence = '-' * len(mutation_sequence)
    consequence = 'D'

    # Return long peptide (created) and information
    return PeptideSequenceInfo(chop_normal_sequence, mutation_sequence, normal_sequence,
                               index.mutation_peptide_position - 1, consequence, relative_germline_positions)


# peptide_extraction > mutation_sequence_creation
def frame_shift_peptide(genome_reference, proteome_reference, mutation_info, peptide_length, PeptideSequenceInfo):
    asserted_genome = reference_assertion(genome_reference, mutation_info, reference_type='genome')
    asserted_proteome = reference_assertion(proteome_reference, mutation_info, reference_type='proteome')
    seq = asserted_genome[0]
    cdna_pos_start = asserted_genome[1]
    cdna_pos_end = asserted_genome[2]
    aaseq = asserted_proteome[0]

    aa_index = index_creator(mutation_info.prot_pos, peptide_length, aaseq, cdna_pos_start=None,
                             index_type='amino_acid', codon=None, frame_type=None)
    normal_sequence = aaseq[aa_index.lower_index:aa_index.higher_index]
    relative_germline_positions = None

    if mutation_info.codon_mut.islower():  # frame shift deletion
        n_index = index_creator(mutation_info.prot_pos, peptide_length, aaseq, mutation_info.codon_normal,
                                cdna_pos_start, index_type='nucleotide', frame_type=None)
        mutation_sequence = seq[
                            n_index.lower_index:n_index.cdna_codon_position] + \
                            mutation_info.codon_mut.lower() + seq[n_index.cdna_codon_position + 3:]
    else:  # frame shift insertion
        n_index = index_creator(mutation_info.prot_pos, peptide_length, aaseq, mutation_info.codon_mut, cdna_pos_start,
                                index_type='nucleotide', frame_type='frameshift_insertion')
        new_codon = re.search(r'([A-Z]+)', mutation_info.codon_mut).group(1).strip()
        if mutation_info.codon_normal == '-':
            mutation_sequence = seq[n_index.lower_index - 3:cdna_pos_start] + new_codon.lower() + seq[cdna_pos_end - 1:]
        else:
            mutation_sequence = seq[n_index.lower_index:cdna_pos_start] + new_codon.lower() + seq[cdna_pos_end - 1:]

    detect_stop_codon(mutation_sequence, mutation_info)

    # BIO PYTHON:
    # ignore biopython warnings (as we intentionally create sequences not multiple by three)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', BiopythonWarning)
        # When mutation sequence is obtained, translate the sequence using biopython
        dna_sequence = Seq(mutation_sequence)  # , generic_dna)
        mutation_aaseq = str(dna_sequence.translate(to_stop=True))

    chop_normal_sequence = '-' * len(mutation_aaseq)
    consequence = 'F'
    frameshift_range = '{}:{}'.format(aa_index.mutation_peptide_position, len(mutation_aaseq))

    # Return long peptides and information
    return PeptideSequenceInfo(chop_normal_sequence, mutation_aaseq, normal_sequence, frameshift_range, consequence, relative_germline_positions)


# peptide_extraction > mutation_sequence_creation > frame_shift_peptide
def detect_stop_codon(mutation_sequence, mutation_info):
    # check if the mutation generates a stop-codon
    pos = 1
    for nucleotide in mutation_sequence:
        # find the mutation - annotated with lowercase
        if nucleotide.islower():
            for i in [0, 1, 2]:
                # i identify if codon is within the reading frame
                codon_value = float(pos + i - 3) / 3
                if codon_value.is_integer():
                    codon = mutation_sequence[pos + i - 3: pos + i]
                    # see if the codon is a stop codon
                    if codon.upper() in ['TAA', 'TAG', 'TGA']:
                        print(
                            '\t\tNOTE:\tFrameshift mutation {}/{} in {} is generating the stop codon '
                            '{}\n\t\t\tNo peptides generated from this mutation'.format(
                                mutation_info.aa_normal, mutation_info.aa_mut, mutation_info.trans_id, codon))
        pos = pos + 1


# peptide_extraction
def long_peptide_fasta_creation(peptide_sequence_info, mutation_info, fasta_printout):
    header = '>DTU|{trans_id}_{pos}|{cons}_{mut_change}'.format(cons=peptide_sequence_info.consequence,
                                                                pos=mutation_info.pos,
                                                                mut_change=mutation_info.aa_normal +
                                                                           str(peptide_sequence_info.mutation_position)
                                                                           + mutation_info.aa_mut,
                                                                trans_id=mutation_info.trans_id)
    seq = peptide_sequence_info.mutation_sequence
    fasta_printout[header] = seq

    return fasta_printout


# peptide_extraction
def peptide_mutation_position_annotation(mutpeps, long_peptide_position, peptide_length):
    peptide_mutation_positions = []
    for i in range(len(mutpeps)):
        if type(long_peptide_position) is int:
            pep_mut_pos = long_peptide_position - i
        else:
            start = int(long_peptide_position.split(':')[0])
            end = int(long_peptide_position.split(':')[1])
            final_end = (end - i) if (end - i) < peptide_length else peptide_length
            final_start = (start - i) if (start - i) > 0 else 1
            pep_mut_pos = '{}:{}'.format(int(final_start), int(final_end))
        peptide_mutation_positions.append(pep_mut_pos)
    return peptide_mutation_positions


# peptide_extraction
def peptide_selection(normpeps, mutpeps, peptide_mutation_position, intermediate_peptide_counters,
                      peptide_sequence_info, peptide_info, mutation_info, p_length, reference_peptides):
    if len(peptide_sequence_info.mutation_sequence) < p_length:
        print('\t\tNOTE:\tMutation peptide sequence {} length is smaller than defined peptide length {}'.format(
            peptide_sequence_info.mutation_sequence, p_length))

    for normpep, mutpep, mutpos in zip(normpeps, mutpeps, peptide_mutation_position):
        if '*' in normpep or '*' in mutpep:  # removing peptides from pseudogenes including a *
            intermediate_peptide_counters['peptide_removal_count'] += 1
            continue
        if 'X' in normpep or 'X' in mutpep:
            intermediate_peptide_counters['peptide_removal_count'] += 1
            continue
        if 'U' in normpep or 'U' in mutpep:  # removing peptides with U - non normal AA
            intermediate_peptide_counters['peptide_removal_count'] += 1
            continue
        if mutpep in reference_peptides:  # Counting peptides matching a 100 % normal - Not removing
            intermediate_peptide_counters['mutation_normal_match_count'] += 1
        intermediate_peptide_counters['mutation_peptide_count'] += 1

        pep_match_info = None  # empty variable
        #### check and annotate if the germline mutation,in case its there, is inside of the chopped peptide
        has_germline = 'No'
        germline_positions = []
        if peptide_sequence_info.consequence == 'M':
            start_pos = p_length - int(mutpos) + 1
            end_pos = len(peptide_sequence_info.mutation_sequence)-int(mutpos) + 1
        if peptide_sequence_info.relative_germline_positions is not None:
            if len(peptide_sequence_info.relative_germline_positions) > 0:
                for germ_pos in peptide_sequence_info.relative_germline_positions:
                    if germ_pos >= start_pos and germ_pos <=end_pos:
                        has_germline = 'Yes'
                        germline_positions.append(germ_pos - start_pos + 1)


        # fill dictionary
        peptide_info[mutpep][normpep] = [mutation_info, peptide_sequence_info, mutpos, pep_match_info, has_germline, germline_positions]

    return peptide_info, intermediate_peptide_counters


# peptide_extraction
def normal_peptide_identification(peptide_info, mutated_peptides_missing_normal, mutpeps, mutation_info):
    if not mutation_info.mutation_consequence == 'missense_variant':
        for mutpep in mutpeps:
            if mutpep in peptide_info:  # check that mutated peptide is in the final dictionary of wanted peptides
                mutated_peptides_missing_normal.add(mutpep)
    return mutated_peptides_missing_normal


# peptide_extraction
def normal_peptide_correction(mutated_peptides_missing_normal, mutation_info, peptide_length,
                              reference_peptide_file_names, peptide_info, peptide_match, tmp_dir, pepmatch_file_names,
                              webserver, print_mismatch, num_mismatches, mode='SNV'):
    """
    make 
    """

    # write input file
    mutpeps_file = NamedTemporaryFile(delete=False, dir=tmp_dir)
    mutpeps_file.write(str.encode('{}\n'.format('\n'.join(mutated_peptides_missing_normal))))
    mutpeps_file.close()

    pepmatch_file = run_peptide_match(mutpeps_file, peptide_length, peptide_match, reference_peptide_file_names,
                                      mutation_info, tmp_dir, webserver, print_mismatch, num_mismatches)
    pep_match = build_pepmatch(pepmatch_file, peptide_length, print_mismatch)
    pepmatch_file_names[peptide_length] = pepmatch_file
    
    # insert normal in peptide info
    for mutated_peptide in pep_match:
        
        if not pep_match[mutated_peptide]:
            continue
        else:
            assert mutated_peptide in peptide_info, 'Mutated peptide "{}" not stated in peptide_info data structure'.format(
                    mutated_peptide)
                
            if mode == 'SNV':
                
                #for normal_peptide in peptide_info[mutated_peptide].keys():
                #    # renaming normal key, thereby inserting the normal peptide
                #    peptide_info[mutated_peptide][pep_match[mutated_peptide].normal_peptide] = peptide_info[
                #        mutated_peptide].pop(normal_peptide)
                #    peptide_info[mutated_peptide][pep_match[mutated_peptide].normal_peptide][3] = pep_match[mutated_peptide]
                
                            # Get the new key value from pep_match
                new_key = pep_match[mutated_peptide].normal_peptide

                # Dynamically get the 'invalid' key and replace it directly in peptide_info[mutated_peptide]
                for old_key in peptide_info[mutated_peptide].keys():
                #    if old_key == '---------' or some_other_condition_for_identifying_the_key:
                    # Copy the value associated with the old key to the new key
                    peptide_info[mutated_peptide][new_key] = peptide_info[mutated_peptide][old_key]

                    # Remove the old key
                    del peptide_info[mutated_peptide][old_key]
                    break  # Exit after finding the first matching key
        
            elif mode == 'FUS':    
                
                fusion_info_i = list(peptide_info[mutated_peptide].values())[0]

                #print('mismatches', pep_match[mutated_peptide].mismatch)
                #print('normal pep', pep_match[mutated_peptide].normal_peptide)

                mismatch = pep_match[mutated_peptide].mismatch
                normal_peptide = pep_match[mutated_peptide].normal_peptide

                fusion_info_i = fusion_info_i._replace(mismatches=mismatch) 
                fusion_info_i = fusion_info_i._replace(normal_peptide=normal_peptide)

                peptide_info[mutated_peptide] = {
                    #pep_match[mutated_peptide].normal_peptide:peptide_info[mutated_peptide][mutated_peptide] # very weird bug here!!
                    pep_match[mutated_peptide].normal_peptide: fusion_info_i
                    }
                
    return peptide_info, pepmatch_file_names

# peptide_extraction > normal_peptide_correction
def run_peptide_match(mutpeps_file, peptide_length, peptide_match, reference_peptide_file_names, mutation_info, tmp_dir,
                      webserver, print_mismatch, num_mismatches):
    reference_peptide_file = reference_peptide_file_names[peptide_length]
    print_ifnot_webserver('\t\tRunning {} aa normal peptide match'.format(peptide_length), webserver)
    reference_peptide_file_name = reference_peptide_file.name if not type(
        reference_peptide_file) == str else reference_peptide_file
    pepmatch_file = NamedTemporaryFile(delete=False, dir=tmp_dir)
    if not print_mismatch is None:
        process_pepmatch = subprocess.Popen(
            [peptide_match, '-mm', '-thr', str(num_mismatches), mutpeps_file.name, reference_peptide_file_name],
            stdout=pepmatch_file)
    else:
        process_pepmatch = subprocess.Popen(
            [peptide_match, '-thr', str(num_mismatches), mutpeps_file.name, reference_peptide_file_name],
            stdout=pepmatch_file)
    process_pepmatch.communicate()  # now wait
    pepmatch_file.close()
    return pepmatch_file


# peptide_extraction > normal_peptide_correction
def build_pepmatch(pepmatch_file, peptide_length, print_mismatch):
    pep_match = defaultdict(dict)  # empty dictionary

    with open(pepmatch_file.name) as f:
        for line in f.readlines():
            if re.search(r'^Hit\s', line):
                line = [x.strip() for x in line.split()]
                # save information
                mutated_peptide = line[2]
                normal_peptide = line[3]
                mismatch_peptide = line[4] if print_mismatch is not None else '-' * len(normal_peptide)
                mismatch = line[5] if print_mismatch is None else line[6]
            elif re.search(r'^No Hit found\s', line):
                line = [x.strip() for x in line.split()]
                mutated_peptide = line[3]
                normal_peptide = '-' * len(mutated_peptide)
                mismatch_peptide = '-' * len(mutated_peptide)
                mismatch = 5
            else:
                continue
            # create named tuple
            PepMatchInfo = namedtuple('pep_match_info', ['normal_peptide', 'mismatch', 'mismatch_peptide'])
            pep_match_info = PepMatchInfo(normal_peptide, int(mismatch), mismatch_peptide)
            # fill dictionary
            pep_match[mutated_peptide] = pep_match_info
    return pep_match


# missense_variant_peptide / insertion_peptide / deletion_peptide / frame_shift_peptide
def index_creator(protein_position, peptide_length, aaseq, codon, cdna_pos_start, index_type, frame_type):
    if index_type == 'amino_acid':
        lower_index = max(protein_position - peptide_length, 0)
        higher_index = min(protein_position + (peptide_length - 1), len(aaseq))
        mutation_peptide_position = protein_position - lower_index
        Index = namedtuple('index', ['lower_index', 'higher_index', 'mutation_peptide_position'])
        index = Index(lower_index, higher_index, mutation_peptide_position)
    elif index_type == 'nucleotide':
        codon_position = find_codon_position(codon, frame_type)
        cdna_codon_position = (cdna_pos_start - 1) - codon_position
        cda_transcription_start_position = cdna_codon_position - ((protein_position - 1) * 3)
        lower_index = max(cdna_codon_position - (peptide_length * 3 - 3), cda_transcription_start_position)
        Index = namedtuple('index', ['lower_index', 'cdna_codon_position'])
        index = Index(lower_index, cdna_codon_position)
    return index


# missense_variant_peptide / insertion_peptide / deletion_peptide / frame_shift_peptide > index_creator
def find_codon_position(codon, frame_type):
    for i in range(len(codon)):
        if codon[i].isupper():
            break
    if frame_type == 'frameshift_insertion':
        i -= 1
    return i


# missense_variant_peptide / insertion_peptide / deletion_peptide / frame_shift_peptide
def reference_assertion(reference, mutation_info, reference_type):
    # Check gene id and transcript id exists in proteome_reference dictionary
    assert mutation_info.gene_id in reference, '{gene_id}: This gene is not in the reference gene set. Check that ' \
                                               'you have used the right reference corresponding ' \
                                               'to the one used when running VEP'.format(
        gene_id=mutation_info.gene_id)
    assert mutation_info.trans_id in reference[
        mutation_info.gene_id], '{trans_id}: This transcript is not in the reference gene set. Check that you have ' \
                                'used the right reference corresponding to the one used when running VEP'.format(
        trans_id=mutation_info.trans_id)
    seq = reference[mutation_info.gene_id][mutation_info.trans_id]
    asserted_output = [seq]

    if reference_type == 'proteome':
        # Check mutation is within length of amino acid sequence
        assert len(
            seq) >= mutation_info.prot_pos, 'amino acid sequence length ({}) less than mutation position {}'.format(
            len(seq), mutation_info.prot_pos)

    elif reference_type == 'genome':
        # Check mutation site is within length of cDNA sequence
        cdna_pos_start = int(mutation_info.cdna_pos.split('-')[0])
        cdna_pos_end = int(mutation_info.cdna_pos.split('-')[1]) if '-' in mutation_info.cdna_pos else cdna_pos_start
        assert len(seq) >= cdna_pos_start, 'cDNA sequence length ({}) less than mutation position start {}'.format(
            len(seq), cdna_pos_start)
        assert len(seq) >= cdna_pos_end, 'cDNA sequence length ({}) less ' \
                                         'than mutation position end {}'.format(len(seq), cdna_pos_end)
        asserted_output.extend([cdna_pos_start, cdna_pos_end])

    return asserted_output


# main
def write_fasta(tmp_dir, fasta_printout, webserver):
    if fasta_printout is not None:
        print_ifnot_webserver('\tWriting Fasta file', webserver)
        fasta_file = NamedTemporaryFile(delete=False, dir=tmp_dir)
        for header in fasta_printout:
            fasta_file.write(
                str.encode('{Header}\n{Sequence}\n'.format(Header=header, Sequence=fasta_printout[header])))
        fasta_file.close()
    else:
        fasta_file = None
    return fasta_file


"""
MuPeI
"""


def write_peptide_file(peptide_info, tmp_dir, webserver, keep_tmp, file_prefix, outdir):
    print_ifnot_webserver('\tWriting temporary peptide file', webserver)

    unique_mutant_peptide_count = 0
    peptide_file_names = defaultdict(dict)  # empty dictionary
    peptides = set()

    for mutant_peptide in peptide_info:
        unique_mutant_peptide_count += 1
        peptides.add(mutant_peptide)
        for normal_peptide in peptide_info[mutant_peptide]:
            peptides.add(normal_peptide)

    # write temporary peptide file
    peptide_file = NamedTemporaryFile(delete=False, dir=tmp_dir)
    peptide_file.write(str.encode('{}\n'.format('\n'.join(peptides))))
    peptide_file.close()

    keep_temp_file(keep_tmp, 'txt', peptide_file.name, file_prefix, outdir, None, 'peptide_netMHCinput')

    return unique_mutant_peptide_count, peptide_file


def run_netMHCpan(HLA_alleles, netMHCpan_path, peptide_file, tmp_dir, webserver, keep_tmp, file_prefix, outdir,
                  netmhc_anal):
    # isolate unique HLAalleles
    unique_alleles_set = set(HLA_alleles.split(','))
    unique_alleles = ','.join(map(str, unique_alleles_set))

    netMHCpan_start = datetime.now()

    print_ifnot_webserver(f'\t{file_prefix} Running NetMHCpan eluted ligand prediction', webserver)
    netMHC_EL_file = NamedTemporaryFile(delete=False, dir=tmp_dir)
    process_netMHC_EL = subprocess.Popen([netMHCpan_path, '-p', '-a', unique_alleles, '-f', peptide_file.name],
                                         stdout=netMHC_EL_file)
    process_netMHC_EL.communicate()  # now wait
    netMHC_EL_file.close()
    keep_temp_file(keep_tmp, 'txt', netMHC_EL_file.name, file_prefix, outdir, None, 'netMHCpan_EL')

    # running binding affinity prediction if user specified all analysis
    if netmhc_anal is not None:
        netMHC_BA_file = NamedTemporaryFile(delete=False, dir=tmp_dir)
        print_ifnot_webserver('\tRunning NetMHCpan binding affinity prediction', webserver)
        process_netMHC = subprocess.Popen([netMHCpan_path, '-p', '-a', unique_alleles, '-BA', '-f', peptide_file.name],
                                          stdout=netMHC_BA_file)
        process_netMHC.communicate()  # now wait
        netMHC_BA_file.close()
        keep_temp_file(keep_tmp, 'txt', netMHC_BA_file.name, file_prefix, outdir, None, 'netMHCpan_BA')
    else:
        netMHC_BA_file = None

    netMHCpan_runtime = datetime.now() - netMHCpan_start
    return netMHCpan_runtime, unique_alleles, netMHC_EL_file, netMHC_BA_file


def build_netMHC(netMHC_file, webserver, affinity):
    netmhc_anal = 'binding affinity prediction' if affinity == 'YES' else 'eluted ligand prediction'
    print_ifnot_webserver('\tCreating NetMHCpan {} file dictionary'.format(netmhc_anal), webserver)
    net_mhc = defaultdict(dict)  # empty dictionary
    NetMHCInfo = namedtuple('NetMHCInfo', ['affinity', 'rank', 'score'])

    # Account for different column output between netMHCpan output when running affinity predictions or not
    af = 12 if affinity == 'YES' else None
    r = 13 if affinity == 'YES' else 12

    # Build dictionary
    if netMHC_file is not None:
        with open(netMHC_file.name) as f:
            for line in f.readlines():
                if re.search(r'^\s+1\s', line):
                    line = [x.strip() for x in line.split()]
                    if 'PEPLIST' in line[10]:
                        line[2] = '-' * len(line[2]) if line[2].startswith('XXXX') else line[2]
                        # save information
                        HLA_allele = line[1].replace('*', '')
                        peptide = line[2]
                        affinity = float(line[af]) if af is not None else None
                        rank = float(line[r])
                        score = float(line[11])
                        # fill tuple
                        netmhc_info = NetMHCInfo(affinity, rank, score)
                        # fill dictionary
                        net_mhc[HLA_allele][peptide] = netmhc_info
    else:
        net_mhc = None

    return net_mhc


def extract_snv_info(snv_info_tuple, mutant_peptide, normal_peptide, proteome_reference, protein_positions,
                     transcript_info,
                     snv_qc, cancer_genes, reference_peptides, webserver, print_mismatch, printed_ids):
    mutation_info = snv_info_tuple[0]  # redo it
    peptide_sequence_info = snv_info_tuple[1]

    pep_match_info = snv_info_tuple[3]
    peptide_position = snv_info_tuple[2]
    has_germline = snv_info_tuple[4]
    germline_positions = snv_info_tuple[5]

    mutation_id_vep = '{}_{}_{}/{}'.format(mutation_info.chr, mutation_info.pos, mutation_info.aa_normal,
                                           mutation_info.aa_mut)

    # Extract mismatches
    Mismatches = pep_match_info.mismatch if not pep_match_info is None else 1
    # Print normal peptide or normal peptide only showing mis matched (...X...XX...)
    Norm_peptide = mismatch_snv_normal_peptide_conversion(normal_peptide, peptide_position,
                                                          peptide_sequence_info.consequence,
                                                          pep_match_info) if not print_mismatch is None else normal_peptide
    # Extract transcript IDs including the peptide
    transcript_ids = extract_transcript_ids(mutation_info.gene_id,
                                            transcript_info[mutation_id_vep][mutation_info.gene_id],
                                            proteome_reference, normal_peptide,
                                            peptide_sequence_info.consequence)
    # Extract protein position
    protein_positions_extracted = extract_protein_position(transcript_ids, mutation_id_vep,
                                                           mutation_info.gene_id, protein_positions)

    snv_qcs = state_snv_qc(snv_qc, mutation_info)

    # Extract cancer genes if file is given
    if not cancer_genes is None:
        cancer_gene = 'Yes' if mutation_info.symbol in cancer_genes else 'No'
    else:
        cancer_gene = '-'

    return {**{
        'Norm_peptide': Norm_peptide, 'Gene_ID': mutation_info.gene_id,
        'mutation_id_vep': mutation_id_vep, 'Transcript_ID': mutation_info.trans_id,
        'Amino_Acid_Change': '{}/{}'.format(mutation_info.aa_normal, mutation_info.aa_mut),
        'peptide_position': protein_positions_extracted, 'Chr': mutation_info.chr,
        'Genomic_Position': mutation_info.pos, 'Protein_position': mutation_info.prot_pos,
        'Mutation_Consequence': mutation_info.mutation_consequence,
        'Gene_Symbol': '-' if mutation_info.symbol is None else mutation_info.symbol,
        'Cancer_Driver_Gene': cancer_gene, 'Mismatches': Mismatches,
        'Proteome_Peptide_Match': 'Yes' if mutant_peptide in reference_peptides else 'No'},
         **snv_qcs,
        'Has_germline':has_germline,
        'Germline_positions':germline_positions}


def extract_fusion_info(fusion_info_tuple, cancer_genes, printed_ids):
    fusion_id = fusion_info_tuple.fusion_id
    Gene_ID = fusion_info_tuple.gene_id1
    Gene_ID2 = fusion_info_tuple.gene_id2

    Transcript_ID = fusion_info_tuple.trans_id1
    Transcript_ID2 = fusion_info_tuple.trans_id2

    Chr = re.sub(':.+', '', fusion_info_tuple.breakpoints)
    Chr2 = re.sub('.+\\||:.+[0-9]$', '', fusion_info_tuple.breakpoints)

    Mutation_Consequence = '{}_fusion'.format(fusion_info_tuple.reading_frame)
    sites = fusion_info_tuple.sites
    left_breakpoint_distance = fusion_info_tuple.left_breakpoint_distance

    junction_coverage = fusion_info_tuple.junction_coverage
    spanning_coverage = fusion_info_tuple.spanning_coverage
    reading_frame = fusion_info_tuple.reading_frame
    breakpoints = fusion_info_tuple.breakpoints
    protein_sequence = fusion_info_tuple.protein_sequence
    left_breakpoint_distance = fusion_info_tuple.left_breakpoint_distance

    Gene_Symbols = fusion_info_tuple.symbols
    Gene_Symbol = re.sub(r'\|.+', '', Gene_Symbols)
    Gene_Symbol2 = re.sub(r'.+\|', '', Gene_Symbols)

    normal_peptide = fusion_info_tuple.normal_peptide
    Mismatches = fusion_info_tuple.mismatches

    # Extract cancer genes if file is given
    if not cancer_genes is None:
        cancer_gene = 'Yes' if any([item in cancer_genes for item in [Gene_Symbol, Gene_Symbol2]]) else 'No'
    else:
        cancer_gene = '-'

    return {'fusion_id': fusion_id, 'Gene_ID': Gene_ID, 'Gene_ID2': Gene_ID2,
            'Transcript_ID': Transcript_ID, 'Transcript_ID2': Transcript_ID2,
            'Chr': Chr, 'Chr2': Chr2, 'Mutation_Consequence': Mutation_Consequence,
            'sites': sites, 'left_breakpoint_distance': left_breakpoint_distance,
            'junction_coverage': junction_coverage, 'spanning_coverage': spanning_coverage,
            'breakpoints': breakpoints, 'protein_sequence': protein_sequence, 'Cancer_Driver_Gene': cancer_gene,
            'Mismatches': Mismatches, 'Norm_peptide': normal_peptide,
            'Gene_Symbol': Gene_Symbol, 'Gene_Symbol2': Gene_Symbol2}  # gotta add mismatches here too


def write_output_file(peptide_info, expression, net_mhc_BA, net_mhc_EL,
                      unique_alleles, cancer_genes, tmp_dir,
                      webserver, print_mismatch, snv_qc, expression_file_type, transcript_info,
                      reference_peptides, proteome_reference, protein_positions,
                      version, mode='SNV'):
    print_ifnot_webserver(f'\tWriting output file for {mode}-derived peptides', webserver)
    printed_ids = set()
    row = 0

    mupexi_core_cols = ['HLA_allele', 'Mut_peptide', 'Norm_peptide', 'Mismatches',
                        'Cancer_Driver_Gene', 'Expression_Level', 'Expression_score',
                        'Chr', 'Gene_ID', 'Gene_Symbol', 'Transcript_ID',
                        'Mutation_Consequence', 'Has_germline', 'Germline_positions']  # mismatches should be in the core

    net_mhc_EL_cols = ['Mut_MHCrank_EL', 'Mut_MHCscore_EL', 'Mutant_affinity_score']
    net_mhc_BA_cols = ['Mut_MHCrank_BA', 'Mut_MHCscore_BA', 'Mut_MHCaffinity']

    normal_mhc_EL_cols = ['Norm_MHCrank_EL', 'Norm_MHCscore_EL', 'Normal_affinity_score']
    normal_mhc_BA_cols = ['Norm_MHCrank_BA', 'Norm_MHCscore_BA', 'Norm_MHCaffinity']

    snv_cols = ['mutation_id_vep', 'Amino_Acid_Change', 'peptide_position',
                'tumor_vaf', 'normal_vaf', 't_depth', 'n_depth', 't_alt_count',
                'Genomic_Position', 'Protein_position', 'priority_Score']  # priority score should also be in core cols

    fusion_cols = ['fusion_id', 'Chr2', 'breakpoints',
                   'Gene_ID2', 'Transcript_ID2', 'Gene_Symbol2',
                   'junction_coverage', 'spanning_coverage', 'left_breakpoint_distance', 'sites']

    # Create data frame
    if net_mhc_BA is None:
        if mode == 'SNV':
            cols = mupexi_core_cols + net_mhc_EL_cols + normal_mhc_EL_cols + snv_cols
            df = pandas.DataFrame(columns=cols, )
        elif mode == 'FUS':
            cols = mupexi_core_cols + net_mhc_EL_cols + normal_mhc_EL_cols + fusion_cols
            df = pandas.DataFrame(columns=cols, )
    else:
        if mode == 'SNV':
            cols = mupexi_core_cols + net_mhc_EL_cols + \
                   net_mhc_BA_cols + normal_mhc_BA_cols + normal_mhc_EL_cols + snv_cols
            df = pandas.DataFrame(columns=cols, )
        elif mode == 'FUS':
            cols = mupexi_core_cols + net_mhc_EL_cols + \
                   normal_mhc_EL_cols + normal_mhc_BA_cols + net_mhc_BA_cols + fusion_cols
            df = pandas.DataFrame(columns=cols, )

        # Extract data
        for mutant_peptide in peptide_info:
            for normal_peptide in peptide_info[mutant_peptide]:
                for hla in unique_alleles.split(','):

                    mutant_netmhc_info = net_mhc_EL[hla][mutant_peptide]
                    normal_netmhc_info = net_mhc_EL[hla][normal_peptide]
                    mutant_netmhc_BA_info = net_mhc_BA[hla][mutant_peptide]
                    normal_netmhc_BA_info = net_mhc_BA[hla][normal_peptide]

                    if mode == 'SNV':
                        # Checking concordance between MHC files  and intermediatepeptide_info file
                        assert hla in net_mhc_EL, 'Allele "{}" not stated in NetMHCpan output'.format(hla)
                        assert mutant_peptide in net_mhc_EL[
                            hla], 'Mutant peptide "{}" not found in NetMHCpan output'.format(
                            mutant_peptide)
                        assert mutant_peptide in net_mhc_EL[
                            hla], 'Normal peptide "{}" not found in NetMHCpan output'.format(
                            mutant_peptide)

                        row_content = extract_snv_info(snv_info_tuple=peptide_info[mutant_peptide][normal_peptide],
                                                       mutant_peptide=mutant_peptide, normal_peptide=normal_peptide,
                                                       proteome_reference=proteome_reference,
                                                       protein_positions=protein_positions,
                                                       transcript_info=transcript_info,
                                                       snv_qc=snv_qc,  # start here
                                                       cancer_genes=cancer_genes, reference_peptides=reference_peptides,
                                                       webserver=webserver, print_mismatch=print_mismatch,
                                                       printed_ids=printed_ids
                                                       )

                        # Extract expression value if file is given expression_file_type, expression,
                        # gene_id, webserver, transcripts, printed_ids
                        if expression is not None:
                            Expression_Level = extract_expression_value(expression_file_type=expression_file_type,
                                                                        expression=expression,
                                                                        transcripts=
                                                                        transcript_info[row_content['mutation_id_vep']][
                                                                            row_content['Gene_ID']],
                                                                        gene_id=row_content['Gene_ID'],
                                                                        webserver=webserver, printed_ids=printed_ids)

                            # calculate priority score rank_normal, rank_mutant
                            expression_score, priority_score = score_creation(rank_normal=normal_netmhc_info.rank,
                                                                              rank_mutant=mutant_netmhc_info.rank,
                                                                              expression_sum=Expression_Level,
                                                                              gene_symbol=row_content['Gene_Symbol'],
                                                                              mismatches=row_content['Mismatches'],
                                                                              tumor_vaf=row_content[
                                                                                  'tumor_vaf'],
                                                                              reference_peptides=reference_peptides,
                                                                              mutant_peptide=mutant_peptide)
                        else:
                            Expression_Level = None
                            expression_score = None
                            priority_score = None

                        row_content.update({'HLA_allele': hla, 'Mut_peptide': mutant_peptide,
                                            'Mut_MHCrank_EL': mutant_netmhc_info.rank,
                                            'Mut_MHCscore_EL': mutant_netmhc_info.score,
                                            'Normal_affinity_score': normal_netmhc_info.affinity,
                                            'Mutant_affinity_score': mutant_netmhc_info.affinity,
                                            'Norm_MHCrank_EL': normal_netmhc_info.rank,
                                            'Norm_MHCscore_EL': normal_netmhc_info.score,
                                            'Expression_Level': '-' if Expression_Level is None else Expression_Level,
                                            'Expression_score': expression_score, 'priority_Score': priority_score})

                        if net_mhc_BA is not None:
                            # Checking concordance between MHC files  and intermediatepeptide_info file
                            assert hla in net_mhc_BA, 'Allele "{}" not stated in NetMHCpan output'.format(hla)
                            assert mutant_peptide in net_mhc_BA[
                                hla], 'Mutant peptide "{}" not found in NetMHCpan output'.format(
                                mutant_peptide)
                            assert normal_peptide in net_mhc_BA[
                                hla], 'Normal peptide "{}" not found in NetMHCpan output'.format(normal_peptide)

                            row_content.update({'Mut_MHCrank_BA': mutant_netmhc_BA_info.rank,
                                                'Mut_MHCscore_BA': mutant_netmhc_BA_info.score,
                                                'Mut_MHCaffinity': mutant_netmhc_BA_info.affinity,
                                                'Norm_MHCrank_BA': normal_netmhc_BA_info.rank,
                                                'Norm_MHCscore_BA': normal_netmhc_BA_info.score,
                                                'Norm_MHCaffinity': normal_netmhc_BA_info.affinity})

                    elif mode == 'FUS':
                        # Checking concordance between MHC files  and intermediatepeptide_info file
                        assert hla in net_mhc_EL, 'Allele "{}" not stated in NetMHCpan output'.format(hla)
                        assert mutant_peptide in net_mhc_EL[
                            hla], 'Mutant peptide "{}" not found in NetMHCpan output'.format(
                            mutant_peptide)

                        row_content = extract_fusion_info(
                            fusion_info_tuple=peptide_info[mutant_peptide][normal_peptide],
                            cancer_genes=cancer_genes, printed_ids=printed_ids)

                        if expression is not None:

                            try:
                                Expression_Level = extract_expression_value(expression_file_type=expression_file_type,
                                                                            expression=expression,
                                                                            gene_id=None, webserver=webserver,
                                                                            printed_ids=printed_ids,
                                                                            transcripts=[row_content['Transcript_ID'],
                                                                                         row_content['Transcript_ID2']])
                            except:
                                Expression_Level = '-'
                        else:
                            Expression_Level = None

                        row_content.update({'HLA_allele': hla, 'Mut_peptide': mutant_peptide,
                                            'Mut_MHCrank_EL': mutant_netmhc_info.rank,
                                            'Mut_MHCscore_EL': mutant_netmhc_info.score,
                                            'Mutant_affinity_score': mutant_netmhc_info.affinity,
                                            'Norm_MHCrank_EL': normal_netmhc_info.rank,
                                            'Normal_affinity_score': normal_netmhc_info.affinity,
                                            'Norm_MHCscore_EL': normal_netmhc_info.score,
                                            'Expression_Level': Expression_Level,
                                            'Expression_score': Expression_Level})

                        if net_mhc_BA is not None:
                            # Checking concordance between MHC files  and intermediatepeptide_info file
                            assert hla in net_mhc_BA, 'Allele "{}" not stated in NetMHCpan output'.format(hla)
                            assert mutant_peptide in net_mhc_BA[
                                hla], 'Mutant peptide "{}" not found in NetMHCpan output'.format(
                                mutant_peptide)

                            row_content.update({'Mut_MHCrank_BA': mutant_netmhc_BA_info.rank,
                                                'Mut_MHCscore_BA': mutant_netmhc_BA_info.score,
                                                'Mut_MHCaffinity': mutant_netmhc_BA_info.affinity,
                                                'Norm_MHCrank_BA': normal_netmhc_BA_info.rank,
                                                'Norm_MHCscore_BA': normal_netmhc_BA_info.score,
                                                'Norm_MHCaffinity': normal_netmhc_BA_info.affinity})

                    df.loc[row] = pandas.Series(row_content)[cols]

                row += 1

    # Sort, round up prioritization score
    print_ifnot_webserver('\tSorting output file', webserver)

    if mode == 'SNV':
        df_sorted = df.sort_values('priority_Score',
                                   ascending=False) if not pandas.__version__ == '0.16.0' else df.sort(
            columns='priority_Score', ascending=False)
        # df_sorted.loc[:, 'priority_Score'] = df_sorted.priority_Score.multiply(100).round().astype(int)
        # df_sorted.loc[:, 'Mismatches'] = df_sorted.Mismatches.astype(int)
    else:
        df_sorted = df

    # Print data frame to intermediate file
    df_file = NamedTemporaryFile(delete=False, dir=tmp_dir)
    df_sorted.to_csv(df_file.name, sep='\t', index=False)
    df_file.close()

    # Print header to output file
    header_file = NamedTemporaryFile(delete=False, dir=tmp_dir)
    header = "# Mode:\t{mode}\t\n# VERSION:\tMuPeXI {version}\n# CALL:\t\t{call}\n# DATE:\t\t{day} " \
             "{date} of {month} {year}\n# TIME:\t\t{print_time}\n# PWD:\t\t{pwd}\n"
    header_file.write(str.encode(header.format(mode=mode,
                                               version=version,
                                               call=' '.join(map(str, sys.argv)),
                                               day=datetime.now().strftime("%A"),
                                               month=datetime.now().strftime("%B"),
                                               year=datetime.now().strftime("%Y"),
                                               date=datetime.now().strftime("%d"),
                                               print_time=datetime.now().strftime("%T"),
                                               pwd=os.getcwd())))
    header_file.close()

    # combining header and data frame file
    output_file = NamedTemporaryFile(delete=False, dir=tmp_dir)
    process_cat_files = subprocess.Popen(['cat', header_file.name, df_file.name], stdout=output_file)
    process_cat_files.communicate()
    output_file.close()

    return output_file


def state_tumor_vaf(allele_fractions, mutation_info):
    # check IDs frim allele_fractions dict are in mutation info, forones missing in mutation_info - fill with NA
    # create ID
    allele_fraction_id = '{chromosome}_{genomic_position}_{altered_allele}'.format(
        chromosome=mutation_info.chr,
        genomic_position=mutation_info.pos if not '-' in mutation_info.pos else mutation_info.pos.split('-')[0],
        altered_allele=mutation_info.alt_allele)
    # extract fraction
    if allele_fractions is not None:
        if allele_fraction_id in allele_fractions:
            tumor_vaf = allele_fractions[allele_fraction_id]
        else:
            tumor_vaf = 0
            print('\tAllele frequency for id {} not found; annotated as 0'.format(allele_fraction_id))
    else:
        tumor_vaf = 0

    return tumor_vaf


def state_snv_qc(snv_qc, mutation_info):
    # fills snv_qc info or a single peptide/mutation while writing the output file
    # checks IDs from allele_fractions dict are in mutation info, for ones missing in mutation_info - fill with NA
    # returns a dict {tumor_vaf, :normal_vaf, t_depth, n_depth, t_alt_count}

    snv_qc_id = '{chromosome}_{genomic_position}_{altered_allele}'.format(
        chromosome=mutation_info.chr,
        genomic_position=mutation_info.pos if not '-' in mutation_info.pos else mutation_info.pos.split('-')[0],
        altered_allele=mutation_info.alt_allele)
    # extract fraction
    qc_none = {'tumor_vaf': None, 'normal_vaf': None,
               't_depth': None, 'n_depth': None, 't_alt_count': None}

    if snv_qc is not None:
        if snv_qc_id in snv_qc:
            snv_qcs = snv_qc[snv_qc_id]
        else:
            snv_qcs = qc_none
            print('\tSNV QC is {} not found; annotated as NA')
    else:
        snv_qcs = qc_none

    return snv_qcs


def extract_expression_value(expression_file_type, expression, gene_id, webserver, transcripts, printed_ids):
    expression_sum = None
    if not expression is None:
        if expression_file_type == 'transcript':
            # should I remove .[0-9]+ here or while adding tr id in the beginning?
            transcripts = [re.sub('\\..+', '', x) for x in transcripts]

            for trans_id in transcripts:
                if not trans_id in expression:
                    if not trans_id in printed_ids:
                        print_ifnot_webserver(
                            '\t\t{} not found in expression file\n\t\t\t"-" annotated '
                            'or value not included in expression sum'.format(trans_id), webserver)
                        printed_ids.add(trans_id)
                else:
                    if expression_sum is None:
                        expression_sum = float(expression[trans_id])
                    else:
                        expression_sum += float(expression[trans_id])
        elif expression_file_type == 'gene':
            if gene_id not in expression:
                if gene_id not in printed_ids:
                    print_ifnot_webserver('\t\t{} not found in expression file, "-" annotated'.format(gene_id),
                                          webserver)
                    printed_ids.add(gene_id)
            else:
                expression_sum = expression[gene_id]

    return expression_sum


def extract_transcript_ids(gene_id, trans_ids, proteome_reference, normal_peptide, consequence):
    peptide_trans_ids = []

    # Find transcript id's including the entire peptide.
    if not consequence == 'M':
        peptide_trans_ids = trans_ids
    else:
        for transID in trans_ids:
            aa_sequence = proteome_reference[gene_id][transID]
            if normal_peptide in aa_sequence:
                peptide_trans_ids.append(transID)

    return peptide_trans_ids


def extract_protein_position(transcript_ids, mutation_id_vep, gene_id, protein_positions):
    protein_positions_extracted = []
    for trans_id in transcript_ids:
        protein_positions_extracted.append(protein_positions[mutation_id_vep][gene_id][trans_id])

    return protein_positions_extracted


def score_creation(rank_normal, rank_mutant, expression_sum, gene_symbol, mismatches, tumor_vaf,
                   reference_peptides, mutant_peptide):
    expression = 0 if expression_sum is None else float(expression_sum)
    normal_match = 0 if mutant_peptide in reference_peptides else 1

    # define constants
    k_expression_adjustment = 0.1  #
    k_rank_adjustment = 0.5  #

    affinity_mutant_score = logistic_funtion(rank_mutant)
    affinity_normal_score = logistic_funtion(rank_normal)
    expression_score = hyperbolic_tangent_function(expression + k_expression_adjustment, expression_sum)

    # calculate priority score
    if tumor_vaf is not None and affinity_normal_score is not None:
        if ',' in tumor_vaf:
            tumor_vaf = tumor_vaf.replace(',','.')
        priority_score = affinity_mutant_score * expression_score * (
                1 - ((2 ** - mismatches) * affinity_normal_score)) * float(tumor_vaf) * normal_match
    else:
        priority_score = None
    # Scores = namedtuple('scores', ['affinity_mutant_score', 'affinity_normal_score', 'expression_score'])
    # scores = Scores(affinity_mutant_score, affinity_normal_score, expression_score)

    return expression_score, priority_score


def logistic_funtion(rank_mutant):
    affinity_score = 1 / (1 + math.exp(5 * (float(rank_mutant) - 2)))
    return affinity_score


def hyperbolic_tangent_function(expression, expression_sum):
    expression_score = 1 if expression_sum is None else math.tanh(expression)
    return expression_score


def mismatch_snv_normal_peptide_conversion(normal_peptide, peptide_position, consequence, pep_match_info):
    if consequence == 'M':
        dot_line = '.' * len(normal_peptide)
        normal_mismatch_peptide = dot_line[0:peptide_position - 1] + normal_peptide[peptide_position - 1] + dot_line[
                                                                                                            peptide_position:]
    else:
        normal_mismatch_peptide = pep_match_info.mismatch_peptide
    return normal_mismatch_peptide


def write_log_file(argv, peptide_length, sequence_count, reference_peptide_counters, vep_counters, peptide_counters,
                   start_time_mupex, start_time_mupei, start_time, end_time_mupex, HLAalleles, netMHCpan_runtime,
                   unique_mutant_peptide_count, unique_alleles, tmp_dir, webserver, version):
    print_ifnot_webserver('\tWriting log file\n', webserver)  # gotta make log with fusions
    log_file = NamedTemporaryFile(delete=False, dir=tmp_dir)

    log = """
        # VERSION:  MuPeXI {version}
        # CALL:     {call}
        # DATE:     {day} {date} of {month} {year}
        # TIME:     {print_time}
        # PWD:      {pwd}

        ----------------------------------------------------------------------------------------------------------
                                                        MuPeX
        ----------------------------------------------------------------------------------------------------------

          Reading protein reference file:            Found {sequence_count} sequences, with {reference_peptide_count} {peptide_length}mers
                                                           of which {unique_reference_peptide_count} were unique peptides
          Reading VEP file:                          Found {non_used_mutation_count} irrelevant (synonymous) mutation consequences which were discarded
                                                     Found {misssense_variant_count} missense variant mutation(s) 
                                                           {inframe_insertion_count} insertion(s)
                                                           {inframe_deletion_count} deletion(s)
                                                           {frameshift_variant_count} frameshift variant mutation(s)
                                                     These non-synonymous mutations where found in {gene_count} genes and {transcript_Count} transcripts
          Checking peptides:                         {normal_match_count} peptides matched a normal peptide. 
                                                     {removal_count} peptides included unsupported symbols (e.g. *, X, U) and were discarded 
          Final Result:                              {peptide_count} potential mutant peptides
          MuPeX Runtime:                             {time_mupex}

        ----------------------------------------------------------------------------------------------------------
                                                        MuPeI
        ----------------------------------------------------------------------------------------------------------

          Reading through MuPex file:                Found {peptide_count} peptides of which {unique_mutant_peptide_count} were unique
          Detecting HLA alleles:                     Detected the following {num_of_HLAalleles} HLA alleles:
                                                        {HLAalleles}
                                                        of which {num_of_unique_alleles} were unique
          Running NetMHCpan 4.0:                     Analyzed {num_of_unique_alleles} HLA allele(s)
                                                     NetMHCpan runtime: {netMHCpan_runtime}
          MuPeI Runtime:                             {time_mupei}
          TOTAL Runtime:                             {time}
          """
    log_file.write(str.encode(
        log.format(version=version,
                   call=' '.join(map(str, argv)),
                   sequence_count=sequence_count,
                   reference_peptide_count=reference_peptide_counters.total_peptide_count,
                   unique_reference_peptide_count=reference_peptide_counters.unique_peptide_count,
                   peptide_length=peptide_length,
                   non_used_mutation_count=vep_counters.non_used_mutation_count,
                   misssense_variant_count=vep_counters.misssense_variant_count,
                   inframe_insertion_count=vep_counters.inframe_insertion_count,
                   inframe_deletion_count=vep_counters.inframe_deletion_count,
                   frameshift_variant_count=vep_counters.frameshift_variant_count,
                   gene_count=vep_counters.gene_count,
                   transcript_Count=vep_counters.transcript_Count,
                   normal_match_count=peptide_counters.normal_match_count,
                   removal_count=peptide_counters.removal_count,
                   peptide_count=peptide_counters.peptide_count,
                   time_mupex=end_time_mupex - start_time_mupex,
                   unique_mutant_peptide_count=unique_mutant_peptide_count,
                   num_of_HLAalleles=len(HLAalleles.split(',')),
                   num_of_unique_alleles=len(unique_alleles.split(',')),
                   HLAalleles=HLAalleles,
                   netMHCpan_runtime=netMHCpan_runtime,
                   time_mupei=datetime.now() - start_time_mupei,
                   time=datetime.now() - start_time,
                   day=datetime.now().strftime("%A"),
                   month=datetime.now().strftime("%B"),
                   year=datetime.now().strftime("%Y"),
                   date=datetime.now().strftime("%d"),
                   print_time=datetime.now().strftime("%T"),
                   pwd=os.getcwd()
                   )))
    log_file.close()
    return log_file


def move_output_files(outdir, log_file, logfile_name, fasta_file, fasta_file_name, output_files,
                      webserver, www_tmp_dir):
    # moving output files to user defined dir or webserver dir
    move_to_dir = outdir if webserver is None else www_tmp_dir

    if log_file is not None:
        shutil.move(log_file.name, '{}/{}'.format(move_to_dir, logfile_name))

    for file, output_file_name in output_files.items():
        if file is None:
            continue
        else:
            shutil.move(file.name, '{}/{}'.format(move_to_dir, output_file_name))

    if fasta_file_name is not None:
        shutil.move(fasta_file.name, '{}/{}'.format(move_to_dir, fasta_file_name))

    if webserver is not None:
        os.system('chmod -R a+rwx {}'.format(www_tmp_dir))


def keep_temp_file(keep_temp, file_extension, file_name, file_prefix, move_to_dir, peptide_lengths, filetype):
    if keep_temp is not None:
        if peptide_lengths is None:
            shutil.copy(file_name, '{}/{}_{}.{}'.format(move_to_dir, file_prefix, filetype, file_extension))
        else:
            for p_length in peptide_lengths:
                file_ = file_name[p_length]
                if not type(file_) == str:
                    shutil.copy(file_.name,
                                '{}/{}_{}_{}.{}'.format(move_to_dir, file_prefix, filetype, p_length, file_extension))
                else:
                    continue


def clean_up(tmp_dir):
    shutil.rmtree(tmp_dir)


def webserver_print_output(webserver, www_tmp_dir, output, logfile, fasta_file_name):
    # If webserver copy log and output files to www_tmp_dir
    if webserver is not None:
        dir_name = os.path.basename(os.path.normpath(www_tmp_dir))

        os.system('cat {}/{}'.format(www_tmp_dir, logfile))

        print('\n-------------------------------------------------------------\n')

        print('\nLink to MuPeXI output file <a href="/services/MuPeXI-1.1/tmp/{}/{}">MuPeXI.out</a>\n'.format(dir_name,
                                                                                                              output))
        print('Link to MuPeXI log file <a href="/services/MuPeXI-1.1/tmp/{}/{}">MuPeXI.log</a>\n'.format(dir_name,
                                                                                                         logfile))
        if fasta_file_name is not None:
            print(
                'Link to Fasta file with peptide info '
                '<a href="/usr/opt/www/pub/CBS/services/MuPeXI-1.1/tmp/{}/{}">Fasta file</a>\n'.format(
                    dir_name, fasta_file_name))


def usage():
    usage = """
        MuPeXI - Mutant Peptide Extractor and Informer

        The current version of this program is available from
        https://github.com/ambj/MuPeXI

        MuPeXI.py accepts a VCF file describing somatic mutations as input, and from this 
        derives a set of mutated peptides of specified length(s). These mutated peptides 
        are returned in a table along with various annotations that may be useful for 
        predicting the immunogenicity of each peptide.

        Usage: {call} -v <VCF-file> [options]
                                                                                    DEFAULT

        Required arguments:
        -v, --vcf-file <file>   VCF file of variant calls, preferably from          SampX.vcf
                                MuTect (only SNVs) or MuTect2 (SNVs and indels)
        -z, --fusion_file       Fusion file, tsv file with chimeric proteins        SampX_fusions.tsv
                                Currently only supports output of Arriba    

        Recommended arguments:
        -a, --alleles           HLA alleles, comma separated.                       HLA-A02:01
        -l, --length            Peptide length, given as single number,             9
                                range (9-11) or comma separated (9,10,11).
        -e, --expression-file   Expression file, tab separated
                                ((ENST*/ENSG*) \t mean)
        -E, --expression-type   Are the expression values in the expression         transcript
                                files determined on transcript or gene level.
                                (transcript/gene)

        Optional arguments affecting output files:

        -o, --output-file       Output file name.                                   <VEP-file>.mupexi
        -d, --out-dir           Output directory - full path.                       current directory
        -p, --prefix            Prefix for output files - will be applied           <VCF-file>
                                to all (.mupexi, .log, .fasta) unless specified 
                                by -o or -L.
        -L, --log-file          Logfile name.                                       <VCF-file>.log
        -m, --mismatch-number   Maximum number of mismatches to search for in       4
                                normal peptide match.
        -A, --assembly          The assembly version to run VEP.                    GRCh38
        -s, --species           Species to analyze (human / mouse / mouse_black6 /  human
                                mouse_balbc)
                                If mouse is set default assembly is GRCm38.
                                Remember to download the correct VEP cache
                                and state species corresponding MHC alleles.

        Optional arguments affecting computational process:
        -F, --fork              Number of processors running VEP.                   2
        --vcf-type              Input VCF caller style (mutect2/merged)            mutect2
        --vcf_type              Alias for --vcf-type
        --germlines             Include germline context variants (true/false)      false
        --rna-edit              Include RNA_EDIT variants (true/false)              false
        --rna_edit              Alias for --rna-edit
        --tumor-sample          Explicit tumor sample column name in VCF
        --normal-sample         Explicit normal sample column name in VCF
        --rnaedit-known-only    Keep only RNA_EDIT variants marked as known         true
        --rnaedit-allow-novel   Include novel RNA_EDIT variants (overrides known-only)
        --rnaedit-known-key     INFO key used to mark known RNA edits               KNOWN_RNAEDIT_DB
        Other options (these do not take values)
        -f, --make-fasta        Create FASTA file with long peptides 
                                - mutation in the middle
        -c, --config-file       Path to the config.ini file                         current directory
        -t, --keep-temp         Retain all temporary files
        -M, --mismatch-only     Print only mismatches in normal peptide sequence 
                                and otherwise use dots (...AA.....)
        -w, --webface           Run in webserver mode
        -g, --liftover          Perform liftover HG19 to GRCh38.
                                Requires local picard installation with paths
                                stated in the config file
        -n, --netmhc-full-anal  Run NetMHCpan 4.0 with the full analysis including 
                                both eluted ligand (EL) and binding affinity (BA) 
                                prediction output (priority score calculated from 
                                EL rank score)
        -h, --help              Print this help information

        REMEMBER to state references in the config.ini file
        """
    print((usage.format(call=sys.argv[0], path='/'.join(sys.argv[0].split('/')[0:-1]))))


def read_options(argv):
    def parse_bool_opt(value, name):
        v = value.strip().lower()
        if v in ('true', 't', '1', 'yes', 'y'):
            return True
        if v in ('false', 'f', '0', 'no', 'n'):
            return False
        sys.exit('ERROR: {} must be true/false'.format(name))

    try:
        optlist, args = getopt.getopt(argv,
                                      'v:z:a:l:o:d:L:e:c:p:E:m:A:F:s:nftMwgh',
                                      ['input-file=', 'fusion-file=', 'alleles=', 'length=', 'output-file=', 'out-dir=',
                                       'log-file=',
                                       'expression-file=', 'config-file=', 'prefix=', 'expression-type=',
                                       'mismatch-number=', 'assembly=', 'fork=', 'species=', 'netmhc-full-anal',
                                       'make-fasta', 'keep-temp', 'mismatch-print', 'webserver', 'liftover', 'help',
                                       'vcf-type=', 'vcf_type=', 'germlines=', 'rna-edit=', 'rna_edit=',
                                       'tumor-sample=', 'normal-sample=', 'rnaedit-known-only',
                                       'rnaedit-allow-novel', 'rnaedit-known-key='])
        if not optlist:
            print('No options supplied')
            usage()
    except getopt.GetoptError as e:
        usage()
        sys.exit(e)

    # Create dictionary of long and short formats
    format_dict = {
        '-h': '--help',
        '-v': '--vcf-file',
        '-z': '--fusion-file',
        '-l': '--length',
        '-e': '--expression-file',
        '-w': '--webface',
        '-p': '--prefix',
        '-o': '--output-file',
        '-d': '--out-dir',
        '-L': '--log-file',
        '-c': '--config-file',
        '-f': '--make-fasta',
        '-t': '--keep-temp',
        '-a': '--alleles',
        '-m': '--mismatch-number',
        '-M': '--mismatch-only',
        '-g': '--liftover',
        '-E': '--expression-type',
        '-A': '--assembly',
        '-F': '--fork',
        '-s': '--species',
        '-n': '--netmhc-full-anal'
    }

    # Create a dictionary of options and input from the options list
    opts = dict(optlist)

    # Use the long format dictionary to change the option to the short annotation, if long is given by the user.
    for short, long_ in six.iteritems(format_dict):
        if long_ in opts:
            opts[short] = opts.pop(long_)

    # Print usage help
    if '-h' in list(opts.keys()):
        usage()
        sys.exit()

    # Define values
    vcf_file = opts['-v'] if '-v' in list(opts.keys()) else None
    fusion_file = opts['-z'] if '-z' in list(opts.keys()) else None

    if vcf_file is None and fusion_file is None:
        sys.exit('At least vcf file or fusion file must be provided, terminating...')
        usage()
    elif fusion_file is None:
        print('Fusion file is missing, chimeric junction antigens will not be calculated')
    elif vcf_file is None:
        print('VCF file is missing, SNV antigens will not be calculated')
    else:
        print('Calculating neoantigens for both chimeric funtions and SNV')

    peptide_length = opts['-l'] if '-l' in list(opts.keys()) else 9
    expression_file = opts['-e'] if '-e' in list(opts.keys()) else None
    expression_type = opts['-E'] if '-E' in list(opts.keys()) else 'transcript'
    webserver = opts['-w'] if '-w' in list(opts.keys()) else None
    prefix = opts['-p'] if '-p' in list(opts.keys()) else vcf_file.split('/')[-1].split('.')[0]
    output = opts['-o'] if '-o' in list(opts.keys()) else prefix + '.mupexi'
    outdir = opts['-d'] if '-d' in list(opts.keys()) else os.getcwd()
    logfile = opts['-L'] if '-L' in list(opts.keys()) else prefix + '.log'
    script_dir = os.path.dirname(os.path.abspath(__file__))
    config = opts['-c'] if '-c' in list(opts.keys()) else os.path.join(script_dir, 'config.ini')
    fasta_file_name = prefix + '.fasta' if '-f' in list(opts.keys()) else None
    keep_temp = 'yes' if '-t' in list(opts.keys()) else None
    HLA_alleles = opts['-a'] if '-a' in list(opts.keys()) else 'HLA-A02:01'
    print_mismatch = 'Yes' if '-M' in list(opts.keys()) else None
    liftover = 'Yes' if '-g' in list(opts.keys()) else None
    num_mismatches = opts['-m'] if '-m' in list(opts.keys()) else 4
    assembly = opts['-A'] if '-A' in list(opts.keys()) else None
    fork = opts['-F'] if '-F' in list(opts.keys()) else 2
    if int(fork) <= 1:
        usage()
        sys.exit('VEP fork number must be greater than 1')
    species = opts['-s'] if '-s' in list(opts.keys()) else 'human'
    netmhc_anal = 'yes' if '-n' in list(opts.keys()) else None

    vcf_type = opts.get('--vcf-type', opts.get('--vcf_type', 'mutect2')).strip().lower()
    if vcf_type not in ('mutect2', 'merged'):
        sys.exit('ERROR: --vcf-type must be one of: mutect2, merged')
    germlines = parse_bool_opt(opts.get('--germlines', 'false'), '--germlines')
    rna_edit = parse_bool_opt(opts.get('--rna-edit', opts.get('--rna_edit', 'false')), '--rna-edit')
    tumor_sample = opts.get('--tumor-sample', None)
    normal_sample = opts.get('--normal-sample', None)
    rnaedit_known_only = True
    if '--rnaedit-known-only' in opts:
        rnaedit_known_only = True
    if '--rnaedit-allow-novel' in opts:
        rnaedit_known_only = False
    rnaedit_known_key = opts.get('--rnaedit-known-key', 'KNOWN_RNAEDIT_DB').strip()
    if not rnaedit_known_key:
        sys.exit('ERROR: --rnaedit-known-key must be non-empty')

    # Create and fill input named-tuple
    Input = namedtuple('input',
                       ['vcf_file', 'fusion_file', 'peptide_length', 'output', 'logfile', 'HLA_alleles', 'config',
                        'expression_file',
                        'fasta_file_name', 'webserver', 'outdir', 'keep_temp', 'prefix', 'print_mismatch', 'liftover',
                        'expression_type', 'num_mismatches', 'assembly', 'fork', 'species', 'netmhc_anal',
                        'vcf_type', 'germlines', 'rna_edit', 'tumor_sample', 'normal_sample',
                        'rnaedit_known_only', 'rnaedit_known_key'])
    inputinfo = Input(vcf_file, fusion_file, peptide_length, output, logfile, HLA_alleles, config, expression_file,
                      fasta_file_name,
                      webserver, outdir, keep_temp, prefix, print_mismatch, liftover, expression_type, num_mismatches,
                      assembly, fork, species, netmhc_anal, vcf_type, germlines, rna_edit, tumor_sample, normal_sample,
                      rnaedit_known_only, rnaedit_known_key)

    return inputinfo


################
# Run analysis #
################

if __name__ == '__main__':
    main(sys.argv[1:])

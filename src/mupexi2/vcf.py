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


def parse_genotype_state(gt_value, ps_value):
    gt = (gt_value or '').strip()
    ps = (ps_value or '').strip()
    if gt in ('.', './.', '.|.', ''):
        return {'gt': gt, 'is_hom_alt': False, 'is_het': False, 'has_pipe': False, 'has_ps': False,
                'is_phased': False, 'hap': None, 'ps': None}

    is_hom_alt = gt in ('1/1', '1|1')
    is_het = gt in ('0/1', '1/0', '0|1', '1|0')
    has_pipe = '|' in gt
    has_ps = ps not in ('', '.')
    is_phased = has_pipe and has_ps
    hap = None
    if is_phased:
        left, right = gt.split('|', 1)
        if left == '1' and right == '0':
            hap = 0
        elif left == '0' and right == '1':
            hap = 1
        elif left == '1' and right == '1':
            hap = 'both'
    return {'gt': gt, 'is_hom_alt': is_hom_alt, 'is_het': is_het, 'has_pipe': has_pipe, 'has_ps': has_ps,
            'is_phased': is_phased, 'hap': hap, 'ps': ps if has_ps else None}


def should_apply_context_variant(primary_state, context_state):
    if context_state is None:
        return False
    if context_state['is_hom_alt']:
        return True
    if primary_state is None:
        return False
    if not primary_state['is_phased'] or not context_state['is_phased']:
        return False
    if primary_state['ps'] is None or context_state['ps'] is None:
        return False
    return primary_state['ps'] == context_state['ps'] and primary_state['hap'] == context_state['hap']


def resolve_sample_indices(header, tumor_sample=None, normal_sample=None):
    sample_names = header[9:]
    sample_index = {name: i + 9 for i, name in enumerate(sample_names)}

    def find_unique(candidates):
        hits = [name for name in sample_names if candidates(name)]
        return hits[0] if len(hits) == 1 else None

    if tumor_sample is not None and tumor_sample not in sample_index:
        raise ValueError('tumor sample "{}" is not present in VCF header'.format(tumor_sample))
    if normal_sample is not None and normal_sample not in sample_index:
        raise ValueError('normal sample "{}" is not present in VCF header'.format(normal_sample))

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
            raise ValueError('VCF must contain at least two sample columns for SNV QC extraction')
        tumor_name = sample_names[0]
        normal_name = sample_names[1]

    return sample_index[tumor_name], sample_index[normal_name], tumor_name, normal_name

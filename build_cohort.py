#!/usr/bin/env python
"""
Assist user in selecting a set of cases and controls that meet certain user
specified criteria for downstream analysis
"""

import argparse
import logging
import MySQLdb
import sys
from base64 import b64decode
from collections import OrderedDict
from operator import le, lt
from functools import partial

def valid_numerical_argument(
    arg, arg_name, arg_type=int, min_value=0, max_value=sys.maxint,
    left_op=lt, right_op=le):
    """Confirm that the specified value is valid in the range
    (minimum_value, maximum_value] (by default)
    :param arg: the value to be tested
    :param arg_name: the name of the parameter
    :param arg_type: the type of the parameter, e.g. int or float
    :param min_value: the minimum value for the parameter, exclusive
    :param max_value: the maximum value for the parameter, inclusive
    :param left_op: the operator for testing left_op(min_value, value)
    :param right_op: the operator testing right_op(value, max_value)
    :return: arg_type(arg) if arg is valid
    """
    try:
        value = arg_type(arg)
        if left_op(min_value, value) and right_op(value, max_value):
            return value
        else:
            raise argparse.ArgumentTypeError(
                "{arg_name} ({arg}) is not in the range "
                "{left_endpoint}{min_value}, {max_value}{right_endpoint}".format(
                    arg_name=arg_name, arg=arg, min_value=min_value,
                    max_value=max_value, left_endpoint="(" if left_op == lt else "[",
                    right_endpoint="]" if right_op == le else ")"))
    except TypeError:
        raise argparse.ArgumentTypeError(
            "{arg_name} ({arg}) is not a valid {arg_type}".format(
                arg_name=arg_name, arg=arg, arg_type=arg_type.__name__))

# generic query format
QUERY = """
        SELECT {selection}{count_clause} 
        FROM SampleT AS SAMP
        INNER JOIN dragen_sample_metadata AS META
        ON META.sample_name = SAMP.CHGVID AND
        SAMP.SeqType = META.sample_type
        INNER JOIN dragen_qc_metrics AS QC
        ON QC.pseudo_prepid = META.pseudo_prepid 
        INNER JOIN dragen_pipeline_step AS STEP
        ON STEP.pseudo_prepid = META.pseudo_prepid
        WHERE STEP.pipeline_step_id = 108 AND
        STEP.step_status = 'completed' AND
        SAMP.BroadPhenotype <> "" AND
        META.sample_type = "{sequencing_type}"
        {univ_condition}
        {seq_spec_condition}
        {avail_condition}
        {pheno_condition}
        GROUP BY {selection}
        """
        #INNER JOIN seqdbClone AS CLN 
        #ON CLN.CHGVID = QC.CHGVID AND
        #CLN.SeqType = QC.SeqType AND
        #CLN.pseudo_prepid = META.pseudo_prepid
        #CLN.Platform <> "GAIIx" AND

DEFAULT_PREF = ['Exome', 'Genome', 'Custom_Capture'] # default SeqType preference order

COUNT_CLAUSE = ', COUNT(*) AS C' # if counts are needed in SQL query

AVAIL_CONDITION = ' AND SAMP.AvaiContUsed = "yes"' # for controls only
d = 'c2VxdWVuY2VEQg=='
h = 'MTAuNzMuNTAuMzg='
u = 'Y29ob3J0X3NlbGVjdGlvbg=='
p = 'Y29ob3J0X3NlbGVjdGlvbg=='

def run_query(cur, selection, count, seq, univ, seq_spec, avail, pheno):
    """
    Build and run query with given parameters:
    cur = db cursor
    selection = what to select from db
    count = clause to count the db selection
    seq = clause to choose certain SeqTypes
    univ = universal conditions for all queries
    seq_spec = SeqType-specific conditions
    avail = control use availability condition
    pheno = clause to choose certain phenotypes
    """
    
    query = QUERY.format(selection=selection,
                         count_clause=count,
                         sequencing_type=seq,
                         univ_condition=univ,
                         seq_spec_condition=seq_spec,
                         avail_condition=avail,
                         pheno_condition=pheno)
    #print query
    cur.execute(query)
    rows = cur.fetchall()
    return rows

def make_seq_spec_condition(seq_type, cov, dbsnp, contam, kits):
    """
    Build part of query that deals with fields specific to different
    sequencing types:
    seq_type = sequencing type
    cov = coverage condition (CCDSBasesCov10x)
    contam = contamination condition (PercentContamination)
    kits = list of selected capture kits (ExomeSamPrepKit)
    """
    
    # BasesCov5x or ExomeSamPrepKit do not apply to genomes
    if seq_type == 'Genome':
        cov_condition = ''
        kit_condition = ''
    # if Exome or Custom_Capture, BasesCov5x and ExomeSamPrepKit do apply
    else:
        cov_condition = ' AND QC.CCDSBasesCov10x >= {cov}'.format(cov=100*cov)
        dbsnp_condition = 'AND QC.DBSNPOverlapSNVs >= {dbsnp}'.format(dbsnp=dbsnp)
        if kits:
            kit_condition = (' AND META.capture_kit IN ("' +
                             '", "'.join([kit for kit in kits]) + '")')
        else: 
            kit_condition = ''
    
    # PercentContamination does not apply to Custom_Capture
    if seq_type == 'Custom_Capture':
        contam_condition = ''
    else:
        contam_condition = ' AND QC.PercentContamination <= ' + str(contam)
    
    return (cov_condition + dbsnp_condition + contam_condition + kit_condition)

def make_pheno_condition(phenos):
    """
    Build part of query that deals with case/control phenotype lists:
    phenos = list of phenotypes
    """
    
    condition = (' AND SAMP.BroadPhenotype IN ("' + 
                 '", "'.join(phenos) + '")')
    
    return condition

def get_input(message, l):
    """
    Prompt user to provide indices for selection of certain sample 
    categories and match those indices with their categories:
    message = prompt that appears on user command line
    l = list of categories and counts sorted alphabetically by category
    """
    
    while True:
        error = False # error indicator
        reply = raw_input(message)
        if reply == 'NA': # interpret 'NA' as blank list
            inds = []
        elif reply == 'all': # interpret 'all' as complete list
            inds = range(len(l))
            #inds = map(str, inds)
        else: # otherwise make list of selected indices
            inds = []
            for ind in reply.replace(' ', '').strip(',').split(','):
                try:
                    if '-' in ind:
                        start, end = [int(value) for value in ind.split('-')]
                    else:
                        start, end = int(ind), int(ind)
                    for x in xrange(start, end + 1):
                        inds.append(x)
                # handle ValueError during int() coercion to avoid program failure
                except ValueError:
                    print('\nYou have made an invalid index selection.\n')
                    error = True
                    break
        # determine if ValueError occurred and restart process if so
        if error:
            continue
        # list of categories (phenotypes or capture kits) with selected indices
        sel_keys = []
        # match selected indices with categories
        for i, item in enumerate(l):
            if i in inds:
                sel_keys.append(item[0])
        # restart process if invalid index (or no index) has been entered
        if len(sel_keys) != len(inds):
            print('\nYou have made an invalid index selection.\n')
            continue
        return sel_keys

def prune_samples(v, seq_prefs, kit_counter):
    """
    If duplicate CHGVIDs, select one from sample list using following criteria
    1. Preferred sequencing types (set by user)
    2. If multiple of same seq type, prefer most common capture kit
    v = list of metadata (seq type, capture kit, gender, prepID, idnum) associated with
        every record under a particular CHGVID
    seq_prefs = list of seq types in order of preference
    kit_counter = dict of capture kits, ordered  by abundance in selected group of samples
    """
    
    # if duplicate CHGVIDs, choose preferred SeqType
    for seq_pref in seq_prefs:
        seq_indices = [i for i, x in enumerate(v) if x[0] == seq_pref] 
        if len(seq_indices) == 0:
            continue
        elif len(seq_indices) == 1:
            return (v[seq_indices[0]])
        else:
            # if duplicate CHGVIDs and SeqTypes, choose most abundant capture kit
            for kit_pref in kit_counter.keys():
                kit_indices = [j for j, y in 
                    enumerate(v) if y[1] == kit_pref]
                if len(kit_indices) == 0:
                    continue
                else:
                    return (v[kit_indices[0]])
    
def build_cohort(sample_file_name, sequencing_types, sequencing_preference, 
                 include_sex_typing_mismatches, min_bases_covered_10x, 
                 dbsnp_overlap_snv, max_contamination, ethnicities,
                 min_genotyping_rate, log_file_name, qc_file_name):
    db = MySQLdb.connect(user=b64decode(u), passwd=b64decode(p),
                         db=b64decode(d), host=b64decode(h))
    try:
        cur = db.cursor()
        
        # sex mismatch condition
        sex_condition = ' AND QC.SeqGender IN ("M", "F")'
        if not include_sex_typing_mismatches:
            sex_condition += ' AND QC.SeqGender = SAMP.SelfDeclGender'
        
        # ethnicity probability and genotyping rate conditions
        if not ethnicities:
            eth_condition = ''
            genotyping_condition = ''
        else:
            genotyping_condition = ' AND QC.genotyping_rate >= ' + str(
                                                   min_genotyping_rate)
            # consider 0.8 probability to be member of certain ethnic group
            eth_condition = (' AND (' + 
                             ' OR '.join(['QC.{eth}_prob >= {pr}'.format(eth=eth,pr=0.8) 
                                    for eth in ethnicities]) + 
                             ')')
        
        # combine three above conditions since they apply to all queries
        univ_condition = sex_condition + eth_condition + genotyping_condition
        
        # print debugging info to log file if option is set
        if log_file_name:
            logging.basicConfig(filename=log_file_name, level=logging.INFO)
            logging.info('\nSequencing type: {}\n'.format(sequencing_types) +
                         'Min. frac. of bases covered at 10x: {}\n'.format(
                                                       min_bases_covered_10x) +
                         'Max. contamination: {}\n'.format(
                                                          max_contamination) +
                         'Include sex typing mismatches: {}\n'.format(
                                              include_sex_typing_mismatches) +
                         'Ethnicity filter: {}\n'.format(
                                     ethnicities if ethnicities else 'None') +
                         'Min. genotyping rate: {}\n'.format(
                               min_genotyping_rate if ethnicities else 'NA') +
                         'Sequencing type preference: {}\n'.format(
                                sequencing_preference))
        
        d_pheno = {} # phenotype-keyed dictionary for counting samples by pheno
        # query database for phenotypes and counts
        for i, seq in enumerate(sequencing_types):
            # build SeqType-specific condition
            seq_spec_condition = make_seq_spec_condition(seq_type=seq, 
                                             cov=min_bases_covered_10x,
                                             dbsnp=dbsnp_overlap_snv, 
                                             contam=max_contamination, 
                                             kits=None)
            rows = run_query(cur=cur, 
                             selection='SAMP.BroadPhenotype', 
                             count=COUNT_CLAUSE,
                             seq=seq, 
                             univ=univ_condition, 
                             seq_spec=seq_spec_condition, 
                             avail='', 
                             pheno='')
            # record counts of each phenotype
            if rows:
                for phenotype, c in rows:
                    if phenotype not in d_pheno:
                        d_pheno[phenotype] = [0]*len(sequencing_types)
                    d_pheno[phenotype][i] = c
            else:
                print('\nThere are no samples that meet the provided critera.\n')
                sys.exit()
        # sort alphabetically and convert to list
        l_pheno = sorted(d_pheno.items())
        # determine longest phenotype name for purposes of output formatting
        max_pheno_len = max([len(pheno[0]) for pheno in l_pheno])
        width = max_pheno_len + 3 if max_pheno_len > 12 else 15 
        
        # display results
        print('\n{index:<10}{phenotype:<{width}}'.format(index='Index',
                                                         phenotype='Phenotype',
                                                         width=width)
            + ''.join(['{seq:<10}'.format(seq=seq) for seq in sequencing_types]))
        for i, pheno in enumerate(l_pheno):
            print('{index:<10}{phenotype:<{width}}'.format(index=i,
                                                           phenotype=pheno[0],
                                                           width=width)
                + ''.join(['{c:<10}'.format(c=count) for count in pheno[1]]))
        
        while True:
            # get user selection of cases and controls, if any
            cases = get_input('Please enter a comma-separated list of one or '
                              'more CASE indices, or "NA" for none: ', 
                              l_pheno)
            controls = get_input('Please enter a comma-separated list of one or '
                                 'more CONTROL indices, or "NA" for none: ', 
                                 l_pheno)
            # do not allow selection of 0 phenotypes
            if len(cases) == 0 and len(controls) == 0:
                print('\nYou must select at least phenotype.\n')
                continue
            # do not allow selection of overlapping phenotypes
            case_set = set(cases)
            control_set = set(controls)
            if len(case_set.intersection(control_set)) != 0:
                print('\nYou cannot select overlapping CASE and CONTROL phenotypes.\n')
                continue
            break
       
        # record selected phenotypes 
        if log_file_name:
            logging.info('\nCase phenotypes: {}\n'.format(cases) +
                         'Control phenotypes: {}\n'.format(controls))
             
        # make list of phenotype conditions for query
        pheno_conditions = [make_pheno_condition(cases), 
                            make_pheno_condition(controls)] 
        
        d_cap = {} # capture-kit-keyed dictinary for counting samples by kit
        
        # query database for capture kits and counts
        for seq in sequencing_types:
            # build SeqType-specific condition
            seq_spec_condition = make_seq_spec_condition(seq_type=seq, 
                                             cov=min_bases_covered_10x,
                                             dbsnp=dbsnp_overlap_snv, 
                                             contam=max_contamination, 
                                             kits=None)
            for i, pheno_condit in enumerate(pheno_conditions):
                rows = run_query(cur=cur, 
                                 selection='META.capture_kit', 
                                 count=COUNT_CLAUSE,
                                 seq=seq, 
                                 univ=univ_condition, 
                                 seq_spec=seq_spec_condition, 
                                 avail='' if i == 0 else AVAIL_CONDITION,
                                 pheno=pheno_condit)
                # record count of each capture kit
                for prep, c in rows:
                    if prep not in d_cap:
                        d_cap[prep] = [0] * 2
                    d_cap[prep][i] += c
        # sort alphabeticaly and convert to list
        l_cap = sorted(d_cap.items())
        # determine longest capture kit name for purposes of output formatting
        max_kit_len = max([len(cap[0]) for cap in l_cap])
        width = max_kit_len + 3 if max_kit_len > 12 else 15 
        
        # display results
        print('\n{index:<10}{capture_kit:<{width}}{cases:<10}{controls:<10}'.format(
              index='Index',capture_kit='Capture Kit',
              width=width,cases='Cases',controls='Controls'))
        for i, cap in enumerate(l_cap):
            print('{index:<10}{kit:<{width}}'.format(
                                        index=i,
                                        kit=cap[0] if cap[0] != 'N/A' else 'N/A (Genome)',
                                        width=width)
                + ''.join(['{c:<10}'.format(c=count) for count in cap[1]]))
        
        while True:
            # get user selection of capture kits
            kits = get_input('Please enter a comma-separated list of CAPTURE '
                             'KIT indices, or "all" for all: ', l_cap)
            # check for failure to select any kits
            if len(kits) == 0:
                print('\nYou must select at least one capture kit.\n')
                continue
            break    
       
        # record selected capture kits 
        if log_file_name:
            logging.info('\nCapture kits: {}\n'.format(kits))
        
        d_samp = {} # CHGVID-keyed dictionary
        
        # final query for samples with all selected criteria
        for seq in sequencing_types:
            # build SeqType-specific condition
            seq_spec_condition = make_seq_spec_condition(seq_type=seq, 
                                             cov=min_bases_covered_10x,
                                             dbsnp=dbsnp_overlap_snv, 
                                             contam=max_contamination, 
                                             kits=kits)
            for i, pheno_condit in enumerate(pheno_conditions):
                rows = run_query(cur=cur, 
                                 selection='SAMP.CHGVID, QC.seqGender, META.capture_kit, '
                                           'META.pseudo_prepid', 
                                 count='', 
                                 seq=seq, 
                                 univ=univ_condition, 
                                 seq_spec=seq_spec_condition, 
                                 avail='' if i == 0 else AVAIL_CONDITION, 
                                 pheno=pheno_condit)
                for chgvid, gender, kit, prepid in rows:
                    if chgvid not in d_samp:
                        d_samp[chgvid] = []
                    d_samp[chgvid].append((seq,kit,gender,
                                           2 if i == 0 else 1,
                                           prepid))
        
        # organize seq types by preference for purpose of breaking ties
        # for duplicate CHGVIDs
        DEFAULT_PREF.remove(sequencing_preference)
        seq_prefs = [sequencing_preference, DEFAULT_PREF[0], DEFAULT_PREF[1]]
        
        # organize capture kits by abundance for purpose of breaking ties
        # for duplicate CHGVIDs
        kit_counter = {}
        for kit in kits:
            kit_counter[kit] = sum(d_cap[kit][1:])
        kit_counter = OrderedDict(sorted(kit_counter.items(), key=lambda x: x[1],
                                         reverse=True))
        
        # create new dict with CHGVIDs removed
        d_samp_clean = {}
        dups = 0
        for k, v in d_samp.iteritems():
            if len(v) == 1: # if only one record for this CHGVID, choose that one
                d_samp_clean[k] = v[0]
            else: # otherwise deal with the duplicate CHGVIDs
                d_samp_clean[k] = prune_samples(v, seq_prefs, kit_counter)
                dups += (len(v) - 1)
                if log_file_name:
                    logging.info('Removed {x} duplicate sample(s) for {sample}.'.format(
                        x=len(v) - 1, sample=k))
        
        # warn user of pruned duplicate samples
        msg = '\n{num} samples removed to avoid duplicates.\n'.format(num=dups)
        print(msg)
        if log_file_name:
            logging.info(msg)
        
        # create sample file
        with open(sample_file_name, 'w') as outfile:
            for k, v in d_samp_clean.iteritems():
                outfile.write('{fam}\t{sample}\t0\t0\t{sex}\t'
                              '{status}\t{seq}\t{kit}\n'.format(
                              fam=k,sample=k,sex=1 if v[2] == 'M' else 2, 
                              status=v[3], seq=v[0], kit=v[1]))
        # optional qc dump
        if qc_file_name:
            print('\nDumping QC data. This may take a few moments.\n')
            # pull down whole DB
            cur.execute("""
                        SELECT *
                        FROM SampleT as SAMP
                        INNER JOIN dragen_sample_metadata AS META
                        ON SAMP.CHGVID = META.sample_name AND
                        SAMP.SeqType = META.sample_type
                        INNER JOIN dragen_pipeline_step as STEP
                        ON META.pseudo_prepid = STEP.pseudo_prepid
                        INNER JOIN dragen_qc_metrics AS QC
                        ON META.pseudo_prepid = QC.pseudo_prepid
                        WHERE STEP.pipeline_step_id = 108 AND
                        STEP.step_status = 'completed' 
                        """)
            desc = cur.description # column names
            rows = cur.fetchall()
            # use seqtype, pseudoprepid as key to check for desired samples
            checklist = [(v[0].lower(),v[4]) for v in d_samp_clean.itervalues()]
            with open(qc_file_name, 'w') as outfile:
                # write column headers
                outfile.write('\t'.join([x[0] for x in desc])+'\n')
                # write all fields if the keys match
                for row in rows:
                    if ((row[4].lower(), row[38]) in checklist):
                        # write to outfile after removing pesky newlines in Notes
                        outfile.write('\t'.join(
                            [str(x).replace('\r','').replace('\n','') for x in row])+'\n')
    finally:
        if db.open:
            db.close()
     
     

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('sample_file_name', help='Provide name of output sample list')
    parser.add_argument('--sequencing_type', nargs='+', default=['Exome'],
                        choices=['Exome', 'Genome', 'Custom_Capture'],
                        help='Specify one or more types of sequencing to include')
    parser.add_argument('--sequencing_preference', default = 'Exome',
                        choices=['Exome', 'Genome', 'Custom_Capture'],
                        help='Specify a preferred sequencing type in case '
                        'a sample has multiple database entries')
    parser.add_argument('--include_sex_typing_mismatches', default=False,
                        action='store_true', help='Include samples where '
                        'self-declared sex does not match that inferred from '
                        'sequencing')
    parser.add_argument('--bases_covered_10x', default=0.9,
                        type=partial(valid_numerical_argument,
                                     arg_name='bases_covered_10x',
                                     min_value=0.0,
                                     max_value=1.0,
                                     arg_type=float),
                        help='Require all samples to have at least this fraction '
                        'of bases in the CCDS regions to have 10x coverage')
    parser.add_argument('--dbsnp_overlap_snv', default=0.95,
                        type=partial(valid_numerical_argument,
                                     arg_name='dbsnp_overlap_snv',
                                     min_value=0.0,
                                     max_value=1.0,
                                     arg_type=float),
                        help='Require this fraction of overlap with dbSNP (SNVs)')
    parser.add_argument('--max_contamination', default=0.03,
                        type=partial(valid_numerical_argument,
                                     arg_name= 'max_contamination',
                                     min_value=0.0, max_value=1.0,
                                     left_op=le, arg_type=float),
                        help='Require all samples to have <= this amount of '
                        'contamination')
    parser.add_argument('--ethnicity', nargs='+',
                        choices=['African', 'Caucasian', 'EastAsian', 
                                 'Hispanic', 'MiddleEastern', 'SouthAsian'],
                        help='Specify one or more permitted ethnicities')
    parser.add_argument('--min_genotyping_rate', default=0.9,
                        type=partial(valid_numerical_argument,
                                     arg_name='genotyping_rate',
                                     min_value=0.0, max_value=1.0, arg_type=float),
                        help='Require all samples to have at least this '
                        'genotyping rate for ethnicity prediction')
    parser.add_argument('--log_file_name', default=None,
                        help='Provide name of log file name, if desired')
    parser.add_argument('--qc_file_name', help='Provide name of output sample QC file')
    args = parser.parse_args()
    build_cohort(args.sample_file_name, args.sequencing_type, args.sequencing_preference,
                 args.include_sex_typing_mismatches, args.bases_covered_10x,
                 args.dbsnp_overlap_snv,
                 args.max_contamination, args.ethnicity, args.min_genotyping_rate,
                 args.log_file_name, args.qc_file_name)

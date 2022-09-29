#! /usr/bin/python3
import sys
import pandas as pd
import numpy as np
import re
import multiprocessing as mp
import os
from os.path import join as pjoin
import argparse
import argcomplete
import warnings
import datetime


class ArgumentError(Exception):
    pass


warnings.simplefilter(action='ignore')


def to_select(sample: str, mode: str) -> list:
    cl2select =["SAMPLE",
                "CHROM",
                "POS",
                "REF",
                "ALT",
                "HET",
                "%s.GQ" % sample,
                "%s.AD" % sample,
                "%s.DP" % sample,
                "%s.GT" % sample,
                "MULTIALLELIC",
                "HOM-REF",
                "HET.1",
                "HOM-VAR",
                "AC",
                "AN",
                "AF",
                "HOM-REF-SAMPLES",
                "HET-SAMPLES",
                "HOM-VAR-SAMPLES",
                "TOTAL-SAMPLES",
                "ANNOVAR_refseq_Effect",
                "ANNOVAR_refseq_Transcript_ID",
                "ANNOVAR_refseq_Gene_ID",
                "ANNOVAR_refseq_Closest_gene(intergenic_only)",
                "ANNOVAR_refseq_HGVSc",
                "ANNOVAR_refseq_HGVSp",
                "ANNOVAR_refseq_Exon_Rank",
                "ANNOVAR_refseq_summary",
                "splicing_consensus_ada_score",
                "splicing_consensus_rf_score",
                "gnomAD_exomes_flag",
                "gnomAD_exomes_AC",
                "gnomAD_exomes_AN",
                "gnomAD_exomes_AF",
                "gnomAD_exomes_AFR_AC",
                "gnomAD_exomes_AFR_AN",
                "gnomAD_exomes_AFR_AF",
                "gnomAD_exomes_AMR_AC",
                "gnomAD_exomes_AMR_AN",
                "gnomAD_exomes_AMR_AF",
                "gnomAD_exomes_ASJ_AC",
                "gnomAD_exomes_ASJ_AN",
                "gnomAD_exomes_ASJ_AF",
                "gnomAD_exomes_EAS_AC",
                "gnomAD_exomes_EAS_AN",
                "gnomAD_exomes_EAS_AF",
                "gnomAD_exomes_FIN_AC",
                "gnomAD_exomes_FIN_AN",
                "gnomAD_exomes_FIN_AF",
                "gnomAD_exomes_NFE_AC",
                "gnomAD_exomes_NFE_AN",
                "gnomAD_exomes_NFE_AF",
                "gnomAD_exomes_SAS_AC",
                "gnomAD_exomes_SAS_AN",
                "gnomAD_exomes_SAS_AF",
                "gnomAD_exomes_OTH_AC",
                "gnomAD_exomes_OTH_AN",
                "gnomAD_exomes_OTH_AF",
                "CADD_phred",
                "SIFT_score",
                "Polyphen2_HDIV_score",
                "Polyphen2_HVAR_score",
                "MutationTaster_pred",
                "MutationAssessor_pred"]
    if mode == 'INDEL':
        return cl2select[:-6]
    elif mode == 'SNP':
        return cl2select
    else:
        raise ValueError('Name of file %s is not in {sample_name}.{INDEL | SNP}.annotated.xlsx format')


def parseopts():
    parser = argparse.ArgumentParser(
        description="""Indicate a single file or a 
        directory containing multiple files.""",
        prefix_chars="--")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-s", "--single", type=str,
                        help="Path to a single INDEL or SNP file.",
                        action="store", required=False)
    group.add_argument("-d", "--directory", type=str,
                        help="Path to directory to elaborate",
                        action="store", required=False)
    parser.add_argument("-p", "--parallel", type=int,
                        help="Multiprocessing: Parallel degree number (default 2)\nunused in single mode",
                        action="store", required=False, default=2)
    parser.add_argument("-a", "--aux1", type=str,
                        help="""Auxiliary data sheet containing SAMPLE as sample
                                name and PHENOTYPE as diagnosis or ethnic labels column.
                                It is a CSV file format.
                                """,
                        action="store", required=("--directory" in sys.argv) or
                                                 ("-d" in sys.argv))
    parser.add_argument("-b", "--aux2", type=str,
                        help="Auxiliary directory containing OMIM genes classification (e.g. D/r).",
                        action="store", required=("--directory" in sys.argv) or
                                                 ("-d" in sys.argv))
    parser.add_argument("-o", "--output", type=str,
                        help="Path to output directory. (default ./)",
                        action="store", required=False, default='./')
    parser.add_argument('-v', '--verification',
                        help="It verifies whether each sample in sample sheet has files (SNP, INDEL) in input directory",
                        action='store_true', required=False, default=False)
    argcomplete.autocomplete(parser)
    return parser


def CLA(g) -> bool:
    g = list(g)
    y = g[-1]
    x = g[0]
    a, b = x.split('-')
    a = int(a)
    b = int(b)
    if y == 'HET':
        if a == 0:
            a = 1

        if b/a >= 0.25:
            return True
        else:
            return False
    elif y == 'HOM':
        if b != 0:
            if a/b <= 0.05:
                return True
            else:
                return False
        else:
            return False
    else:
        raise NotImplementedError('il caso in %s non è implementato' % str(g))


def c234_transform(slot: str) -> np.float:
    try:
        spt = re.split('\||;', slot)
    except TypeError:
        return float(slot)
    score: np.float = 0
    fold: np.float = 0
    for el in spt:
        if el != '.':
            score += np.float32(el)
            fold += 1
    if score == 0 and fold == 0:
        return np.NaN
    else:
        return score/fold


def MutationTaster_transf(slot: str) -> np.float:
    spt = re.split("\||;", slot)
    l = list(pd.Series(spt).unique())
    if len(l) == 1 and l[0] == '.':
        return np.NaN
    elif 'A' in l or 'D' in l:
        return 1
    else:
        return 0


def MutationAssessor_transf(slot: str) -> np.float:
    grading = {'H': 1, 'M': 0.5, 'L': 0.25, 'N': 0}
    spt = re.split("\||;", slot)
    l = [grading[el] for el in list(pd.Series(spt).unique()) if el != '.']
    if len(l) == 0:
        return np.NaN
    else:
        return min(l)


def pathogenic_index(table: pd.DataFrame):
    table = table.astype({'CADD_phred': np.int32})
    aux = pd.DataFrame()
    aux['EN'] = table['CADD_phred']
    aux['EN'] = np.NaN

    aux.loc[table.CADD_phred >= 20, 'EN'] = 1
    aux.loc[table.CADD_phred < 20, 'EN'] = 0

    aux['FEm'] = table.SIFT_score.apply(c234_transform)
    aux['FE'] = np.NaN
    aux.loc[aux['FEm'] <= 0.05, 'FE'] = 1
    aux.loc[aux['FEm'] > 0.05, 'FE'] = 0
    aux.drop('FEm', axis=1, inplace=True)

    aux['FJm'] = table.Polyphen2_HDIV_score.apply(c234_transform)
    aux['FJ'] = np.NaN
    aux.loc[aux['FJm'] > 0.908, 'FJ'] = 1
    aux.loc[((aux['FJm'] <= 0.908) & (aux['FJm'] >= 0.446)), 'FJ'] = 0.5
    aux.loc[aux['FJm'] < 0.446, 'FJ'] = 0
    aux.drop('FJm', axis=1, inplace=True)

    aux['FLm'] = table.Polyphen2_HVAR_score.apply(c234_transform)
    aux['FL'] = np.NaN
    aux.loc[aux['FLm'] > 0.908, 'FL'] = 1
    aux.loc[((aux['FLm'] <= 0.908) & (aux['FLm'] >= 0.446)), 'FL'] = 0.5
    aux.loc[aux['FLm'] < 0.446, 'FL'] = 0
    aux.drop('FLm', axis=1, inplace=True)

    aux['FR'] = table.MutationTaster_pred.apply(MutationTaster_transf)

    aux['FT'] = table.MutationAssessor_pred.apply(MutationAssessor_transf)

    aux['PI'] = aux.mean(axis=1)
    return aux.mean(axis=1)


def return_time(message: str = '') -> str:
    return "[%s] %s" % (str(datetime.datetime.now()), message)


def multifilter(pathfile: str, outpath: str, phen: str = '_') -> list:
    # It performs all filters except AF
    log = list()
    if pathfile.endswith('.xlsx') or pathfile.endswith('.xls'):
        log.append(return_time("%s loading..." % pathfile.split('/')[-1]))
        ID = pd.read_excel(pathfile)
    elif pathfile.endswith('.csv') or pathfile.endswith('.csv.gz') or pathfile.endswith('.csv.gzip'):
        log.append(return_time("%s loading..." % pathfile.split('/')[-1]))
        ID = pd.read_csv(pathfile)
    else:
        raise NotImplementedError('The %s file format is not implemented in this tool.' %
                                  pathfile.split('/')[-1])
    if '/' in pathfile:
        f_name = pathfile.split('/')[-1]
    else:
        f_name = pathfile
    sample_name = f_name.split('.')[0]
    mode = f_name.split('.')[1]
    ID = ID[to_select(sample_name, mode)]
    ID = ID[(~ID.SAMPLE.isna()) & (~ID.CHROM.isna()) & (~ID.POS.isna())]
    before = ID.shape[0]
    try:
        DP = ID[ID['%s.DP' % sample_name] >= 10]
        after = DP.shape[0]
        log.append(return_time('[%s %s] The number of variants before filter DP >= 10 is %s. After is %s' %
                               (sample_name, mode, before, after)))
        del ID
    except TypeError as err:
        print(sample_name, err)

    DP['gnomAD_exomes_AF'].replace('.', np.NaN, inplace=True)
    DP['gnomAD_exomes_AFR_AF'].replace('.', np.NaN, inplace=True)
    DP['gnomAD_exomes_AMR_AF'].replace('.', np.NaN, inplace=True)
    DP['gnomAD_exomes_ASJ_AF'].replace('.', np.NaN, inplace=True)
    DP['gnomAD_exomes_EAS_AF'].replace('.', np.NaN, inplace=True)
    DP['gnomAD_exomes_FIN_AF'].replace('.', np.NaN, inplace=True)
    DP['gnomAD_exomes_NFE_AF'].replace('.', np.NaN, inplace=True)
    DP['gnomAD_exomes_OTH_AF'].replace('.', np.NaN, inplace=True)
    DP['gnomAD_exomes_SAS_AF'].replace('.', np.NaN, inplace=True)
    DP = DP.astype({'gnomAD_exomes_AF': np.float,
                                'gnomAD_exomes_AFR_AF': np.float,
                                'gnomAD_exomes_AMR_AF': np.float,
                                'gnomAD_exomes_ASJ_AF': np.float,
                                'gnomAD_exomes_EAS_AF': np.float,
                                'gnomAD_exomes_FIN_AF': np.float,
                                'gnomAD_exomes_NFE_AF': np.float,
                                'gnomAD_exomes_OTH_AF': np.float,
                                'gnomAD_exomes_SAS_AF': np.float},)
    before = after
    MAF = DP[(DP['gnomAD_exomes_AF'] <= 0.01) | (DP['gnomAD_exomes_AF'].isna())]
    after = MAF.shape[0]
    log.append(return_time('[%s %s] The number of variants before filter MAF <= 0.01 or n.a. is %s. After is %s' %
          (sample_name, mode, before, after)))
    del DP
    before = after
    MAF = MAF[((MAF['gnomAD_exomes_AFR_AF'] <= 0.01) | (MAF['gnomAD_exomes_AFR_AF'].isna())) &
              ((MAF['gnomAD_exomes_AMR_AF'] <= 0.01) | (MAF['gnomAD_exomes_AMR_AF'].isna())) &
              ((MAF['gnomAD_exomes_ASJ_AF'] <= 0.01) | (MAF['gnomAD_exomes_ASJ_AF'].isna())) &
              ((MAF['gnomAD_exomes_EAS_AF'] <= 0.01) | (MAF['gnomAD_exomes_EAS_AF'].isna())) &
              ((MAF['gnomAD_exomes_FIN_AF'] <= 0.01) | (MAF['gnomAD_exomes_FIN_AF'].isna())) &
              ((MAF['gnomAD_exomes_NFE_AF'] <= 0.01) | (MAF['gnomAD_exomes_NFE_AF'].isna())) &
              ((MAF['gnomAD_exomes_OTH_AF'] <= 0.01) | (MAF['gnomAD_exomes_OTH_AF'].isna())) &
              ((MAF['gnomAD_exomes_SAS_AF'] <= 0.01) | (MAF['gnomAD_exomes_SAS_AF'].isna()))]
    after = MAF.shape[0]
    log.append(return_time('[%s %s] The number of variants before filter MAF <= 0.01 or n.a on the other populations is %s. After is %s' %
          (sample_name, mode, before, after)))
    NR = MAF['ANNOVAR_refseq_Effect'].apply(lambda x: 'ncRNA' not in x)
    EX = MAF['ANNOVAR_refseq_Effect'].apply(lambda x: 'exonic' in x)
    SP = MAF['ANNOVAR_refseq_Effect'].apply(lambda x: 'splicing' in x)
    if mode == 'INDEL':
        FS_NFS = MAF['ANNOVAR_refseq_Effect'].apply(lambda x: 'frameshift' in x)
        before = after
        Coding = MAF[(FS_NFS | EX | SP) & NR]
        after = Coding.shape[0]
        log.append(return_time('[%s %s] The number of variants before filter coding is %s. After is %s' %
              (sample_name, mode, before, after)))
    elif mode == 'SNP':
        NS = MAF['ANNOVAR_refseq_Effect'].apply(lambda x: 'nonsynonymous' in x)
        ST = MAF['ANNOVAR_refseq_Effect'].apply(lambda x: 'stop' in x)
        before = after
        Coding = MAF[(EX | SP | NS | ST) & NR]
        after = Coding.shape[0]
        log.append(return_time('[%s %s] The number of variants before filter coding is %s. After is %s' %
              (sample_name, mode, before, after)))
    else:
        raise ValueError('Name of file %s is not in {sample_name}.{INDEL | SNP}.annotated.xlsx format')

    before = after
    Ratio = Coding[Coding[['%s.AD' % sample_name, 'HET']].apply(lambda g: CLA(g), axis=1)]
    after = Ratio.shape[0]
    log.append(return_time('[%s %s] The number of variants before HOM/HET filter is %s. After is %s' %
          (sample_name, mode, before, after)))
    if phen != '_':
        store = pjoin(outpath, phen, mode, '%s.%s.filtered.csv' % (sample_name, mode))
    else:
        store = pjoin(outpath,'%s.%s.filtered.csv' % (sample_name, mode))

    if mode == 'SNP':
        Ratio['PIK'] = pathogenic_index(Ratio)
        before = Ratio.shape[0]
        Ratio = Ratio[Ratio['PIK'] >= 0.7]
        after = Ratio.shape[0]
        log.append(return_time('[%s %s] The number of variants before IPK filter is %s. After is %s' %
              (sample_name, mode, before, after)))

        Ratio.to_csv(store)
        log.append(return_time('[%s %s] Filtered sample stored at %s.' % (sample_name, mode, store)))
    elif mode == 'INDEL':
        Ratio.to_csv(store)
        log.append(return_time('[%s %s] Filtered sample stored at %s.' % (sample_name, mode, store)))

    return log


def read_columns(x):
    a = pd.read_csv(x, usecols=['CHROM', 'POS', 'REF', 'ALT', 'HET'])
    # print(a)
    return a


def hom_het_transform(statement: str) -> int:
    if statement == 'HOM':
        return 2
    elif statement == 'HET':
        return 1
    elif statement == 'pC-HET':
        return 1
    else:
        raise NotImplementedError("The statement management for %s is not implemented" % statement)


def countAF(to_import: list, processes: int) -> pd.DataFrame:
    # to_import = [os.path.join(directory, el) for el in os.listdir(directory)]
    with mp.Pool(processes=processes, maxtasksperchild=1) as pl:
        cbk = pl.map(func=read_columns, iterable=to_import,
                     chunksize=len(to_import) // processes)

    # print(cbk)
    max_df = pd.concat(cbk, axis=0)
    max_df['count'] = max_df['HET'].apply(lambda x: hom_het_transform(x))
    count = pd.DataFrame(max_df.groupby(['CHROM', 'POS', 'REF', 'ALT'])['count'].sum())

    count.columns = ['count']
    count.reset_index(drop=False, inplace=True)
    count['AF'] = count['count'] / (data_sheet.shape[0] * 2)
    return count


def filterAF(outdir: str, cond: str, mode: str, fl: str) -> pd.DataFrame:
    new = pd.read_csv(pjoin(outdir, cond, mode, fl), index_col=0)
    before = new.shape[0]
    new.drop('AF', axis=1, inplace=True)
    new = new.merge(af_recalc, on=['CHROM', 'POS', 'REF', 'ALT'])
    new.columns = [el.replace(fl.split('.')[0], 'SAMPLE') for el in new.columns]
    # print(new.columns)
    new = new[new['AF'] <= 0.01]
    after = new.shape[0]
    new.to_csv(pjoin(outdir, cond, mode, fl))
    print('[%s %s] Before AF <= 0.01 is %s. After is %s' %
                (fl.split('_')[0], mode, before, after))
    return new


def samples_verif(samples: pd.Series, directory: str):
    print('''
          
                 *** SAMPLES CHECKER MODULE EXECUTION ***
          
          ''')
    dict_f = dict()
    for f in sorted(os.listdir(directory)):
        sample = f.split('.')[0]
        dict_f.setdefault(sample, list())
        if 'SNP' in f:
            dict_f[sample].append(0)
        elif 'INDEL' in f:
            dict_f[sample].append(1)
    c = 0
    for sample in samples:
        if sample in dict_f:
            if (0 in dict_f[sample]) and (1 in dict_f[sample]):
                pass
            elif (0 not in dict_f[sample]) and (1 in dict_f[sample]):
                c += 1
                print('WARNING: SNP file for %s sample is not in %s directory' % (sample, directory))
            elif (0 in dict_f[sample]) and (1 not in dict_f[sample]):
                c += 1
                print('WARNING: INDEL file for %s sample is not in %s directory' % (sample, directory))
            elif (0 not in dict_f[sample]) and (1 not in dict_f[sample]):
                c += 1
                print('WARNING: SNP and INDEL files for %s sample are not in %s directory or file names do not match standards.'
                      % (sample, directory))
            else:
                raise NotImplementedError('for sample %s' % sample)
        else:
            c += 1
            print('WARNING: SNP and INDEL files for %s sample are not in %s directory' % (sample, directory))
    if c == 0:
        print('Samples check is ok. You can launch the analysis without -v option\n')
    else:
        print('Samples check is not ok. You should manually check samples and files in directory\n')


if __name__ == '__main__':
    print(datetime.datetime.now())
    parser_obj = parseopts()
    args = parser_obj.parse_args()
    single, directory, processes, aux1, aux2, outdir, verif = \
        args.single, args.directory, args.parallel, args.aux1,args.aux2, \
        args.output, args.verification

    print("""
    ####################################################################
    # SNPs & INDELs filtering tool, resuming Krausz filtering pipeline #
    #             code_developers: Defazio,G; Farnetani,G              #
    ####################################################################
    
    """)

    if directory is not None and os.path.exists(directory) and os.path.isdir(directory):
        data_sheet = pd.read_csv(aux1)
        if verif:
            samples_verif(data_sheet.SAMPLE, directory)
        else:
            print("""
            Input directory: %s
            Processes: %s
            Output directory: %s
            _________________
    
            """ % (directory, processes, outdir))
            conds = data_sheet.PHENOTYPE.unique()
            #silenziato solo per continuare la pipeline
            try:
                os.mkdir(outdir)
                for cond in conds:
                    os.mkdir(pjoin(outdir, cond))
                    os.mkdir(pjoin(outdir, cond, 'SNP'))
                    os.mkdir(pjoin(outdir, cond, 'INDEL'))
            except FileExistsError:
                print("The directory %s already exists and contains %s objects" %
                      (outdir, os.listdir(outdir).__len__()))
            # print(int(len(os.listdir(directory))/processes))
            iterable = list()
            for el in data_sheet.SAMPLE:
                for f in os.listdir(directory):
                    if f.split('.')[0] == el:
                        # print(el, f, data_sheet[data_sheet.SAMPLE == el.split('.')[0]]['PHENOTYPE'].iloc[0])
                        iterable.append((pjoin(directory, f), outdir, data_sheet[data_sheet.SAMPLE == el.split('.')[0]]['PHENOTYPE'].iloc[0]))
            try:
                pl = mp.Pool(processes=processes, maxtasksperchild=1)
                cbk = pl.starmap(func=multifilter,
                                 iterable=iterable,
                                 chunksize=len(iterable)//processes)
                pl.close()
                pl.join()
            except IndexError:
                for el in os.listdir(directory):
                    try:
                        print(data_sheet[data_sheet.SAMPLE == el.split('.')[0]]['PHENOTYPE'].iloc[0])
                    except IndexError:
                        print(el)
                raise Exception
            for el in cbk:
                print('\n'.join(el))
            return_time('AF recalculation START')
            to_conc4af = list()
            for cond in conds:
                for mode in ['SNP', 'INDEL']:
                    to_concat = list()
                    for fl in os.listdir(pjoin(outdir, cond, mode)):
                        to_conc4af.append(pjoin(outdir, cond, mode, fl))
            af_recalc = countAF(to_conc4af, processes)
            af_recalc.to_csv(os.path.join(outdir, 'AF_recalculation.csv'))
            return_time('AF recalculation STOP & AF filtering START')
            for cond in conds:
                for mode in ['SNP', 'INDEL']:
                    # adornd_filterAF = lambda x: filterAF(outdir, cond, mode, x)
                    with mp.Pool(processes=processes,maxtasksperchild=1) as afpl:
                        pts = [[outdir, cond, mode, x] for x in os.listdir(pjoin(outdir, cond, mode))]
                        to_concat = afpl.starmap(func=filterAF,
                                                 iterable=pts,
                                                 chunksize= len(pts)//processes)
                    afpl.close()
                    afpl.join()
                    resume_var = pd.concat(to_concat, axis=0)
                    resume_var.to_csv(pjoin(outdir, "%s_%s_filtered.csv" % (cond, mode)), index=False)
            return_time('AF filtering STOP')
            for cond in conds:
                ctn = pd.concat([pd.read_csv(pjoin(outdir, "%s_%s_filtered.csv" % (cond, 'SNP'))),
                                 pd.read_csv(pjoin(outdir, "%s_%s_filtered.csv" % (cond, 'INDEL')))],
                                axis=0)
                ctn.reset_index(drop=True, inplace=True)

                ctn['gene'] = ctn['ANNOVAR_refseq_Gene_ID'].apply(lambda x: x.split('|')[0])
                comphet = list()
                clmns = ['SAMPLE', 'CHROM', 'gene', 'POS', 'REF', 'ALT', 'HET']
                for el in ctn['SAMPLE'].unique():
                    sub = ctn[ctn['SAMPLE'] == el]
                    ct = sub.gene.value_counts()
                    genes = list(ct[ct > 1].index)
                    comphet.append(sub[sub.gene.isin(genes)][clmns].sort_values(by='gene'))
                comp_het = pd.concat(comphet, axis=0)
                comp_het['HET'] = 'pC-HET'
                ctn.loc[comp_het.index, 'HET'] = 'pC-HET'

                hom = ctn[ctn['HET'] == 'HOM'][clmns]

                X = ctn[ctn['CHROM'] == 'chrX'][clmns]
                pre_recessive = pd.concat([comp_het, hom, X], axis=0).drop_duplicates(ignore_index=True)

                path_omim_recessive = pjoin(aux2, 'OMIM_recessive_genes_approved_symbol.pkl')
                omim_recessive = pd.read_pickle(path_omim_recessive)
                recessive = pre_recessive[pre_recessive['gene'].isin(omim_recessive.symbol)]

                recessive.to_csv(pjoin(outdir, "%s_RECESSIVE.csv" % cond), index=False)

                dominant = ctn[ctn['HET'] == 'HET'][clmns]
                path_omim_dominant = pjoin(aux2, 'OMIM_dominant_genes_approved_symbol.pkl')
                omim_dominant = pd.read_pickle(path_omim_dominant)
                dominant = dominant[dominant['gene'].isin(omim_dominant.symbol)]
                dominant = dominant[~dominant.set_index(['SAMPLE', 'gene']).index.isin(
                    pre_recessive.set_index(['SAMPLE', 'gene']).index)]

                dominant.to_csv(pjoin(outdir, "%s_DOMINANT.csv" % cond), index=False)
                ctn = ctn[clmns]
                remaining = ctn[(~ctn.set_index(['SAMPLE','gene','POS']).index.isin(dominant.set_index(
                    ['SAMPLE','gene','POS']).index)) & (~ctn.set_index(['SAMPLE','gene','POS']).index.isin(
                    recessive.set_index(['SAMPLE','gene','POS']).index))]
                remaining.to_csv(pjoin(outdir,'%s_noDOM_noREC.csv' % cond))
                genes = pd.Series(recessive.gene.to_list() + dominant.gene.to_list()).unique()
                with open(pjoin(outdir,'./%s_genes4GO.txt') % cond, 'wt') as gg:
                    for g in genes:
                        gg.write('%s\n' % g)

    elif os.path.exists(single) and os.path.isdir(single) is False:
        multifilter(single, outdir)
    print(datetime.datetime.now())

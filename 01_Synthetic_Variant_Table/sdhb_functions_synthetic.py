from Bio import SeqIO
import pandas as pd
import numpy as np
import re
from gtfparse import read_gtf
import math
import os
import pickle
import matplotlib.pyplot as plt
import matplotlib.patches as pplt
from functools import reduce

# search for specific column 
def search_col(inp, df):
    return [i for i in df.columns if inp.lower() in i.lower()]


# create synthetic variants
def synth_vars(sdhb_fa_path = '00_Andres_SDHB/00_sdhb_seq_dfs/sdhb_fasta',
               sdhb_syn_out = '00_Andres_SDHB/00_sdhb_seq_dfs/sdhb_syn.csv'):
    # open new document, referred to as 'o'
    # either every transcript or only reference transcript
    with open(sdhb_syn_out,'w') as o:
        # result in dictionary format
        resultDict={}
        # write first line in new document 'o'
        #o.write('Identifier,CHROM,POS,REF,ALT,Ref_Genome,Strand,Ref_Transcript,Trans_Version_Syn\n')
        o.write('Identifier,CHROM,POS,REF,ALT,Strand,Transcript\n')
        # open document with FASTA sequence, referred to as 'f'
        # FASTA sequence by UCSC Genome Browser
        with open(sdhb_fa_path) as f:
            # parse: all single features as a list iterator
            for record in SeqIO.parse(f, 'fasta'):
                # .split()[].split(): split by (), then take the piese [] and split it by ()
                # take first number from position range as variable 's', with iteration through sequence s+=1 = POS
                s=int(record.description.split('range')[1].split(':')[1].split('-')[0])
                # take reference genome

                #refg = record.description.split('_')[0]
                # take chromosome number
                chrom=record.description.split('range')[1].split('chr')[1].split(':')[0]
                # take strand +/-
                strand=record.description.split('strand=')[1].split(' repeat')[0]
                # take transcript/record of fasta file
                trans=record.id
                # for transcript version
                #trans_vers = record.description.split('.')[1].split('_')[0]
                # for '-' strands: reverse complement and reverse exon order + 1
                if strand == '-':
                    record.seq = record.seq.reverse_complement()
                else:
                    pass                
                # only for reference transcript
                #if trans == refseqtr:
                # iterate through sequence by nucleotide
                for nuc in str(record.seq):
                    # all nucleotides in upper case
                    nuc = str(nuc.upper())
                    # for specific nucleotide, print list of alternative nucleotides
                    alt = 'ACGT'.replace(str(nuc), '')
                    # for nuc in alternative nucleotides entry for dictionary
                    for cr in alt:
                        # entr = specific position of each nucleotide (for key in dictionary)
                        entr = (f'{chrom}_{str(s)}_{nuc}_{cr}_{strand}_{trans}')
                        # if nucleotide with specific position is not yet in resulting dictionary
                        if entr not in resultDict:
                            # append specific transcript as value of the specific position (= key) in dictionary
                            resultDict[entr]=trans
                        else:
                            # if specific nucleotide position is already in dictionary, append other transcript(s)
                            resultDict[entr]=resultDict[entr]#+' '+trans
                    s+=1

        # address dictionary with both features (key and value), sort the dictionary
        #i = 0
        for key,value in sorted(resultDict.items()):
            # every key element (specific position element) is connected by '_' --> split by '_' to access each key
            x=key.split('_')
            # for identifier to integrate information to specific position later on
            identi = str(x[0]) + ':g.' + str(x[1]) + str(x[2]) + '>' + str(x[3])
            # write in new document 'o'
            o.write(f'{str(identi)},{str(x[0])},{str(x[1])},{str(x[2])},{str(x[3])},{str(x[4])},{value}\n')


# variant in exonic or intronic region
def exon_or_intron(sdhb_df, gtf_sdbh):
    sdhb_df['exon_intron'] = np.nan
    # exonic regions
    for start,stop,exon in zip(gtf_sdbh[gtf_sdbh['feature']=='exon']['start'], 
                            gtf_sdbh[gtf_sdbh['feature']=='exon']['end'], 
                            gtf_sdbh[gtf_sdbh['feature']=='exon']['exon_number']):
        sdhb_df.loc[((sdhb_df['POS']>=start)&(sdhb_df['POS']<=stop)),'exon_intron'] = f'exon_{exon}'

    # intronic regions
    intr_start = []
    intr_end = []
    intr_name = []
    for i,e_end, e_start in zip(list(gtf_sdbh[gtf_sdbh['feature']=='exon']['exon_number'][:-1]),
                                list(gtf_sdbh[gtf_sdbh['feature']=='exon']['end'][1:]), 
                                list(gtf_sdbh[gtf_sdbh['feature']=='exon']['start'][:-1])):
        sdhb_df.loc[((sdhb_df['POS']>=(e_end+1))&(sdhb_df['POS']<=(e_start-1))),'exon_intron'] = f'intron_{i}'
        intr_start.append(e_end+1)
        intr_end.append(e_start-1)
        intr_name.append(f'intron_{i}')

    intron_gtf = pd.DataFrame()
    for col,n in zip([intr_start, intr_end, intr_name], ['start','end','name']):
        intron_gtf[n] = col
    intron_gtf['len_intron'] = intron_gtf['end']-intron_gtf['start']

    # upstream and downstream 1000 bases
    gene_start = gtf_sdbh[gtf_sdbh['feature']=='gene'].iloc[0]['start']
    gene_stop = gtf_sdbh[gtf_sdbh['feature']=='gene'].iloc[0]['end']
    sdhb_df.loc[(sdhb_df['POS']>gene_stop),'exon_intron'] = 'promoter_upstream'
    sdhb_df.loc[(sdhb_df['POS']<gene_start),'exon_intron'] = 'downstream'

    # merge exonic and intronic gtf
    intr_ex_gtf_sdhb = pd.merge(gtf_sdbh[gtf_sdbh['feature']=='exon'], intron_gtf, how='outer')
    intr_ex_gtf_sdhb['intr_ex'] = np.nan
    intr_ex_gtf_sdhb.loc[intr_ex_gtf_sdhb['exon_number'].notna(),'intr_ex'] = intr_ex_gtf_sdhb['exon_number']
    intr_ex_gtf_sdhb.loc[intr_ex_gtf_sdhb['name'].notna(),'intr_ex'] = intr_ex_gtf_sdhb['name'].str.strip('intron_')
    intr_ex_gtf_sdhb.loc[intr_ex_gtf_sdhb['feature'].isna(),'feature'] = 'intron'
    intr_ex_gtf_sdhb.loc[intr_ex_gtf_sdhb['seqname'].isna(),'seqname'] = intr_ex_gtf_sdhb['seqname'][
        intr_ex_gtf_sdhb['exon_number'].notna()][0]
    intr_ex_gtf_sdhb.loc[intr_ex_gtf_sdhb['strand'].isna(),'strand'] = intr_ex_gtf_sdhb['strand'][
        intr_ex_gtf_sdhb['strand'].notna()][0]
    intr_ex_gtf_sdhb = intr_ex_gtf_sdhb[['feature','seqname','start','end','strand','intr_ex']]

    return intr_ex_gtf_sdhb


# include certain region into intron
# num = up to which position into intron, last_num = last highest num (to not have unnecessary duplicates)
def intron_df_vep(num, last_num):
    if last_num > num:
        print('Error: last_num > num')
    else:
        intron_df = pd.DataFrame()
        for i_start, i_end, i_len in zip(intron_gtf['start'], intron_gtf['end'], intron_gtf['len_intron']):
            half_intr_len = math.ceil(i_len/2)
            if num >= half_intr_len:
                use_num = half_intr_len
            else:
                use_num = num
            intron_df = intron_df.append(sdhb_df[(sdhb_df['POS']>=i_start+last_num)&
                                                 (sdhb_df['POS']<i_start+use_num)])
            intron_df = intron_df.append(sdhb_df[(sdhb_df['POS']>i_end-use_num)&
                                                 (sdhb_df['POS']<=i_end-last_num)])
        intron_df = intron_df.drop_duplicates().reset_index(drop=True)
        intron_df['Identifier'].to_csv(f'00_Andres_SDHB/01_syn_tables_out/synth_sdhb_intron_{num}_variants.txt', 
                                       index=False, header=False)
        return intron_df


# adjust VEP output file
def vep_output(file_path):
    df = pd.read_table(file_path)
    df = df[df['Feature']=='NM_003000.3']
    
    for i in search_col('SpliceAI_pred_DS', df):
        df[i] = df[i].astype(float)
    for i in search_col('SpliceAI_pred_DP', df):
        df[i] = df[i].astype(int)
    
    df[['#CHROM','loc1']] = df['Location'].str.split(':', expand=True)
    df[['POS','pos2']] = df['loc1'].str.split('-', expand=True)

    df = df.rename(columns={'USED_REF':'REF', 'Allele':'ALT'})
    df['POS'] = df['POS'].astype(int)
    df['#CHROM'] = df['#CHROM'].astype(str)

    cols_df = ['#Uploaded_variation','#CHROM','POS','REF','ALT', 'Consequence', 'IMPACT', 'SYMBOL', 'Feature_type', 
               'Feature','EXON','INTRON', 'cDNA_position','Codons','STRAND','SIFT','PolyPhen','SpliceAI_pred_DP_AG',
               'SpliceAI_pred_DP_AL', 'SpliceAI_pred_DP_DG', 'SpliceAI_pred_DP_DL','SpliceAI_pred_DS_AG', 
               'SpliceAI_pred_DS_AL', 'SpliceAI_pred_DS_DG','SpliceAI_pred_DS_DL', 'SpliceAI_pred_SYMBOL', 
               'MaxEntScan_alt','MaxEntScan_diff', 'MaxEntScan_ref', 'CADD_PHRED', 'CADD_RAW','ada_score', 
               'rf_score']

    return df[cols_df].reset_index(drop=True)


# for all files in 02_vep_results -> run vep_output()
def run_all_vep(file_path = '00_Andres_SDHB/02_vep_results/'):
    vep_files_last = list(pd.read_table(file_path+'vep_dict_list')['files'])
    vep_files_now = [file for file in os.listdir(file_path) if '_sdhb_vep' in file]
    if set(vep_files_last) == set(vep_files_now):
        print('All results already in file "vep_dict.pickle"')
        with open(file_path+'vep_dict.pickle', 'rb') as o:
            vep_dfs = pickle.load(o)
    else:
        vep_dfs = {}
        for file in os.listdir(file_path):
            if file.endswith('.txt'):
                vep_dfs[file.strip('.txt')] = vep_output(file_path+file)
        with open(file_path+'vep_dict.pickle', 'wb') as o:
            pickle.dump(vep_dfs, o)
        pd.DataFrame(vep_files_now, columns=['files']).to_csv(file_path+'vep_dict_list', index=False)
    return vep_dfs


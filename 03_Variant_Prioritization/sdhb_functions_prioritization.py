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


### SpliceAI cutoff functions ###
def splai_cutoff(cutoff, df):
    cols_spl_rel = ['#CHROM','POS','REF','ALT','SpliceAI_pred_DP_AG','SpliceAI_pred_DP_AL', 'SpliceAI_pred_DP_DG', 
                    'SpliceAI_pred_DP_DL','SpliceAI_pred_DS_AG', 'SpliceAI_pred_DS_AL', 'SpliceAI_pred_DS_DG',
                    'SpliceAI_pred_DS_DL', 'SpliceAI_pred_SYMBOL', 'MaxEntScan_alt','MaxEntScan_diff', 
                    'MaxEntScan_ref', 'CADD_PHRED','EXON','INTRON']
    df_spl = df[(df['SpliceAI_pred_DS_AG']>cutoff)|(df['SpliceAI_pred_DS_AL']>cutoff)|
                (df['SpliceAI_pred_DS_DG']>cutoff)|(df['SpliceAI_pred_DS_DL']>cutoff)].copy()
    return df_spl[cols_spl_rel]

def splai_is0(df):
    cols_spl_rel = ['#CHROM','POS','REF','ALT','SpliceAI_pred_DP_AG','SpliceAI_pred_DP_AL', 'SpliceAI_pred_DP_DG', 
                    'SpliceAI_pred_DP_DL','SpliceAI_pred_DS_AG', 'SpliceAI_pred_DS_AL', 'SpliceAI_pred_DS_DG',
                    'SpliceAI_pred_DS_DL', 'SpliceAI_pred_SYMBOL', 'MaxEntScan_alt','MaxEntScan_diff', 
                    'MaxEntScan_ref', 'CADD_PHRED','EXON','INTRON']
    df_spl = df[(df['SpliceAI_pred_DS_AG']==0.0)&(df['SpliceAI_pred_DS_AL']==0.0)&
                (df['SpliceAI_pred_DS_DG']==0.0)&(df['SpliceAI_pred_DS_DL']==0.0)&(df['EXON']!='-')&
                (df['Consequence']=='synonymous_variant')].copy()
    return df_spl[cols_spl_rel]

def splai_spec_cutoff(cutoff, ds, df):
    cols_spl_rel = ['#CHROM','POS','REF','ALT','SpliceAI_pred_DP_AG','SpliceAI_pred_DP_AL', 'SpliceAI_pred_DP_DG', 
                    'SpliceAI_pred_DP_DL','SpliceAI_pred_DS_AG', 'SpliceAI_pred_DS_AL', 'SpliceAI_pred_DS_DG',
                    'SpliceAI_pred_DS_DL', 'SpliceAI_pred_SYMBOL', 'MaxEntScan_alt','MaxEntScan_diff', 
                    'MaxEntScan_ref', 'CADD_PHRED','EXON','INTRON']
    df_spl = df[df['SpliceAI_pred_DS_'+ds]>cutoff].copy()
    return df_spl[cols_spl_rel]

def splai_count_df(vep_dict):
    splai_df = pd.DataFrame()
    splai_df['name'] = vep_dict.keys()
    for i in np.arange(0,1,0.1):
        lst = []
        for name in vep_dict.keys():
            lst.append(len(splai_cutoff(i, vep_dict[name])))
        splai_df[round(i, ndigits=1)] = lst
    return splai_df.sort_values('name').reset_index(drop=True)

def splai_spec_count_df(vep_dict):
    splai_df = pd.DataFrame()
    splai_df['name'] = vep_dict.keys()
    for i in np.arange(0,1,0.1):
        for ds in ['AL','AG','DL','DG']:
            lst = []
            for name in vep_dict.keys():
                lst.append(len(splai_spec_cutoff(i, ds, vep_dict[name])))
            col_num = str(round(i, ndigits=1))
            splai_df[f'{ds}_{col_num}'] = lst
    return splai_df.sort_values('name').reset_index(drop=True)

# columns for interpreting SpliceAI results
def spl_columns(df,intr_ex_gtf_sdhb):
    spl_vars = df.copy()
    spl_vars['feature'] = np.nan
    spl_vars.loc[spl_vars['EXON']=='-','feature'] = 'intron'
    spl_vars.loc[spl_vars['INTRON']=='-','feature'] = 'exon'

    spl_vars['intr_ex'] = np.nan
    spl_vars[['intr_ex1','tr1']] = spl_vars['EXON'].str.split('/', expand=True)
    spl_vars.loc[~spl_vars['intr_ex1'].isin(['-',np.nan]),'intr_ex'] = spl_vars['intr_ex1']
    
    if len(spl_vars['INTRON'].unique()) > 1:
        spl_vars[['intr_ex2','tr2']] = spl_vars['INTRON'].str.split('/', expand=True)
        spl_vars.loc[~spl_vars['intr_ex2'].isin(['-',np.nan]),'intr_ex'] = spl_vars['intr_ex2']
        spl_vars = spl_vars.drop(columns=['intr_ex2','tr2'])
    
    spl_vars = spl_vars.drop(columns=['EXON','INTRON','intr_ex1','tr1'])

    acc_don = []
    dist_exin = []

    strand = intr_ex_gtf_sdhb['strand'][0]
    for pos, feat, num in zip(spl_vars['POS'],spl_vars['feature'],spl_vars['intr_ex']):
        start = int(intr_ex_gtf_sdhb['start'][(intr_ex_gtf_sdhb['feature']==feat)&(intr_ex_gtf_sdhb['intr_ex']==num)])
        end = int(intr_ex_gtf_sdhb['end'][(intr_ex_gtf_sdhb['feature']==feat)&(intr_ex_gtf_sdhb['intr_ex']==num)])
        pos = int(pos)

        if strand == '+':
            if feat == 'exon':
                acc_dist = pos - start
                don_dist = end - pos
            elif feat == 'intron':
                start = start-1
                end = end+1
                don_dist = pos - start
                acc_dist = end - pos
        else:
            if feat == 'exon':
                don_dist = pos - start
                acc_dist = end - pos
            elif feat == 'intron':
                start = start-1
                end = end+1
                acc_dist = pos - start
                don_dist = end - pos
        if acc_dist < don_dist:
            acc_don.append('acceptor')
            dist_exin.append(acc_dist)
        else:
            acc_don.append('donor')
            dist_exin.append(don_dist)

    spl_vars['acc_don'] = acc_don
    spl_vars['dist_exin'] = dist_exin
    
    return spl_vars


def make_splai_result_table(splai_out_vcf, used_vep_dicts, pickle_out):
    splai_new = pd.read_table(splai_out_vcf, comment='#', 
                              names=['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO'], 
                              dtype={'#CHROM':str, 'POS':int})

    col_list_info = ['ALLELE', 'SYMBOL', 'DS_AG', 'DS_AL', 'DS_DG', 'DS_DL', 'DP_AG', 'DP_AL', 'DP_DG', 'DP_DL',
                     'DS_AG_REF','DS_AG_ALT','DS_AL_REF','DS_AL_ALT','DS_DG_REF','DS_DG_ALT','DS_DL_REF','DS_DL_ALT']
    splai_new[col_list_info] = splai_new['INFO'].str.split('=', expand=True)[1].str.split('|', expand=True)
    splai_new = splai_new.drop(columns=['ID','QUAL','FILTER','INFO','ALLELE','SYMBOL'])
    
    vep_dict = pd.read_pickle('00_Andres_SDHB/02_vep_results/vep_dict.pickle')
    vep_dict_list = [vep_dict[d] for d in used_vep_dicts]

    vep = reduce(lambda left, right: pd.merge(left, right, how='outer'), vep_dict_list)

    vep_vars = spl_columns(vep)[['#CHROM', 'POS', 'REF', 'ALT', 'Consequence', 'feature', 'intr_ex', 'acc_don', 
                                 'dist_exin', 'MaxEntScan_alt', 'MaxEntScan_diff', 'MaxEntScan_ref', 'ada_score', 
                                 'rf_score', 'CADD_PHRED', 'CADD_RAW']].copy()

    annot_vars = pd.merge(splai_new, vep_vars, how='left', on=['#CHROM','POS','REF','ALT'])
    
    annot_vars[[i for i in col_list_info if 'DS' in i]] = annot_vars[[i for i in col_list_info if 'DS' in i]
                                                                    ].astype(float)
    annot_vars[[i for i in col_list_info if 'DP' in i]] = annot_vars[[i for i in col_list_info if 'DP' in i]
                                                                    ].astype(int)
    annot_vars[['dist_exin','intr_ex']] = annot_vars[['dist_exin','intr_ex']].astype(int)
    
    annot_vars.to_pickle(pickle_out)



### Visualization ###
# figure for exons and introns
def overview_fig(gene, all_vars, gtf, splice_vars):
    #ind = gene_names.index(str(gene))
    all_variants = all_vars.drop_duplicates('POS').copy().reset_index(drop=True).reset_index(drop=False)
    fig,ax = plt.subplots(figsize=(20, 7))
    plt.ylim([-1, 1])
    plt.axis('off')
    gene_start = 0
    real_start = '{:,.0f}'.format(min(all_variants['POS']))
    real_end = '{:,.0f}'.format(max(all_variants['POS']))
    gene_end = max(all_variants['index'])
    strand = gtf['strand'][0]
    chrom = gtf['seqname'][0]
    plt.xlim(0, gene_end)    
    fig.suptitle('Gene model of $%s$ gene in reading direction with SpliceAI scores' 
                 %(gene.upper()), fontsize=25, fontweight='bold')
    
    if strand=='+':
        ax.add_patch(pplt.Rectangle((gene_start, -.1), gene_end, .2, 
                                    fc = 'dodgerblue', ec = 'dodgerblue', alpha = 1, label = 'Introns'))
        ax.arrow(gene_start, -0.7, abs(gene_end-gene_start)-1000, 0, 
                 head_width = 0.05, head_length = 1000, fc = 'k', ec = 'k')
        for start,end,num in zip(gtf['start'][gtf['feature']=='exon'], gtf['end'][gtf['feature']=='exon'],
                                 gtf['intr_ex'][gtf['feature']=='exon']):
            ax.add_patch(pplt.Rectangle((start, -.3), abs(end-start), .6, fc = 'orange', ec = 'orange', alpha = 1))
            plt.text(int(start-200), -0.4, str(num), size = 22)
        #for tup,l in zip(tuple_exons[ind], range(len(tuple_exons[ind]))):
        #    ax.add_patch(pplt.Rectangle((int(tup[0]), -.3), abs(int(tup[0])-int(tup[1])), .6, fc = 'orange', 
        #                                ec = 'orange', alpha = 1))
        #    plt.text(int(tup[0]-300), txt_pos[ind][int(l)], str((int(l)+2)), size = 22)
        #plt.text(tuple_exons[ind][-1][1]-23000, -0.6, 'End: CHROM ' + str(dfs[ind].iloc[0]['CHROM']) + 
        #         ', POS ' + str('{:,.0f}'.format(tuple_exons[ind][-1][1])), size = 25)
        #plt.text((tuple_exons[ind][0][0] + abs(tuple_exons[ind][-1][1]-tuple_exons[ind][0][0])/2)-8000, -0.6, 
        #         'Gene Length: ' + str(abs(tuple_exons[ind][-1][1]-tuple_exons[ind][0][0])/1000) + ' kb', 
        #         size = 25)
    elif strand=='-':
        ax.add_patch(pplt.Rectangle((gene_start, -.1), gene_end, .2, 
                                    fc = 'dodgerblue', ec = 'dodgerblue', alpha = 1, label = 'Introns'))
        ax.arrow(gene_end, -0.7, -abs(gene_end-gene_start)+500, 0, 
                 head_width = 0.05, head_length = 500, fc = 'k', ec = 'k')
        for start,end,num in zip(gtf['start'][gtf['feature']=='exon'], gtf['end'][gtf['feature']=='exon'],
                                 gtf['intr_ex'][gtf['feature']=='exon']):
            pic_start = all_variants[all_variants['POS']==start].iloc[0]['index']
            pic_end = all_variants[all_variants['POS']==end].iloc[0]['index']
            ax.add_patch(pplt.Rectangle((pic_start, -.3), abs(pic_end-pic_start), .6, fc = 'orange', ec = 'orange', 
                                        alpha = 1))
            plt.text(pic_start+abs(pic_end-pic_start)/2, -0.4, str(num), size = 22, ha='center')
        for pos,pic_col,y_pos in zip(splice_vars['POS'], splice_vars['pic_col'], splice_vars['y_pos']):
            new_pos = all_variants[all_variants['POS']==pos].iloc[0]['index']
            plt.plot(new_pos, y_pos, marker='o', markersize=7, color=pic_col, alpha=0.5)
    
    for cutoff,col,x_pos in zip([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1], 
                                ['w','w','gray','k','b','g','c', 'y','m','orangered','firebrick'], 
                                [0,600,1200,1800,2400,3000,3600,4200,4800,5400,6000]):
        n_cutoff = '{:.1f}'.format(cutoff)
        #plus_cutoff = float(n_cutoff)+0.1
        #plus_cutoff = '{:.1f}'.format(plus_cutoff)
        plt.text(x=155+x_pos, y=0.8, s=f'SplAI>{n_cutoff}', color=col, size=20)
        
        #plt.text(tuple_exons[ind][-1][1]+23000, -0.6, 'End: CHROM ' + str(dfs[ind].iloc[0]['CHROM']) + 
        #         ', POS ' + str('{:,.0f}'.format(tuple_exons[ind][-1][1])), size = 25)
        #plt.text((tuple_exons[ind][0][0] - abs(tuple_exons[ind][-1][1]-tuple_exons[ind][0][0])/2)+8000, -0.6, 
        #         'Gene Length: ' + str(abs(tuple_exons[ind][-1][1]-tuple_exons[ind][0][0])/1000) + ' kb', 
        #         size = 25)

    #for spl_pos in splice_vars['POS']:
    #    plt.plot(int(spl_pos), 0.4, marker="o", markersize=2, markerfacecolor="red")
    
    plt.text(gene_start, 0.83, 'Exons', c = 'orange', size = 25, fontweight = 'bold')
    plt.text(gene_start, 0.7, 'Introns', c = 'dodgerblue', size = 25, fontweight = 'bold')

    plt.text(gene_end/2, -0.85, (f'Reading Direction'), size = 20, ha='center')

    plt.text(gene_end/2, y=-0.6, size = 20, ha='center',
             s=(f'Chromosome {chrom}, position {real_start} to {real_end}, '+
                'intronic regions up to 500 bp from exon-intron boundary'))


    plt.gcf().subplots_adjust(bottom = 0.05, top = 0.95, left = 0.01, right = 0.95);

    #plt.savefig(r'figures/' + date + '_' + gene_names[ind].lower() + '_exins_ov.pdf')
    #plt.savefig(r'figures/' + date + '_' + gene_names[ind].lower() + '_exins_ov.jpg')
    plt.savefig(f'00_Andres_SDHB/{gene}_overview_splice_variants.svg')





# short splice overview figure
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

def search_col(inp, df):
    return [i for i in df.columns if inp.lower() in i.lower()]

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

def splai_count_df():
    splai_df = pd.DataFrame()
    splai_df['name'] = vep_dict.keys()
    for i in np.arange(0,1,0.1):
        lst = []
        for name in vep_dict.keys():
            lst.append(len(splai_cutoff(i, vep_dict[name])))
        splai_df[round(i, ndigits=1)] = lst
    return splai_df.sort_values('name').reset_index(drop=True)

def splai_spec_count_df():
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

gtf_sdbh = pd.read_pickle('00_sdhb_seq_dfs/sdhb_gtf.gtf')
gtf_sdbh = gtf_sdbh[['feature','seqname','start','end','strand','frame','gene_id','transcript_id','exon_number']]

sdhb_df = pd.read_csv('00_sdhb_seq_dfs/sdhb_syn.csv')

# variant in exonic or intronic region
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




def spl_columns(df):
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

#synthetic variants saved as pickle
vep_dict = pd.read_pickle('02_vep_results/vep_dict.pickle')

all_vars = pd.DataFrame()
for key in vep_dict.keys():
    all_vars = all_vars.append(vep_dict[key])
all_vars = all_vars.sort_values('POS').reset_index(drop=True)

# define cutoff
cutoff = 0.69

spl_vars = pd.DataFrame()
for key in vep_dict.keys():
    spl_vars = spl_vars.append(splai_cutoff(cutoff,vep_dict[key])).reset_index(drop=True)

spl_vars = spl_columns(spl_vars)

print(len(spl_vars))

spl_vars['pic_col'] = np.nan
spl_vars['y_pos'] = np.nan

for cutoff,col,y_pos in zip([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1], 
                            ['w','w','gray','k','b','g','c', 'y','m','orangered','firebrick'], 
                            [0.42,0.44,0.46,0.48,0.50,0.52,0.54,0.56,0.58,0.6,0.62]):
    spl_vars.loc[((spl_vars['SpliceAI_pred_DS_AG']>=cutoff)|(spl_vars['SpliceAI_pred_DS_AL']>=cutoff)|
                  (spl_vars['SpliceAI_pred_DS_DG']>=cutoff)|(spl_vars['SpliceAI_pred_DS_DL']>=cutoff))&
                 ((spl_vars['SpliceAI_pred_DS_AG']<cutoff+0.1)&(spl_vars['SpliceAI_pred_DS_AL']<cutoff+0.1)&
                  (spl_vars['SpliceAI_pred_DS_DG']<cutoff+0.1)&(spl_vars['SpliceAI_pred_DS_DL']<cutoff+0.1)), 
                 'pic_col']=col
    spl_vars.loc[((spl_vars['SpliceAI_pred_DS_AG']>=cutoff)|(spl_vars['SpliceAI_pred_DS_AL']>=cutoff)|
                  (spl_vars['SpliceAI_pred_DS_DG']>=cutoff)|(spl_vars['SpliceAI_pred_DS_DL']>=cutoff))&
                 ((spl_vars['SpliceAI_pred_DS_AG']<cutoff+0.1)&(spl_vars['SpliceAI_pred_DS_AL']<cutoff+0.1)&
                  (spl_vars['SpliceAI_pred_DS_DG']<cutoff+0.1)&(spl_vars['SpliceAI_pred_DS_DL']<cutoff+0.1)), 
                 'y_pos']=y_pos
spl_vars['y_pos'] = spl_vars['y_pos'].astype(float)

# figure for exons and introns
def overview_fig(gene, cutoff, all_vars = all_vars, gtf=intr_ex_gtf_sdhb, splice_vars = spl_vars):
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
    fig.suptitle('Gene model of $%s$ gene in reading direction with SpliceAI scores (>= %s)' 
                 %(gene.upper(), cutoff), fontsize=25, fontweight='bold')
    
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
    
    x_pos = 150
    for coff,col in zip([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1], 
                                ['w','w','gray','k','b','g','c', 'y','m','orangered','firebrick']):
        if coff >= cutoff:
            x_pos += 750
            n_coff = '{:.1f}'.format(coff)
            #plus_cutoff = float(n_cutoff)+0.1
            #plus_cutoff = '{:.1f}'.format(plus_cutoff)
            plt.text(x=150+x_pos, y=0.8, s=f'SplAI>={n_coff}', color=col, size=20)
        
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
    plt.savefig(f'{gene}_{cutoff}_overview_splice_variants.svg')

overview_fig(gene = 'sdhb', cutoff = 0.7)
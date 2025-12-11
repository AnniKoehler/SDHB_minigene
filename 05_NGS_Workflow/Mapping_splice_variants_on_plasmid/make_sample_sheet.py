from glob import glob
import pandas as pd

AN_folder = '/FASTQ/'
AN_id_t = 'AN123456'
AN_ids = AN_id_t.split(',')
reads1 = []
reads2 = []

for anid in AN_ids:
    reads1.append(glob(AN_folder+f'{anid}*R1*.fastq.gz')[0])
    reads2.append(glob(AN_folder+f'{anid}*R2*.fastq.gz')[0])

samples_df = pd.DataFrame()
samples_df['ID'] = AN_ids
samples_df['forward reads'] = reads1
samples_df['reverse reads'] = reads2
samples_df['Type'] = 'A'
samples_df['Concentration'] = 1

samples_df.to_csv('samples.tsv', sep='\t', index=False)

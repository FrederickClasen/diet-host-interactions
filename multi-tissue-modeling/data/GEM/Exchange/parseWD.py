import sys,csv
import pandas as pd

df = pd.read_csv(sys.argv[1],sep='\t')

for index,row in df.iterrows():
	if 'BOTH' not in row['rxn']:
		rxn = 'EXC_BOTH_'+row['rxn'].split('_')[2]
	else:
		rxn = row['rxn']
	if float(row['UB']) != 0:\
		print(rxn + '\t-' + str(row['UB']) + '\t' + '0')

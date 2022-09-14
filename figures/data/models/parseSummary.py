import sys,os
import pandas as pd

sheets = ['HMRdatabase3_00.xml','Human-GEM.xml','MMR.xml','MMRNetwork.xml','Recon3D_301.xml','iMM1865.xml']

for sheet in sheets:
	df = pd.read_excel('Summary.xlsx',sheet_name=sheet)
	unbalanced = 0
	missing = 0
	for index,row in df.iterrows():
		if row['Status'] == 0:
			if row['C'] != 0 or row['N'] != 0 or row['O'] != 0 or row['S'] != 0 or row['P'] != 0:
				unbalanced += 1
		elif row['Status'] == -1 or row['Status'] == -2:
				missing += 1
		else:
			pass
	print(sheet + '\tunbalanced:' + str(unbalanced) + '\tmissing:' + str(missing))

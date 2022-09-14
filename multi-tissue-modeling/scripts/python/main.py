import os
from collections import defaultdict
import pandas as pd
import numpy as np
import seaborn

######### FUNCTIONS TO PROCESS METAGENOMICS ############

def getFiles(path):
    for file in os.listdir(path):
        if os.path.isfile(os.path.join(path, file)):
            yield file

def getAllTaxa(d):
    taxa = []
    for sample in d:
        for t in d[sample]:
            taxa.append(t)
    return list(set(taxa))

def preprocessMetaphlan(path):
    files = getFiles(path)
    samples = []
    d = defaultdict(dict)
    for file in files:
        sample = file.split('.')[0]
        samples.append(sample)
        f = open(path+file,'r')
        lines = f.readlines()[4:]  # skip first lines
        f.close()
        for line in lines:
            line = line.strip()
            line = line.split('\t')
            d[sample][line[0]] = line[2]
    alltaxa = getAllTaxa(d)
    out = open(path + 'combined.tsv','w')
    out.write('taxa\t')
    for sample in samples:
        out.write(sample + '\t')
    out.write('\n')
    for taxa in alltaxa:
        out.write(taxa + '\t')
        for sample in samples:
            try:
                out.write(str(d[sample][taxa]) + '\t')
            except KeyError:
                out.write(str('0') + '\t')
        out.write('\n')
    out.close()

def subsetMetagenomics(df,annot):
    annot = pd.read_excel(annot, sheet_name='Metagenomics')
    d = dict(zip(annot['ID'], annot['MouseID']))
    select_samples = []
    for index,row in annot.iterrows():
        if row['Tissue'] == 'Faeces':
            select_samples.append(row['ID'])
    df = df[select_samples]
    df = df.rename(columns=d, inplace=False)
    return df

def getTaxonomyLevel(df,level):
    ren = {}
    for index in df.index:
        parts = index.split('|')
        player = parts[-1].split('__')[0]
        if player == level:
            ren[index] = parts[-1].split('__')[1]
    df = df.loc[list(ren.keys())]
    df = df.rename(index=ren, inplace=False)
    return df

######### FUNCTIONS TO PROCESS METABOLOMICS ############
def quantileCutoff(df,percentile):
    q = np.quantile(df,percentile)
    cols = []
    for index,row in df.iterrows():
        for col in df.columns:
            if df.at[index,col] >= q:
                #df.at[index,col] == 0
                cols.append(col)
            else:
                pass
    df = df.drop(list(set(cols)), axis=1)
    return df

def metaboliteAnnot(f):
    dat = pd.read_excel(f,sheet_name='Chemical Annotation',index_col='CHEM_ID')
    d = {}
    for index, row in dat.iterrows():
        d[index] = row['CHEMICAL_NAME']
    return d

def preprocessMetabolomics(f):
    dat = pd.read_excel(f,sheet_name='Log Transformed Data',index_col='PARENT_SAMPLE_NAME')
    annot = metaboliteAnnot(f)
    dat = dat.rename(columns=annot, inplace=False)
    dat = quantileCutoff(dat,0.99)
    return dat

def renameSamples(df,annot):
    annot = pd.read_excel(annot,sheet_name='bloodMetabolomics')
    d = dict(zip(annot['SampleID'],annot['MouseID']))
    df = df.rename(index=d,inplace=False)
    return df

###### FUNCTIONS TO GENERATE SPECIES AND METABOLITE ANNOTATIONS ##########
def generateTaxonomyTable(lineages,out):
    out = open(out,'w')
    out.write('taxa\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies\t\n')
    done = []
    for i in lineages:
        parts = i.split('|')
        if str(parts[-1].split('__')[1]) not in done:
            done.append(str(parts[-1].split('__')[1]))
            out.write(parts[-1].split('__')[1] + '\t')
            for item in parts:
                item = item.split('__')[1].strip()
                out.write(str(item) + '\t')
            out.write('\n')
        else:
            done.append(str(parts[-1].split('__')[1]))
    out.close()

def generateMetAnnotationTable(df,out):
    df = pd.read_excel(df,sheet_name='Chemical Annotation')
    df2 = df[['CHEMICAL_NAME','SUPER_PATHWAY','SUB_PATHWAY']]
    df2.to_csv(out)

if __name__ == '__main__':
    preprocessMetaphlan('/Volumes/lab-anastasioud/home/users/clasenf/Microbiome/metaphlan3/out/')
    metagenomics = pd.read_csv('/Volumes/lab-anastasioud/home/users/clasenf/Microbiome/metaphlan3/out/combined.tsv',sep='\t',index_col='taxa')
    metagenomics = subsetMetagenomics(metagenomics,'designs.xlsx')
    generateTaxonomyTable(metagenomics.index,'taxonomy.tsv')

    metabolomics = preprocessMetabolomics('/Volumes/lab-anastasioud/home/users/clasenf/Metabolomics/Results/KGCO-04-20MD MOUSE SERUM DATA TABLES.XLSX')
    metabolomics = renameSamples(metabolomics,'designs.xlsx')
    # only get metabolomics samples that is in metagenomics
    metabolomics = metabolomics.loc[list(metagenomics.columns)]
    metabolomics = metabolomics.T
    generateMetAnnotationTable('/Volumes/lab-anastasioud/home/users/clasenf/Metabolomics/Results/KGCO-04-20MD MOUSE SERUM DATA TABLES.XLSX','metAnnot.csv')

    #print(metabolomics.head())
    #print(metagenomics.head())
    combined = pd.concat([metagenomics,metabolomics])
    combined = combined.T
    print('running correlation ....')
    correlation_df = combined.corr(method='spearman')
    correlation_df = correlation_df.loc[correlation_df.columns.str.contains('Bacteria'), ~correlation_df.columns.str.contains('Bacteria')]
    correlation_df = correlation_df.fillna(0)
    correlation_df = getTaxonomyLevel(correlation_df,'s')
    correlation_df.to_csv('corr.csv')
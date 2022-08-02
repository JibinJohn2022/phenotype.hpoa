import pandas as pd
import numpy as np
import argparse
import sys
import warnings

warnings.simplefilter(action='ignore')


parser=argparse.ArgumentParser(description="It is pre processing phenotype.hpoa file for kg creation \n")
parser.add_argument('-data','--data', help="phenotype.hpoa file", required=True)
parser.add_argument('-output','--output', help="out put file name", required=True)

args=parser.parse_args()
data=args.data
outpout=args.output

df=pd.read_csv(data,sep="\t",skiprows=4)

#Column description : https://hpo-annotation-qc.readthedocs.io/en/latest/smallfile.html

#Identify all disease associated HPO
Disease_HPO=pd.DataFrame(df['HPO_ID'].unique()).rename(columns={0:"Disease_HPO"})


#Low weitage should give to the hpo terms where 'Qualifier'=='NOT'
df[df['Qualifier']=='NOT']

#Frequency: df['Frequency']
Frequency={'HP:0040284':"Very rare(1-4%)",'HP:0040283':"Occasional(5-29%)",
           'HP:0040282':"Frequent(30-79%)",'HP:0040281':"Very frequent(80-99%)",
           'HP:0040280':"Obligate(100%)"}

###########-------------Frequency synchronisation---------------------------------- 
Frequency_1=df[df['Frequency'].fillna("0").str.contains("/")]

Frequency_1[["P_present","Total_reported"]]=Frequency_1['Frequency'].str.split("/",expand=True).astype(int)
Frequency_1['ReportedPercentage']=Frequency_1["P_present"]/Frequency_1['Total_reported']*100

conditions = [(Frequency_1['ReportedPercentage']<5),
              ((Frequency_1['ReportedPercentage']>=5) &(Frequency_1['ReportedPercentage']<30)),
              ((Frequency_1['ReportedPercentage']>=30) &(Frequency_1['ReportedPercentage']<80)),
              ((Frequency_1['ReportedPercentage']>=80) &(Frequency_1['ReportedPercentage']<100)),
              (Frequency_1['ReportedPercentage']==100)]
choices = ["HP:0040284","HP:0040283","HP:0040282","HP:0040281","HP:0040280"]
Frequency_1['HPO_Reported_Frequency'] = np.select(conditions, choices)


#HPO terms
Frequency_2=df[df['Frequency'].fillna("0").str.contains("HP")]
Frequency_2['HPO_Reported_Frequency']=Frequency_2['Frequency']

#Percentage
Frequency_3=df[df['Frequency'].fillna("0").str.contains("%")]
Frequency_3['TotalPercentage']=Frequency_3['Frequency'].str.replace("\%","").astype(float)
conditions = [(Frequency_3['TotalPercentage']<5),
              ((Frequency_3['TotalPercentage']>=5) &(Frequency_3['TotalPercentage']<30)),
              ((Frequency_3['TotalPercentage']>=30) &(Frequency_3['TotalPercentage']<80)),
              ((Frequency_3['TotalPercentage']>=80) &(Frequency_3['TotalPercentage']<100)),
              (Frequency_3['TotalPercentage']==100)]
choices = ["HP:0040284","HP:0040283","HP:0040282","HP:0040281","HP:0040280"]
Frequency_3['HPO_Reported_Frequency'] = np.select(conditions, choices)

#No frequency
Frequency_4=df[(df['Frequency'].fillna("0").str.contains("0")) & ~(df['Frequency'].fillna("0").str.contains("HP|/|%"))]
Frequency_4['HPO_Reported_Frequency']=Frequency_4['Frequency']

df2=pd.concat([Frequency_1,Frequency_2,Frequency_3,Frequency_4])
df2.drop(['TotalPercentage','Total_reported','ReportedPercentage','P_present'],axis=1,inplace=True)

###########-------------Frequency synchronisation----------------------------------

df2['SubId']=df2['#DatabaseID']+"_"+df2['HPO_ID'].str.replace("\:","_")
df2['SubId']=df2['SubId'].str.replace(":","_")


#----------------------####################################################################

df2['Reference']=df2['Reference'].str.split(";")
df2=df2.set_index(['#DatabaseID','DiseaseName','Qualifier','HPO_ID','Evidence','Onset','SubId','Frequency', 'Sex', 'Modifier', 
         'Aspect','Biocuration','HPO_Reported_Frequency']).apply(pd.Series.explode).reset_index()

for column in ['PMID', 'OMIM','ORPHA','DECIPHER','ISBN','http']:
    df2[column]=np.where(df2['Reference'].str.contains(column),df2['Reference'],np.nan)

df3=df2[['#DatabaseID','DiseaseName','Qualifier','HPO_ID','Evidence','Onset','Sex','Biocuration',
        'Modifier','Aspect','Biocuration','Frequency','HPO_Reported_Frequency','Reference','PMID',
        'OMIM','ORPHA','DECIPHER','ISBN','http','SubId']]

df3['OMIM']=df3['OMIM'].str.replace("OMIM:","").str.strip()
df3['ORPHA']=df3['ORPHA'].str.replace("ORPHA:","Orphanet_").str.strip()
df3['PMID']=df3['PMID'].str.replace("PMID:","").str.strip()
df3['ISBN']=df3['ISBN'].str.split(":",expand=True)[1]



#############----------Datbase id ------------------------------
df3['DatabaseID-OMIM']=np.where(df3['#DatabaseID'].str.contains('OMIM'),df3['#DatabaseID'],np.nan)
df3['DatabaseID-ORPHA']=np.where(df3['#DatabaseID'].str.contains('ORPHA'),df3['#DatabaseID'],np.nan)
df3['DatabaseID-DECIPHER']=np.where(df3['#DatabaseID'].str.contains('DECIPHER'),df3['#DatabaseID'],np.nan)

df3['DatabaseID-OMIM']=df3['DatabaseID-OMIM'].str.replace("OMIM:","").str.strip()
df3['DatabaseID-ORPHA']=df3['DatabaseID-ORPHA'].str.replace("ORPHA:","Orphanet_").str.strip()
df3=df3.replace(r'^\s*$', np.nan, regex=True)
df3.drop("#DatabaseID",axis=1,inplace=True)


#####------------Aspect-----------------------------------------
df3['Aspect_Ontology']=df3['Aspect'].replace({'P': 'HP_0000118', 'I': 'HP_0000005', 'C': 'HP_0031797', 'M': 'HP_0012823'})


###----------Onset---------------------------------------------------
df3['Onset']=df3['Onset'].str.replace(":","_")

## ---------'Modifier'-----------------
df3['Modifier']=df3['Modifier'].str.replace(":","_")


###---------Sex-------------------------------------------------
df3['Sex']=df3['Sex'].replace({'male':'MALE','female':'FEMALE'})

###----------HPO_ID------------------------------------------
df3['HPO_ID']=df3['HPO_ID'].str.replace("\:","_")

df3.to_csv("Disease_to_HPO.csv",index=None)

###---------HPO_Reported_Frequency--------------------
df3['HPO_Reported_Frequency']=df3['HPO_Reported_Frequency'].str.replace("\:","_")

#--------------------------------------------------------

def split_dataframe(df, chunk_size = 10000): 
    chunks = list()
    num_chunks = len(df) // chunk_size + 1
    for i in range(num_chunks):
        chunks.append(df[i*chunk_size:(i+1)*chunk_size])
    return chunks

df_split2df_split2=split_dataframe(df3,chunk_size = 20000)

for i in range(len(df_split2df_split2)):
        df_split2df_split2[i].to_csv(outpout+"_"+str(i)+".csv",index=None)

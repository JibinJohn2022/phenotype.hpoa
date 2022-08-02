
#!/usr/bin/python3

import pandas as pd
import numpy as np
import argparse
import sys
from rdflib import Namespace
from rdflib.namespace import OWL,RDF, RDFS,SKOS
from rdflib import Graph, URIRef, Literal, BNode
import warnings
warnings.filterwarnings('ignore')

parser=argparse.ArgumentParser(description="It is for KG creaation from phenotype.hpoa file\n")
parser.add_argument('-data','--data', help="processed phenotype.hpoa file in csv format", required=True)
parser.add_argument('-schema','--schema', help="schema file", required=True)
parser.add_argument('-output','--output', help="schema file", required=True)


args=parser.parse_args()
data=args.data
schema=args.schema
outpout=args.output

g = Graph()
nm = g.namespace_manager

gvaont = Namespace("http://semanticwebindia.in/GVA/ont/");prefix = "gvaont";nm.bind(prefix, gvaont)
obo = Namespace("http://purl.obolibrary.org/obo/");prefix = "obo";nm.bind(prefix, obo)

hgnc=Namespace("http://identifiers.org/hgnc/");prefix = "hgnc";nm.bind(prefix, hgnc)
ncbi=Namespace("https://www.ncbi.nlm.nih.gov/gene/") ;prefix = "ncbi" ;nm.bind(prefix, ncbi)

omimbio=Namespace("http://purl.bioontology.org/ontology/OMIM/") ;prefix = "omimbio" ;nm.bind(prefix, omimbio)
orphabio=Namespace("http://www.orpha.net/ORDO/") ;prefix = "orphabio" ;nm.bind(prefix, orphabio)


Schema=pd.read_csv(schema,encoding='cp1252')
#Schema.head()
Data=pd.read_csv(data,sep=",",dtype=str)

ListcOlumns=list(Schema["Subjectcolumns"].unique())+list(Schema["ObjectColumn"].unique())
ListcOlumns=[x for x in ListcOlumns if x != 'NotApplicable']


for column in Data.columns:
    Data[column] = Data[column].str.strip()

for column in Schema.columns:
    Schema[column] = Schema[column].str.strip()


Data=Data.fillna("NA")
Schema=Schema.fillna("NA")
#Data=Data.head(n=10000)
Data=Data[ListcOlumns]
g = Graph()
Subjeclist=set(Data[",".join(Schema['Subjectcolumns'].unique())].to_list())

n=0
for subject in Subjeclist:
    n=n+1
    print(n)
    for i, j in Schema.iterrows():
        #Subject
        SPrefix=j['Subject']
        SubjectValue=subject
        Subject = URIRef(SPrefix+SubjectValue)
        
        #Predicate
        Predicate=j['Predicate']
        
        #Object identification
        #1a) if ObjectColumn equal to NotApplicable
        if  j['ObjectColumn']=="NotApplicable" and j['object_ObjectPrefix']!="Literal()":
            if j['object_ObjectPrefix']!="NA":
                if  (j['ObjectColumn']=="NotApplicable") and ("http" in j['object_ObjectPrefix']):
                    object_ObjectPrefix=j['object_ObjectPrefix']
                    g.add((Subject, URIRef(Predicate), URIRef(object_ObjectPrefix)))
                elif (j['ObjectColumn']=="NotApplicable") and ~("http" in j['object_ObjectPrefix']):
                    object_ObjectPrefix=j['object_ObjectPrefix']
                    g.add((Subject, URIRef(Predicate), Literal(object_ObjectPrefix)))
            elif j['object_ObjectPrefix']=="NA":
                if  "http" in j['object_ObjectPrefix']:
                    g.add((Subject, URIRef(Predicate), URIRef(ObjectColumn)))
        elif  j['ObjectColumn']=="NotApplicable" and j['object_ObjectPrefix']=="Literal()":
            object_ObjectPrefix=j['object_ObjectPrefix']
            object_ObjectPrefix=Literal(object_ObjectPrefixSufix)
            g.add((Subject, URIRef(Predicate), object_ObjectPrefix))
        #1a) if ObjectColumn not equal to NotApplicable
        if  j['ObjectColumn']!="NotApplicable" and j['object_ObjectPrefix']=="Literal()":
            objects=set(Data[Data[j['Subjectcolumns']]==subject][j['ObjectColumn']].to_list())
            objects=[x for x in objects if str(x) != 'NA']
            if len(objects) >0:
                for Object in objects:
                    object_ObjectPrefix=Literal(Object)
                    g.add((Subject, URIRef(Predicate), object_ObjectPrefix))
        if  j['ObjectColumn']!="NotApplicable" and j['object_ObjectPrefix']!="Literal()":
            objects=set(Data[Data[j['Subjectcolumns']]==subject][j['ObjectColumn']].to_list())
            #If object_objectPrefix not equal to "NA"
            if j['object_ObjectPrefix']!="NA":
                objects=[x for x in objects if str(x) != 'NA']
                if len(objects) >0:
                    for Object in objects:
                        object_ObjectPrefix = URIRef(j['object_ObjectPrefix']+str(Object))
                        g.add((Subject, URIRef(Predicate), object_ObjectPrefix))
            elif j['object_ObjectPrefix']=="NA":
                objects=[x for x in objects if str(x) != 'NA']
                if len(objects) >0:
                    for Object in objects:
                        object_ObjectPrefix = URIRef(str(Object))
                        g.add((Subject, URIRef(Predicate), object_ObjectPrefix))



#for s, p, o in g:
#    print(s, p, o)

s = g.serialize(format="ttl",destination=outpout+".ttl")
#print(s)
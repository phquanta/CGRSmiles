# -*- coding: utf-8 -*-
"""
Created on Wed Sep 30 22:16:03 2020

@author: Andrei
"""


import pandas as pd
from html.parser import HTMLParser
from bs4 import BeautifulSoup
import numpy as np
import copy
import rdkit.Chem.rdChemReactions as cR
from rdkit.Chem import AllChem
from rdkit import Chem
from pubchempy import *
import CGRtools as cgr
from io import StringIO

#col_list = ["reaction"]
dbs_to_read=["CID_784_rhea.csv","CID_784_pathwayreaction.csv"]
seps=["a> =","⟶"]
cols=["htmlequation","reaction"]
unmapped_smarts=set()
last3Chars='</a>'
def getMolFromCID(cid):
     #print("CID:",cid)
     cid_smile= (Compound.from_cid(cid)).canonical_smiles
     #print(cid_smile,type(cid_smile))
     #print("\n")
     return Chem.MolFromSmiles(cid_smile)     

for ind,fn in enumerate(dbs_to_read):
    print(cols[ind])
    df = pd.read_csv(fn, usecols=[cols[ind]])
    cntt=0
    for x in df[cols[ind]]:
     try:   
        #print("STARTRTTTTTTTTTTTTTTTTTTTTTTTTTTTTtt")
#        print(x)
        #print("SEPPPPPP:",seps[ind])
        #print("\n")
        
        All=x.split(seps[ind])
        #print("SEEEEEEEEEEEEEEE")
        #print("SEEEEEEEEEEEEEEE")
        #print("SEEEEEEEEEEEEEEE")
        #print(All)
        
        reactants=All[0].split('+')
        products=All[1].split('+')
        lenAll=len(reactants)+len(products)
        xx=copy.copy(x)
        if (xx.strip())[-4:]!=last3Chars:
            print("###################################")   
            
            print("########### Chopped reaction STRING ################\n\n")
            
            print(x+"\n\n")
            
            continue
        if 'protein' in x:
            #print("FOIUND")
            continue
        cntt+=1
        cntR=0
        cntP=0
        chem_reaction=cR.ChemicalReaction()
        
        for z in reactants:
            soup = BeautifulSoup(z, 'html.parser')
            for link in soup.find_all('a'):
                a=link.get('href')
                cid=(a.split('compound/'))
#                print(a,"AAA")
                mol=getMolFromCID(cid[1])
#                print(mol)
                chem_reaction.AddReactantTemplate(mol)
#                print(a,cid[1],len(cid))
                cntR+=1
    
    
        for z in products:
            soup = BeautifulSoup(z, 'html.parser')
            for link in soup.find_all('a'):
                a=link.get('href')
                cid=(a.split('compound/'))
                mol=getMolFromCID(cid[1])
#                print(mol)

                chem_reaction.AddProductTemplate(mol)
#                print(a,cid[1],len(cid))
                cntP+=1
        
        counter1=len(BeautifulSoup(x, 'html.parser').find_all('a'))
        counter3=len(BeautifulSoup(x, 'html.parser').find_all('.a'))
        counter2=cntR+cntP
        
        if counter1!=counter2 or counter1!=lenAll:
            print("Errrorneous reaction SKIPPPING")
            continue
        smarts_string=Chem.rdChemReactions.ReactionToSmiles(chem_reaction)
        
        #nextCGRObj=next(cgr.files.SMILESRead(StringIO(smarts_string)))
        
        if len(smarts_string)<150:
          if smarts_string not in unmapped_smarts:
            print("SM_String:",smarts_string)    
            unmapped_smarts.add(smarts_string)
        
        
    
    
        #print(All[0],"##################33", All[1])
        #print("SEEEEEEEEEEEee")
     except Exception as e:
         print("Error Reading ---- Continuing")
         continue
    
    print("CountT",cntt)
    
unmapped_smarts=list(unmapped_smarts)    
with open('Extra_OO_Reactions_Unmapped_FineTuning.smi', 'w') as f:
            for item in unmapped_smarts:
                f.write("%s\n" % item)            
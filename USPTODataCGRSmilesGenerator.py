#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  8 19:47:32 2020

@author: phquanta
"""


import rdkit

from rdkit.Chem import AllChem

from rdkit import Chem, RDConfig
from rdkit.Chem import MolFromSmiles, AllChem
from rdkit.Chem.rdChemReactions import PreprocessReaction
import CGRtools as cgr
from pathlib import Path
import numpy as np
from io import StringIO
import os

#reactionsFile="/mnt/hgfs/H/DataSets/USPTO/USPTO_AllReactions_SMARTS.smi"
reactionsFile="H:\\DataSets\\USPTO\\USPTO_AllReactions_SMARTS.smi"
#reactionsFile="/mnt/hgfs/H/DataSets/JinsData/test1.smi"
maxLenCGR=300
CGRSmileFileOut=f"USPTOCGRSmiles{maxLenCGR}.smi"

fname = Path(CGRSmileFileOut)
if fname.exists():
    os.remove(fname)
    # file exists


generator=cgr.files.SMILESRead(reactionsFile)




def GetAllfromRDkit(reactionStr,i):
  #try:  
    
    rxn = AllChem.ReactionFromSmarts(str(nextReactionObj))
    #reactantsRD1=[rdkit.Chem.MolToSmiles(i) for i in rxn.GetReactants()]
    #print(reactantsRD1)
    
    rm=[]
    ra=[]
    pm=[]
    ReactantsMol= rxn.GetReactants()
    
    ReagentstsMol= rxn.GetAgents()
    ProductsMol= rxn.GetProducts()
    
    
    for mol in ReactantsMol:
         Chem.SanitizeMol(mol,Chem.SanitizeFlags.SANITIZE_FINDRADICALS|Chem.SanitizeFlags.SANITIZE_KEKULIZE|
                          Chem.SanitizeFlags.SANITIZE_SETAROMATICITY|Chem.SanitizeFlags.SANITIZE_SETCONJUGATION|
                          Chem.SanitizeFlags.SANITIZE_SETHYBRIDIZATION|Chem.SanitizeFlags.SANITIZE_SYMMRINGS
                          ,catchErrors=True)
         rm.append(mol)
         
         
    for mol in ReagentstsMol:     
         Chem.SanitizeMol(mol,Chem.SanitizeFlags.SANITIZE_FINDRADICALS|Chem.SanitizeFlags.SANITIZE_KEKULIZE|
                          Chem.SanitizeFlags.SANITIZE_SETAROMATICITY|Chem.SanitizeFlags.SANITIZE_SETCONJUGATION|
                          Chem.SanitizeFlags.SANITIZE_SETHYBRIDIZATION|Chem.SanitizeFlags.SANITIZE_SYMMRINGS
                          ,catchErrors=True)
         ra.append(mol)
    for mol in ProductsMol:     
         Chem.SanitizeMol(mol,Chem.SanitizeFlags.SANITIZE_FINDRADICALS|Chem.SanitizeFlags.SANITIZE_KEKULIZE|
                          Chem.SanitizeFlags.SANITIZE_SETAROMATICITY|Chem.SanitizeFlags.SANITIZE_SETCONJUGATION|
                          Chem.SanitizeFlags.SANITIZE_SETHYBRIDIZATION|Chem.SanitizeFlags.SANITIZE_SYMMRINGS
                          ,catchErrors=True)
         pm.append(mol)
      
         
        
        
    
    
    
    reactantsRD=[rdkit.Chem.MolToSmiles(i) for i in rm]
    reagentsRD=[rdkit.Chem.MolToSmiles(i) for i in ra]
    productsRD=[rdkit.Chem.MolToSmiles(i) for i in pm]
    
#    if i==234:
#         print("AAAAAAAAAAAAa",reactantsRD,reagentsRD,productsRD)
#         raise Exception()
    if reactantsRD==productsRD or reactantsRD==None or productsRD==None or reactantsRD==[] or productsRD==[]:
        return [[],[],[]]
   
    
  #except:
  #    return [[],[],[]]
    return [reactantsRD,reagentsRD,productsRD]



ll=set()

strs=[]
ll1=[]



length=[]

with open(CGRSmileFileOut, 'a') as the_file:
    #the_file.write('Hello\n')

 
#for i in range(1000000):
# try:
    #reactor = cgr.CGRReactor(next(SMILESread(StringIO('[C:1](=[O:2])-[O:3]>>[C:1](=[O:2])'))),delete_atoms=True)
   #
   #print("HERE0")
   
  while True:
   try:
    nextReactionObj=next(generator)
    if type(nextReactionObj)!=cgr.containers.cgr.CGRContainer:
        nextReactionObj.standardize()
        
    # isStn=nextReactionObj.standardize()
        cgrContainer=~nextReactionObj
    else:
        cgrContainer=nextReactionObj
        
   #aa,bb,dd=GetAllfromRDkit(str(nextReactionObj),i)
   #print(i,aa,bb,dd)
   # if str(cgrContainer)  not in ll1:
   #     ll1.append(str(cgrContainer))
   #     if '>-^>' in str(cgrContainer):
   #         print(cgrContainer)
       
   
    #nextReactionStr=str(nextReactionObj)
    nextReactionStr=str(cgrContainer)
    #print("reaction STR:",nextReactionStr)
    if len(nextReactionStr)>maxLenCGR:
        continue
    
    if nextReactionStr not in strs:
        strs.append(nextReactionStr)
        length.append(len(nextReactionStr))
        print("Reaction CGR:",nextReactionStr)
    #if nextReactionStr not in 
    
        the_file.write(nextReactionStr+'\n')
   
   
   #reactor = cgr.CGRReactor(next(cgr.SMILESRead(StringIO(nextReaction))),delete_atoms=True)
   #reactor(reactants)
   #print("HERE1")
   #reactantsCGR=[str(i) for i in nextReactionObj.reactants()]
   #reactantsRD=[rdkit.Chem.MolToSmiles(i) for i in rxn.GetReactants()]
   #print("HERE2")
   #reagentsRD=[rdkit.Chem.MolToSmiles(i) for i in rxn.GetAgents()]
   #print("HERE3")
   #productsRD=[rdkit.Chem.MolToSmiles(i) for i in rxn.GetProducts()]
   #print("HERE")
   #print(reactantsRD)
   #print(reagentsRD)
   #print(productsRD)
   
   except StopIteration:
       print("########### STOP ITERATION ##################")
       break
   
   except Exception as e:
       print("########### SOME ERROR ##################")
       print()
       
       print(e)
       
       
       print("Error")


print("mean of Seq. lenght:", np.mean(length))
print("min of Seq. lenght:", np.min(length))
print("max of Seq. lenght:", np.max(length))
   
   
 #    pass
#Chem.SmilesMolSupplier

# -*- coding: utf-8 -*-
"""
Created on Sun Oct 18 10:22:29 2020

@author: Andrei
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

reactionsFile="CGRsGenSMILENew_3_Clean.smi"
oldReactions="FineTuneCGR_OO.smi"
oldCGRs=[]


def isValidReaction(reactionObj):
    valid=False
    reactionMols=list(reactionObj.reactants)
    reactionMols.extend(reactionObj.products)
    reactionMols.extend(reactionObj.reagents)
    
    for i in reactionMols:
        try:
            i.kekule()
            
            lst=i.check_valence()
            if len(lst)==0 and i.check_thiele():
                valid=True
            
        except Exception as e:
            print("ERROR IN isValidReactrion")
            return False
    return valid   



def getCGRs(reactionsFile):

  generator=cgr.files.SMILESRead(reactionsFile)
  ll=set()
  cgrObjs=set()
  ll1=[]

  length=[]
  cnt=0

  cntAll=0 
  while True:
   try:
    nextReactionObj=next(generator)
    #print(str(nextReactionObj),type(nextReactionObj))
    cntAll+=1
    #print("All, len of current list:",cntAll,print(len(ll)))
    if type(nextReactionObj)!=cgr.containers.cgr.CGRContainer:
         nextReactionObj.standardize()
         cgrContainer=~nextReactionObj
         decomposed = cgr.ReactionContainer.from_cgr(cgrContainer)
 #        print("decomposed",decomposed)
         if isValidReaction(decomposed):
             cnt+=1
             if cnt%2==0:
                 print(f"done {cnt} SMILES From TRAIN.TXT")
             nextReactionStr=str(cgrContainer)
             ll.add (nextReactionStr)
             cgrObjs.add(cgrContainer)

         print()
    else:
         cgrContainer=nextReactionObj
         decomposed = cgr.ReactionContainer.from_cgr(cgrContainer)
     #    print("decomposed",decomposed)
         if isValidReaction(decomposed):
             cnt+=1
             if cnt%2==0:
                 print(f"done {cnt} SMILES FROM CONVERTED CGR FILE ALREADY")
             nextReactionStr=str(cgrContainer)
             ll.add (nextReactionStr)
             cgrObjs.add(cgrContainer)
         else:
             continue
   except StopIteration:
       print("########### STOP ITERATION ##################")
       print(cnt)
       break
   
   except Exception as e:
       print("########### SOME ERROR ##################")
       print()
       
       print(e)
       
       
       print("Error")
  return [ll,cgrObjs]


[CGRs,objs]=getCGRs(reactionsFile)
[oldCGRs,olb_objs]=getCGRs(oldReactions)

#for i in oldCGRs:
    

decomposed=[]
Drawing=True
from rdkit.Chem import Draw   
cntt=0
cnttAll=0
for x in objs:
    if x not in oldCGRs:
    #  try:
        
        decomposed.append(str(cgr.ReactionContainer.from_cgr(x)))
        
        cnttAll+=1
    #    decomposed.explicify_hydrogens();e
        if Drawing:
            print(str(decomposed))
            m = AllChem.ReactionFromSmarts(decomposed[-1])
        
            if m is not None:
                cntt+=1  
                #Draw.ReactionToImage(m)
                print(cntt)
    #    decomposed.clean2d()
        #decomposed
    #  except:
     #     pass
    
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

import pickle
#reactionsFile="/mnt/hgfs/H/DataSets/USPTO/USPTO_AllReactions_SMARTS.smi"6
#reactionsFile="train.txt"
#reactionsFile="/mnt/hgfs/H/DataSets/JinsData/test1.smi"
reactionsFile="CGRSmiles156.smi"
#CGRSmileFileOut="CGRSmiles156.smi"
#CGRSSmilesGenerated=["RNNGeneratedReactioins.dat","RNNGeneratedReactioinsStartingOO.dat","BIMODAL_Fixed_Aug1_T_07.smi"]
CGRSSmilesGenerated=["RNNGeneratedReactioins.dat","RNNGeneratedReactioinsStartingOO.dat","..\\CleanData\\RNN_5LSTM_512Each_138Epochs_And_BiMODAL_Fixed\\General_WithBIMODAL_CGRsGenSMILE_2_Clean.smi", "GeneralGeneratedReactioinsStarting_AfterFineTuningOn_OO_Reactions.smi","BiModaL_Aug1_FineTuned_OO.smi","BiModaL_Aug1_FineTuned_OO_1f_ep4.smi","BiModal_30K_Generated_FromGuelphComputer_256.smi","BiModaL_Aug5_FineTuned_OO_Guelph_1f_ep19.smi"]
CGRFineTuneOO="FineTuneCGR_OO.smi"
#CGRSSmilesGenerated=["BiModal_30K_Generated_FromGuelphComputer_256.smi"]


#CGRSSmilesGenerated=["RNNGeneralGeneratedReactioinsNew.dat","RNNGeneralGeneratedReactioinsStartingOONew.dat","..\\CleanData\\RNN_5LSTM_512Each_138Epochs_And_BiMODAL_Fixed\\General_WithBIMODAL_CGRsGenSMILE_2_Clean.smi"]

#CGRSSmilesGenerated=["BIMODAL_Fixed_Aug1_T_07.smi"]

lengthCorrect=[0 for x in range(len(CGRSSmilesGenerated))]
lengthCorrectPercentage=[0 for x in range(len(CGRSSmilesGenerated))]

lengthCorrect_l=[0 for x in range(len(CGRSSmilesGenerated))]
lengthCorrectPercentage_l=[0 for x in range(len(CGRSSmilesGenerated))]



notPickled=True
notPickledTrained=False
n_of_generatedCGRs=30000.

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


def getReactionCentersDictionary(CGRSmiles,cgrs,Verbose=False):
  rcs={}
  for n,cgr_obj in enumerate(cgrs):
   try:
        #print(n,cgr_obj) 
        reactionObj = cgr.ReactionContainer.from_cgr(cgr_obj)
       #reactionObj=preparer.decompose(nextCGRObj)
       #reactionObj=preparer.decompose(nextCGRObj)

        isValid=isValidReaction(reactionObj)
        if isValid:
            rC=cgr_obj.centers_list
            
            for i in range(len(rC)):
                rc1 = cgr_obj.substructure(rC[i], as_query=True) 
                #if (';-' in str(rc1)) and (';+' in str(rc1)):
                #    CntChargedRCs+=1 
                 #hashes.add(rc1.__hash__())
                rcs[rc1.__hash__()]=str(rc1)
                if len(rC)>1:
                    if Verbose:    
                        print(rC," ----- ", rC[i],str(rc1))

                
                #if rc1.__hash__() in  rcs:
                #    print("FOUNDDDDDDDDDDDDDDDDDDDDDD:",str(rc1),rcs[rc1.__hash__()])

            #print("Reaction Centers Number OF:",len(rcs))
            #print(rc1)
            if (n%10)==0 and n>1:
                print(f"Done : {n} Smiles, len(RCS): {len(rcs)} ")
            #Cnt_validCGR+=1
        else:
            print("invalidDDDD")

        
        
        if Verbose:
         print("###################")
         print("")
         
         print(str(cgr_obj))
         print(str(cgr_obj) in CGRSmiles)
         print("#############")
   except StopIteration:
       print("########### STOP ITERATION or END of LOOP ##################")
   except Exception as e:
       print("########### SOME ERROR Detected ##################")
       print(e)
       print("Error")
  return rcs     





def getCGRs(reactionsFile):

  generator=cgr.files.SMILESRead(reactionsFile)
  ll=set()
  cgrObjs=set()
  ll1=[]
  cgrObjs1=[]

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
             ll1.append(nextReactionStr)
             cgrObjs1.append(cgrContainer)

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
             ll1.append(nextReactionStr)
             cgrObjs1.append(cgrContainer)

             print(nextReactionStr)
         else:
             continue
   except StopIteration:
       print("########### STOP ITERATION ##################")
       print(cnt)
       #continue
       break
   
   except Exception as e:
       print("########### SOME ERROR ##################")
       print()
       
       print(e)
       
       
       print("Error")
       continue
  return [ll,cgrObjs,ll1,cgrObjs1]


cgrGen=[]
objGen=[]
Unique=[]

cgrGen_l=[]
objGen_l=[]
Unique_l=[]


if notPickledTrained:
    [CGRs,objs,_,_]=getCGRs(reactionsFile)
    
       
    with open('CGRsTrainNew.pkl', 'wb') as f:
            pickle.dump(CGRs, f)
    with open('CGRsTrainObjsNew.pkl', 'wb') as f:
            pickle.dump(objs, f)
            
    with open('CGRTrainObjsSMILENew_Clean.smi', 'w') as f:
        for item in CGRs:
            f.write("%s\n" % item)        
            
else:
    CGRs = pickle.load( open( "CGRsTrainNew.pkl", "rb" ) )
    objs = pickle.load( open( "CGRsTrainObjsNew.pkl", "rb" ) )
    
    

if notPickled:
      for i,fn in  enumerate(CGRSSmilesGenerated):
        print(fn)
        [cgrR,obj,cgrR_l,obj_l]=getCGRs(fn)
        with open(f'CGRsGenNew_{i}.pkl', 'wb') as f:
                pickle.dump(cgrR, f)
        with open(f'CGRsGenObjsNew_{i}.pkl', 'wb') as f:
                pickle.dump(obj, f)        
        with open(f'CGRsGenSMILENew_{i}_Clean.smi', 'w') as f:
            for item in cgrR:
                f.write("%s\n" % item)        
        cgrGen_l.append(cgrR_l)
        objGen_l.append(obj_l)
        cgrGen.append(cgrR)
        objGen.append(obj)
        
        
        lengthCorrect[i]=len(cgrR)
        lengthCorrect_l[i]=len(cgrR_l)

else:
    for i,fn in  enumerate(CGRSSmilesGenerated):
         cgrR=pickle.load( open( f'CGRsGenNew_{i}.pkl', "rb" ) )
         obj=pickle.load( open( f'CGRsGenObjsNew_{i}.pkl', "rb" ) )
         
         cgrGen.append(cgrR)
         objGen.append(obj)
         lengthCorrect[i]=len(cgrR)
         lengthCorrect_l[i]=len(cgrR_l)

CGRs=list(CGRs)
objs=list(objs)

[CGRsOO,objsOO,_,_]=getCGRs(CGRFineTuneOO)
CGRsOO=list(CGRsOO)
objsOO=list(objsOO)



    
for i,(cg,ob)  in  enumerate(zip(cgrGen,objGen)):
     print(i)
     cgrGen[i]=list(cg)
     objGen[i]=list(ob)
     


lengthCorrectPercentage=[x/n_of_generatedCGRs*100. for x in lengthCorrect]
lengthCorrectPercentage_l=[x/n_of_generatedCGRs*100. for x in lengthCorrect_l]

for lst  in  cgrGen:
        Unique.append([x for x in lst  if x not in CGRs])



rcsAll=getReactionCentersDictionary(CGRs,objs)
rcsAllOO=getReactionCentersDictionary(CGRsOO,objsOO)


rcsGen=[{} for i in  range(len(CGRSSmilesGenerated))]



####################### For Depicting Pics from Generated OO reactions #########################
rcsOOreactionsDataset=[]

for i in objsOO:
     if isValidReaction(cgr.ReactionContainer.from_cgr(i)):
         
         rcsOOreactionsDataset.append([str(i),str(cgr.ReactionContainer.from_cgr(i))])
         

rcsGenOO={}
reactionGensOO=[]
reactionGensOOAll=[]
reactionGensOOAll1=[]
cgrSGenOO=[]
#pp=7
pp=7
elementsOO=[]
for elem,sm in zip(objGen[pp],cgrGen[pp]):
        #print(sm,str(elem))
        reaction_center=getReactionCentersDictionary([sm],[elem])
        keys=[z for z in reaction_center.keys()]
        if isValidReaction(cgr.ReactionContainer.from_cgr(elem)):
            if 'OO'  in str(cgr.ReactionContainer.from_cgr(elem)):
                if elem not in objsOO:
                    reactionGensOOAll.append(str(cgr.ReactionContainer.from_cgr(elem)))
                    reactionGensOOAll1.append([str(cgr.ReactionContainer.from_cgr(elem)),sm])
        if     len(keys)==1 and keys[0] not in rcsAllOO:
            rcsGenOO.update(reaction_center)
            if isValidReaction(cgr.ReactionContainer.from_cgr(elem)):
                #print(type(elem))
                #elem.clean2d()
                elementsOO.append(elem)
                reactionGensOO.append([keys[0],str(cgr.ReactionContainer.from_cgr(elem)),sm])
                #reactionGensOO.append([keys[0],str(cgr.ReactionContainer.from_cgr(elem))])
                
                cgrSGenOO.append([keys[0],sm])
            print(len(rcsGenOO))
##################### End of Depicting

pp=2
reactionGensOO_fromAllGenerated=[]
reactionGensOO_fromAllGenerated1=[]
reactionGensOO_fromAllGenerated_RCS={}
reactionGensOO_fromAllGenerated_RCS1=[]
elementsOO_fromAllGenerated=[]
for elem,sm in zip(objGen[pp],cgrGen[pp]):
        #print(sm,str(elem))
        reaction_center=getReactionCentersDictionary([sm],[elem])
        keys=[z for z in reaction_center.keys()]
        if isValidReaction(cgr.ReactionContainer.from_cgr(elem)):
            if 'OO'  in str(cgr.ReactionContainer.from_cgr(elem)):
                if elem not in objsOO:
                    reactionGensOO_fromAllGenerated.append(str(cgr.ReactionContainer.from_cgr(elem)))
                    reactionGensOO_fromAllGenerated1.append([str(cgr.ReactionContainer.from_cgr(elem)),sm])
                    print(len(reactionGensOO_fromAllGenerated))
        if   len(keys)==1 and keys[0] not in rcsAllOO and 'OO'  in str(cgr.ReactionContainer.from_cgr(elem)):
            reactionGensOO_fromAllGenerated_RCS.update(reaction_center)
            if isValidReaction(cgr.ReactionContainer.from_cgr(elem)):
                #print(type(elem))
                #elem.clean2d()
                elementsOO_fromAllGenerated.append(elem)
                reactionGensOO_fromAllGenerated_RCS1.append([keys[0],str(cgr.ReactionContainer.from_cgr(elem)),sm])
                #reactionGensOO.append([keys[0],str(cgr.ReactionContainer.from_cgr(elem))])
                
                #cgrSGenOO.append([keys[0],sm])
            print(len(reactionGensOO_fromAllGenerated_RCS1))


reactionsOO_withS=[]
for i in reactionGensOOAll1:
    if 'S' in i[0] and   len(i[0])<70:
        reactionsOO_withS.append(i)




reactionsOO_withI=[]
for i in reactionGensOOAll1:
        #if  '1' in i[0] and 'O=O' not in i[0] and len(i[0])<150:
        # if  '.O'  not in i[0] and '.O.' not in i[0] and len(i[0])<150:
         if  'NH4+'  in i[0] and 'NH3+' not in i[0] and len(i[0])<50:
           reactionsOO_withI.append(i)


for i in  range(len((CGRSSmilesGenerated))):
    for elem,sm in zip(objGen[i],cgrGen[i]):
        #print(sm,str(elem))
        reaction_center=getReactionCentersDictionary([sm],[elem])
        keys=[z for z in reaction_center.keys()]
        if     len(keys)==1 and keys[0] not in rcsAll:
            rcsGen[i].update(reaction_center)
            print(len(rcsGen[i]))
            




if False:
    # for lst  in  cgrGen:
    #    Unique.append([x for x in lst  if x not in CGRs])
    from rdkit.Chem import Draw   
    cntt=0
    cnttAll=0
    for x in objs:
    #  try:
        decomposed = cgr.ReactionContainer.from_cgr(x)
        cnttAll+=1
    #    decomposed.explicify_hydrogens();e
        print(str(decomposed))
        m = AllChem.ReactionFromSmarts(str(decomposed))
        
        if m is not None:
          cntt+=1  
          Draw.ReactionToImage(m)
          print(cntt)
    #    decomposed.clean2d()
        #decomposed
    #  except:
     #     pass
    
    
    
    
    #if notPickled:
        
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
from CGRtools import CGRpreparer #


CGRSmileFileIn="GeneratedData/fwdRNN.dat"
generator=cgr.files.SMILESRead(CGRSmileFileIn)





def isValidCGR(nextReactionObj):
    if type(nextReactionObj)!=cgr.containers.cgr.CGRContainer:
        return False
    else:
        return True



preparer = CGRpreparer()






with open(CGRSmileFileIn, 'r') as the_file:
#   try:
        CGRs=[]
        for line in the_file:
            line=line.rstrip('\n')
            CGRs.append(line)
Cnt_validCGR=0


while True:
   try:
    nextCGRObj=next(generator)
    if isValidCGR(nextCGRObj):
        reactionObj=preparer.decompose(nextCGRObj)
        print(type(reactionObj))
        

        
        #        print(type(nextReactionObj))
#        nextReactionObj.standardize()

    # isStn=nextReactionObj.standardize()
#        cgrContainer=~nextReactionObj

#    else:
#        nextReactionObj.standardize()
#        cgrContainer=~nextReactionObj
        print("")
        reactant_part, product_part = ~nextCGRObj
        print("Reactant_part:",reactant_part)
        print("product_part:",product_part)

#        cgrContainer=nextReactionObj



#    nextReactionStr=str(cgrContainer)
    print("CGR:",str(nextReactionObj))

#            print(line)
#  while True:
#   try:
#    nextReactionObj=next(generator)
#    if type(nextReactionObj)!=cgr.containers.cgr.CGRContainer:
#        nextReactionObj.standardize()
        
#    # isStn=nextReactionObj.standardize()
#        cgrContainer=~nextReactionObj
#    else:
#        cgrContainer=nextReactionObj
        
   #aa,bb,dd=GetAllfromRDkit(str(nextReactionObj),i)
   #print(i,aa,bb,dd)
   # if str(cgrContainer)  not in ll1:
   #     ll1.append(str(cgrContainer))
   #     if '>-^>' in str(cgrContainer):
   #         print(cgrContainer)
       
   
    #nextReactionStr=str(nextReactionObj)
#    nextReactionStr=str(cgrContainer)
#    if len(nextReactionStr)>156:
#        continue
#    length.append(len(nextReactionStr))
#    print(nextReactionStr)
#    the_file.write(nextReactionStr+'\n')
   
   
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
   
   #break

   
   except Exception as e:
       print("########### SOME ERROR ##################")
#       print()
       
       print(e)
       
       
       print("Error")



#print(len(CGRs))


#print("mean of Seq. lenght:", np.mean(length))
#print("min of Seq. lenght:", np.min(length))
#print("max of Seq. lenght:", np.max(length))
   
   
 #    pass
#Chem.SmilesMolSupplier

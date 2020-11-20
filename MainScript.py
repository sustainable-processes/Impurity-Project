# -*- coding: utf-8 -*-
"""
Created on Fri Nov 20 11:00:43 2020

@author: AADA01
"""
from PredictImpurities import predict_impurities,visualizeoutput
import itertools
from collections import Counter
from MainFunctions import getlist

#%%
rxnlib,smles,analogue_rxns,template_dict,error_dict=predict_impurities('Case2',screeningdone=True)
visualizeoutput('Case2', rxnlib, analogue_rxns, template_dict)

#%%
impuritylist=getlist(template_dict,'Impurity Smiles')
count_impurity=Counter(map(tuple,impuritylist))

# for rxnid,rxn in template_dict.items():
#     if rxn.get('Impurity Smiles'):
#         if rxn['Impurity Smiles'] !=[]:
#             if 'CC(=O)Oc1ccc(N)cc1' in [spec for speclist in rxn['Impurity Smiles'] for spec in speclist]:
#                 print(rxnid)
#%%
error_count=Counter(getlist(error_dict,'Reason for screen out'))
rxnsearch=[rxnid for rxnid in error_dict if error_dict[rxnid]['Reason for screen out']=='Likely mismatch in number of template and query reactants. Check template.']

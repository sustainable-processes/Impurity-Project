# -*- coding: utf-8 -*-
"""
Created on Fri Oct  9 11:24:30 2020

@author: AADA01
"""

from MainFunctions import molfromsmiles,openpickle, getfragments,maprxn,rdMolDraw2D
from DataProcessing import rxnlib,idlist




#%% DATA PROCESSING
# Data is in pickle format (need to deserialize). No need pandas or numpy, can 
# directly read. This cell builds a substance dictionary with key as reactant/product Reaxys
# ID and values being the smiles string (supplied already as Smiles.pkl),
 # mol file, mole drawing, and molecular formula. JUST RUN ONCE, UNLESS SMILES.PKL
 # changes.

# smles=openpickle('smiles') #Reference dictionary with keys as rxs id and values as smiles strings
# for rxsid,smiles in smles.items():
#     prod,fig=molfromsmiles(smiles)
#     mol_formula=Chem.rdMolDescriptors.CalcMolFormula(prod)
#     smles[rxsid]={'Smiles': smiles, 'Mol':prod,'Struct':fig, 'Formula':mol_formula}
    
    
    
    
    
    # reac=getfragments(rctids,smles) #Calls getfragments() to generate reaction string containing reactants (smiles strings from reference substance dictionary)
    # prod=getfragments(prdids,smles) ##Calls getfragments() to generate reaction string containing products (smiles strings from reference substance dictionary)
    # if rgtids[0]!='': #Code doesn' work yet since substance dictionary doesn' have solvents, catalysts and reagents
    #     rgt=getfragments(rgtids, smles)
    # else:
    #     rgt=''
    # rxnstr='{}>>{}'.format(reac, prod) #Assembling final reaction string

# Mapping reactions using IBM Rxn mapper

# results4=maprxn(rxn_fin)

#Extracting reaction center
candidate_rxns=[]
screened_rxns=[]
rxn_smarts=[]


# for maps in results4:
#     curr_rxn=rdChemReactions.ReactionFromSmarts(maps['mapped_rxn'],useSmiles=True)
#     curr_rxn.RemoveUnmappedReactantTemplates()
#     candidate_rxns.append(curr_rxn)
#     rxn_smarts.append(drawReaction(curr_rxn))
#     res,rmap=parsemap(curr_rxn)
#     changed_atoms, changed_mapidx=get_changed_atoms(mapdict) #Works now
    


# directory=os.path.join(os.getcwd(),'Images')

# for rxnimg in rxn_smarts:
#     writetofile(rxnimg)
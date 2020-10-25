# -*- coding: utf-8 -*-
"""
Created on Fri Oct  9 10:39:17 2020

This script processes candidate reactions from the upstream workflow into a dictionary. Also 
produces an ID list which is used to retrieve SMILES strings and mole files from the server. 

@author: AADA01
"""

from MainFunctions import Chem, rdChemReactions,molfromsmiles,openpickle, getfragments,gethelpfragments,maprxn,rdMolDraw2D,drawReaction,parsemap,get_changed_atoms, os, writetofile,getlist,convSVGtoPNG, balance_stoichiometry,writepickle,json,isbalanced,hc_smilesDict,hc_molDict
import copy

#%% CANDIDATE REACTION PROCESSING (0 to 14; refer to excel sheet)
#  Data is in pickle format (need to deserialize). No need pandas or numpy, can 
# directly read.This cell assembles a dictionary based on candidate reactions extracted from
# Reaxys. Requires file candidate_rxns.pkl in working directory.
 # ONLY RUN ONCE IF candidate_rxns.pkl changes

candidate_rxns=openpickle('candidate_rxns.pkl')
rxn_list=candidate_rxns['rxns'][1:] #Removing header (Contains formatting)
rxt_assigns=candidate_rxns['rxtsAssigns']

delim='\t'  #Tab delimiter based on structure of candidate_rxns (change if different)

rxnlib={}
idlist=[]

for rxninfo in rxn_list:
    splitstring=rxninfo.split(delim)
    reacid=splitstring[0]            #Reaxys ID of reaction
    rctids=splitstring[1].split(",") #Reaxys ID of reactants, split based on ,
    if rctids[0]!='':
        idlist+=rctids
    prdids=splitstring[2].split(",") #Reaxys ID of products, split based on ,
    if prdids[0]!='':
        idlist+=prdids
    numref=splitstring[3]            #Number of references to the reaction (literature)
    yld=splitstring[4].split(";")    #Yield of reaction products (yield (%): reaxys ID)
    temp=splitstring[5]              #Temperature of reaction
    press=splitstring[6]             #Pressure of reaction
    pH=splitstring[7]                #pH of reaction
    restime=splitstring[8]           # Residence time of reaction
    numstep=splitstring[9]           # Number of steps of reaction 
    rgtids=splitstring[11].split(",") #Reaxys ID of reagents, split based on ,
    if rgtids[0]!='':
        idlist+=rgtids
    solvids=splitstring[10].split(",") #Reaxys ID of solvents, split based on ,
    if solvids[0]!='':
        idlist+=solvids
    catids=splitstring[12].split(",")#Reaxys ID of catalysts, split based on ,
    if catids[0]!='':
        idlist+=catids
    year=splitstring[14][:-1]        #Year published
    rxnlib[reacid]={'Reactants': rctids, 'Products': prdids, 'Reagents': rgtids, 'Solvent': solvids, 'Catalyst': catids, 'Temperature': temp,'Pressure': press,'pH': pH,'Residence Time': restime,'Yield': yld,'Steps': numstep,'Num Ref': numref,'Year': year}


analogue_compds=openpickle('analogue_compounds.pkl')



#%% Substance dictionary
# Data is in pickle format (need to deserialize). No need pandas or numpy, can 
# directly read. This cell builds a substance dictionary with key as reactant/product Reaxys
# ID and values being the smiles string (supplied already as Smiles.pkl, required in working directory),
 # mol file, mole drawing, and molecular formula. JUST RUN ONCE, UNLESS SMILES.PKL
 # changes.

smles=openpickle('smiles.pkl') #Reference dictionary with keys as rxs id and values as smiles strings
#%% Temp code to generate update smiles of solvents, reagents, catalysts [automate when given access to database]
boollist=[i in list(smles.keys()) for i in idlist]
check=[idlist[i] for i,val in enumerate(boollist) if not boollist[i]]
smles['11323289']='O=[Mn]=O' #Wrong in Reaxys (O.O.[Mn])
smles['635760']='CC1=CC=CC=C1'
smles['1098229']='CO'
smles['1098293']='S=C=S'
smles['11342937']='[OH-].[K+]'
smles['11343689']='[K+].[O-]S(=O)(=O)OOS([O-])(=O)=O'
smles['13720881']='OS(=O)(=O)O' #May not be in the database
smles['1921286']='OP(O)(O)=O'
smles['21093063']='[Rh].[C]' #May not be in database (catalyst so doesn't matter)
smles['3587189']='[HH]'
smles['3587155']='O'
smles['102551']='C1COCCO1'
smles['472690']='CC1=CC=C(C=C1)S(=O)(=O)O'
smles['605283']='CCN(CC)CC'
smles['1730800']='C(Cl)Cl'
smles['11460447']='[H-].[H-].[H-].[H-].[Li+].[Al+3]'
smles['102391']='C1CCOC1'
smles['1718733']='CCO'
smles['20782673']='[Pd]'
smles['3659978']='[OH-].[OH-].[Pd+2]'
smles['3587194']='II'
smles['1098295']='C(Cl)(Cl)(Cl)Cl'

#Continue susbtance dictionary creation

for rxsid,smiles in smles.items():
    carrier_frag=''
    # prod,fig=molfromsmiles(smiles)
    prod=molfromsmiles(smiles)
    mol_formula=Chem.rdMolDescriptors.CalcMolFormula(prod)
    for dic in rxt_assigns:
        if dic.get(rxsid):
            carrier_frag=dic.get(rxsid)[0][1]
            query_compd=dic.get(rxsid)[0][0]
    # for carrier, analoguelist in analogue_compds.items():
    #     if rxsid in analoguelist:
    #         carrier_frag=carrier
    
    # smles[rxsid]={'Smiles': smiles, 'Mol':prod,'Struct':fig, 'Formula':mol_formula}
    smles[rxsid]={'Smiles': smiles, 'Mol':prod, 'Formula':mol_formula}
    if carrier_frag:
        smles[rxsid].update({'Carrier Fragment': carrier_frag,'Query Compound': query_compd}) #Updating carrier fragment to which analogue compound belongs. Reagents/solvents don't have this

#%% Building reaction strings and mapping (WITHOUT SCREENING)

for rxnid,rxn in rxnlib.items():
    reacstrl=[] #Empty list of chemicals on reaction LHS
    reac=getfragments(rxn['Reactants'],smles) #Calls getfragments() to generate reaction string containing reactants (smiles strings from reference substance dictionary)
    reacstrl.append(reac)
    if rxn['Reagents'][0]!='': 
        reag=getfragments(rxn['Reagents'], smles)
        reacstrl.append(reag)
    if rxn['Solvent'][0]!='':
        solv=getfragments(rxn['Solvent'],smles)
        reacstrl.append(solv)
    reacstr='.'.join(reacstrl)
    prodstr=getfragments(rxn['Products'],smles)  #Calls getfragments() to generate reaction string containing products (smiles strings from reference substance dictionary)
    currrxnstr='{}>>{}'.format(reacstr,prodstr) #reacstr if want to include reagents, solvents 
    currrxnstr2='{}>>{}'.format(reac,prodstr) #reacstr2 contained just reactants
    rxn.update({'RSmiles': currrxnstr}) #Updating rxnlib dictionary with smarts string
    rxn.update({'RSmiles2': currrxnstr2}) #Updating rxnlib dictionary with smarts string (only reactants)


# Mapping reactions using IBM Rxn mapper

#Reactants, reagents, solvents
rxnstr=getlist(rxnlib,'RSmiles')
results4=maprxn(rxnstr)

# Reactants only

rxnstr2=getlist(rxnlib,'RSmiles2')
results5=maprxn(rxnstr2)

#Visualizing reactions and extracting reaction center

#Reactants, reagents, solvents

for maps,rxn in zip(results4,list(rxnlib.values())):
    curr_rxn=rdChemReactions.ReactionFromSmarts(maps['mapped_rxn'],useSmiles=True)
    curr_rxn.RemoveUnmappedReactantTemplates()
    curr_rxndraw=drawReaction(curr_rxn)
    rxn.update({'Mapping': maps['mapped_rxn'], 'Confidence': maps['confidence'], 'RDKit Rxn': curr_rxn, 'Sketch': curr_rxndraw})
    res,rmap=parsemap(curr_rxn)
    changed_atoms, changed_mapidx=get_changed_atoms(res) #Works now
    rxn.update({'Reaction Center': changed_mapidx})


# Reactants only

for maps,rxn in zip(results5,list(rxnlib.values())):
    curr_rxn2=rdChemReactions.ReactionFromSmarts(maps['mapped_rxn'],useSmiles=True)
    curr_rxn2.RemoveUnmappedReactantTemplates()
    curr_rxndraw2=drawReaction(curr_rxn2)
    rxn.update({'Mapping2': maps['mapped_rxn'], 'Confidence2': maps['confidence'], 'RDKit Rxn2': curr_rxn2, 'Sketch2': curr_rxndraw2})
    res2,rmap2=parsemap(curr_rxn2)
    changed_atoms2, changed_mapidx2=get_changed_atoms(res2) #Works now
    rxn.update({'Reaction Center2': changed_mapidx2})

#%%Building reaction strings and mapping (WITH SCREENING). RUN ONCE

analogue_rxns={} #To store only screened reaction data with updated reaction string and data
for rxnid,rxn in rxnlib.items():
    balance=isbalanced(rxnid,rxnlib,smles)
    if balance: #Only reactions where isbalance returns something
        analogue_rxns.update({rxnid:rxn}) 
        if type(balance)==tuple and balance[1]!='Warning': #This means reaction has been balanced using chempy, returns tuple
            reacst=balance[1] # Reactants and stoich coefficients
            prodst=balance[2] # Products and stoich coefficients
            reaclist=[] #Empty list of chemicals on reaction LHS
            rhelplist=[] #Empty list of help compounds on reaction LHS
            prodlist=[] #Empty list of chemicals on reaction RHS
            phelplist=[] #Empty list of help compounds on reaction RHS 
            
            for reac,coeff in reacst.items():
                helpcomp=True
                for reactant in rxn['Reactants']:
                    if smles[reactant]['Formula']==reac:
                        reaclist.extend([reactant for _ in range(coeff)])
                        helpcomp=False
               
                if helpcomp:
                    rhelplist.extend([reac for _ in range(coeff)])
                    
            reacstr=getfragments(reaclist,smles) 
            if rhelplist:
                reacstr=reacstr+'.'+gethelpfragments(rhelplist,hc_smilesDict)
                rxn.update({'Help Reactants': rhelplist})
        
            for prod,coeff in prodst.items():
                helpcomp=True
                for product in rxn['Products']:
                    if smles[product]['Formula']==prod:
                        prodlist.extend([product for _ in range(coeff)])
                        helpcomp=False
                if helpcomp:
                    phelplist.extend([prod for _ in range(coeff)])
                    
            prodstr=getfragments(prodlist,smles) 
            if phelplist:
                prodstr=prodstr+'.'+gethelpfragments(phelplist,hc_smilesDict)
                rxn.update({'Help Products': phelplist})
            currrxnstr3='{}>>{}'.format(reacstr,prodstr)
            rxn.update({'Balanced RSmiles': currrxnstr3})
            
# Mapping
balance_ids=[key for key,rxn in analogue_rxns.items() if 'Balanced RSmiles' in rxn.keys()]
rxnstr3=[analogue_rxns[key]['Balanced RSmiles'] for key in balance_ids]
results6=maprxn(rxnstr3)
                
for maps,rxnid in zip(results6,balance_ids):
    curr_rxn3=rdChemReactions.ReactionFromSmarts(maps['mapped_rxn'],useSmiles=True)
    curr_rxn3.RemoveUnmappedReactantTemplates()
    curr_rxndraw3=drawReaction(curr_rxn3)
    analogue_rxns[rxnid].update({'Balanced Mapping': maps['mapped_rxn'], 'Balanced Confidence': maps['confidence'], 'Balanced RDKit Rxn': curr_rxn3, 'Balanced Sketch': curr_rxndraw3})
    res3,rmap3=parsemap(curr_rxn3)
    changed_atoms3, changed_mapidx3=get_changed_atoms(res3) #Works now
    analogue_rxns[rxnid].update({'Balanced Reaction Center': changed_mapidx3})
    
#%% Further screening based on reaction center
for rxnid,rxn in analogue_rxns.items():
    if 'Balanced Reaction Center' in rxn.keys():
        RC=rxn['Balanced Reaction Center']
        rdrxn=rxn['Balanced RDKit Rxn']
    else:
        RC=rxn['Reaction Center2']
        rdrxn=rxn['RDKit Rxn2']
    clean_rxn=copy.copy(rdrxn)
    rdChemReactions.RemoveMappingNumbersFromReactions(clean_rxn)
    mapset=set()
    for idx,reactant in enumerate(clean_rxn.GetReactants()):
        rdreactant=rdrxn.GetReactants()[idx]
        for reacid in rxn['Reactants']:
            if smles[reacid]['Smiles']==Chem.MolToSmiles(reactant):
                carrier_frag=smles[reacid]['Carrier Fragment']
                break
        matches=rdreactant.GetSubstructMatches(Chem.RemoveAllHs(Chem.MolFromSmarts(carrier_frag))) #tuples of atom indices
        if mapset:
            mapset=mapset.union({rdreactant.GetAtomWithIdx(atomidx).GetAtomMapNum() for match in matches for atomidx in match})
        else:
            mapset={rdreactant.GetAtomWithIdx(atomidx).GetAtomMapNum() for match in matches for atomidx in match}
    if set(RC).issubset(mapset) and RC:
        template_dict={rxnid: rxn}
        
        
        
        
        
                
            
            
            
            
 #%% Write libraries to pickle file

# writepickle(rxnlib,'rxnlib')
# writepickle(smles,'smles')

#%% Write libraries to json file
# with open('rxnlib.json', 'w') as handle:
#     json.dump(rxnlib, handle)
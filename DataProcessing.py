# -*- coding: utf-8 -*-
"""
Created on Fri Oct  9 10:39:17 2020

This file processes candidate reactions from the upstream workflow into a dictionary. Also 
produces an ID list which is used to retrieve SMILES strings and mole files from the server, building 
a substance dictionary, smles. This code should only be run on the server. IMPORTANT: This code
needs candidate_rxns.pkl (generated from upstream workflow) in the working directory\\Input\\casenum folder

@author: AADA01
"""

from MainFunctions import getMols, Chem, rdChemReactions,molfromsmiles,openpickle, getfragments,gethelpfragments,maprxn,rdMolDraw2D,drawReaction,parsemap,get_changed_atoms, os, writetofile,getlist,convSVGtoPNG, balance_stoichiometry,writepickle,json,isbalanced,hc_smilesDict,hc_molDict, valid_rxn_center
from FindFunctionalGroups import identify_functional_groups as IFG

def processdata(casenum):
    sep=os.sep
    inputdir=os.path.join(os.getcwd(),'Input'+sep+casenum) 
    #%% CANDIDATE REACTION PROCESSING (0 to 14; refer to excel sheet)
    #  Data is in pickle format (need to deserialize). No need pandas or numpy, can 
    # directly read.This cell assembles a dictionary based on candidate reactions extracted from
    # Reaxys. Requires file candidate_rxns.pkl in input directory of case.
    
    
    candidate_rxns=openpickle(os.path.join(inputdir,'candidate_rxns.pkl'))
    rxn_list=candidate_rxns['rxns'][1:] #Removing header (Contains formatting)
    rxt_assigns=candidate_rxns['rxtsAssigns']
    # analogue_compds=openpickle(os.path.join(inputdir,'analogue_compounds.pkl'))
    
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
    
    #%% Substance dictionary
    # Data is in pickle format (need to deserialize). No need pandas or numpy, can 
    # directly read. This cell builds a substance dictionary with key as reactant/product Reaxys
    # ID and values being the smiles string (supplied already as Smiles.pkl, required in working directory),
     # mol file, mole drawing, and molecular formula. JUST RUN ONCE, UNLESS SMILES.PKL
     # changes.
    
    # if os.path.isdir(os.path.join(inputdir,'smiles.pkl')): #If a smiles list is already provided
    #     smles=openpickle(os.path.join(inputdir,'smiles.pkl')) #Reference dictionary with keys as rxs id and values as smiles strings IF IT EXISTS
    #     boollist=[i in list(smles.keys()) for i in idlist]
    #     check=[idlist[i] for i,val in enumerate(boollist) if not boollist[i]]
    #     
    # In case code doesn't work for case 1
    # smles['11323289']='O=[Mn]=O' #Wrong in Reaxys (O.O.[Mn])
    # smles['635760']='CC1=CC=CC=C1'
    # smles['1098229']='CO'
    # smles['1098293']='S=C=S'
    # smles['11342937']='[OH-].[K+]'
    # smles['11343689']='[K+].[O-]S(=O)(=O)OOS([O-])(=O)=O'
    # smles['13720881']='OS(=O)(=O)O' #May not be in the database
    # smles['1921286']='OP(O)(O)=O'
    # smles['21093063']='[Rh].[C]' #May not be in database (catalyst so doesn't matter)
    # smles['3587189']='[HH]'
    # smles['3587155']='O'
    # smles['102551']='C1COCCO1'
    # smles['472690']='CC1=CC=C(C=C1)S(=O)(=O)O'
    # smles['605283']='CCN(CC)CC'
    # smles['1730800']='C(Cl)Cl'
    # smles['11460447']='[H-].[H-].[H-].[H-].[Li+].[Al+3]'
    # smles['102391']='C1CCOC1'
    # smles['1718733']='CCO'
    # smles['20782673']='[Pd]'
    # smles['3659978']='[OH-].[OH-].[Pd+2]'
    # smles['3587194']='II'
    # smles['1098295']='C(Cl)(Cl)(Cl)Cl'
    
    smles={}
    mols=getMols(idlist)
    for idx,mol in enumerate(mols):
        if mol:
            subst_id=idlist[idx]
            for dic in rxt_assigns:
                if dic.get(subst_id):
                    carrier_frag=dic.get(subst_id)[0][1]
                    query_compd=dic.get(subst_id)[0][0]
            subst_info={'Smiles': Chem.MolToSmiles(mol),'Mol': mol,'Formula': Chem.rdMolDescriptors.CalcMolFormula(mol)}
            subst_info.update({'Carrier Fragment': carrier_frag,'Query Compound': query_compd})
            smles.update({subst_id: subst_info})
        else:
            continue
    writepickle(rxnlib,os.path.join(inputdir,'rxnlib'))
    writepickle(smles,os.path.join(inputdir,'smles'))
        
        
        

    
    # for rxsid,smiles in smles.items():
    #     carrier_frag=''
    #     prod=molfromsmiles(smiles)
    #     mol_formula=Chem.rdMolDescriptors.CalcMolFormula(prod)
    #     for dic in rxt_assigns:
    #         if dic.get(rxsid):
    #             carrier_frag=dic.get(rxsid)[0][1]
    #             query_compd=dic.get(rxsid)[0][0]
    #     smles[rxsid]={'Smiles': smiles, 'Mol':prod, 'Formula':mol_formula}
    #     if carrier_frag:
    #         smles[rxsid].update({'Carrier Fragment': carrier_frag,'Query Compound': query_compd}) #Updating carrier fragment to which analogue compound belongs. Reagents/solvents don't have this
    
    
    

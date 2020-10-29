# -*- coding: utf-8 -*-
"""
Created on Thu Oct 29 17:16:58 2020

Irrelevant reactions (requiring non-analogue compounds,unspecified/unknown 
chemicals to balance reaction, or react somewhere other than the functional group in the carrier fragment)
are filtered out. Template_dict is the final output which contains balanced, relevant reactions for
template generation. IMPORTANT: If DataProcessing has never been run before, RUN IT FIRST on the server.
It will populate the input folder with the substance dictionary and adapt the candidate reaction dictionary,
both of which are needed for screening.


@author: E0014
"""

from MainFunctions import getMols, Chem, rdChemReactions,molfromsmiles,openpickle, getfragments,gethelpfragments,maprxn,rdMolDraw2D,drawReaction,parsemap,get_changed_atoms, os, writetofile,getlist,convSVGtoPNG, balance_stoichiometry,writepickle,json,isbalanced,hc_smilesDict,hc_molDict, valid_rxn_center
from FindFunctionalGroups import identify_functional_groups as IFG
import copy


def screening(casenum):
    
    inputdir=os.path.join(os.getcwd(),'Input\\'+casenum)
    rxnlib=openpickle(os.path.join(inputdir,'rxnlib.pickle'))
    smles=openpickle(os.path.join(inputdir,'smles.pickle'))

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
        rxn.update({'Mapping Dictionary': res,'Reaction Center': changed_mapidx})
    
    
    # Reactants only
    
    for maps,rxn in zip(results5,list(rxnlib.values())):
        curr_rxn2=rdChemReactions.ReactionFromSmarts(maps['mapped_rxn'],useSmiles=True)
        curr_rxn2.RemoveUnmappedReactantTemplates()
        curr_rxndraw2=drawReaction(curr_rxn2)
        rxn.update({'Mapping2': maps['mapped_rxn'], 'Confidence2': maps['confidence'], 'RDKit Rxn2': curr_rxn2, 'Sketch2': curr_rxndraw2})
        res2,rmap2=parsemap(curr_rxn2)
        changed_atoms2, changed_mapidx2=get_changed_atoms(res2) #Works now
        rxn.update({'Mapping Dictionary2': res2,'Reaction Center2': changed_mapidx2})
    
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
                rxn.update({'Balanced Reactants': reaclist, 'Balanced Products': prodlist})
                
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
        analogue_rxns[rxnid].update({'Balanced Mapping Dictionary': res3,'Balanced Reaction Center': changed_mapidx3})
        
    #%% Further screening based on reaction center
    template_dict={}
    for rxnid,rxn in analogue_rxns.items():
        valid=valid_rxn_center(rxnid,analogue_rxns,smles)
        if valid:
            rxn.update({'Fragment Location': valid[1]})
            rxn.update({'Clean Reaction': valid[2]})
            template_dict.update({rxnid:rxn})
    
    return rxnlib,smles,analogue_rxns,template_dict
        
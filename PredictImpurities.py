# -*- coding: utf-8 -*-
"""
Created on Fri Oct  9 11:24:30 2020

This is the main script, and uses outputs from DataProcessing.py as well as
main functions defined in MainFunctions.py. Extracts relevant templates
and generates impurities

@author: AADA01
"""

from MainFunctions import Chem, rdChemReactions,molfromsmiles,openpickle, getfragments,maprxn,rdMolDraw2D,drawReaction,parsemap,get_changed_atoms, os, writetofile,getlist,convSVGtoPNG, balance_stoichiometry,writepickle,copy,gen_template

# #%% Importing libraries directly from Screening.py
from Screening import screening


#%%  OR Reading libraries from file
# directory=os.path.join(os.getcwd(),'Dictionaries\\Case1')
# rxnlib=openpickle(os.path.join(directory,'rxnlib.pickle'))
# smles=openpickle(os.path.join(directory,'smles.pickle'))
# analogue_rxns=openpickle(os.path.join(directory,'analogue_rxns.pickle'))
# template_dict=openpickle(os.path.join(directory,'template_dict.pickle'))

def predict_impurities(casenum, visualizeoutput=True):

    rxnlib,smles,analogue_rxns,template_dict=screening(casenum)
    
    #%% Extracting template
    for rxnid,rxn in template_dict.items():
        fragloc=rxn['Fragment Location']
        clean_rxn=rxn['Clean Reaction']
        if 'Balanced Reaction Center' in rxn.keys():
            RC=rxn['Balanced Reaction Center']
            rdrxn=rxn['Balanced RDKit Rxn']
            mapdict=rxn['Balanced Mapping Dictionary']
            rctlist=rxn['Balanced Reactants']
        else:
            RC=rxn['Reaction Center2']
            rdrxn=rxn['RDKit Rxn2']
            mapdict=rxn['Mapping Dictionary2']
            rctlist=rxn['Reactants']
        
        oldreacid='' #To keep track of repeat reactants
        reacfrag='' #Store smart strings of fragments
        reacmap=set() #Cumulatively store mapping indices of fragments
        reacidx=set() #Cumulatively store atom indices of fragments
        for idx,rdreactant in enumerate(rdrxn.GetReactants()):
            recur=0
            clean_reactant=clean_rxn.GetReactants()[idx]
            for reacid in rctlist:
                if smles[reacid]['Smiles']==Chem.MolToSmiles(clean_reactant): #Redundant...just to cross-check
                    break
            if not oldreacid:
                oldreacid=reacid
            elif reacid==oldreacid:
                recur+=1
            else:
                oldreacid=reacid
                
            fragidx=fragloc[reacid][1][recur] #Local (reactant) store of atom indices of fragments
            fragmap=fragloc[reacid][0][recur] #Local (reactant) store of mapping indices of fragments
            if reacfrag:
                reacfrag=reacfrag+'.'+Chem.MolFragmentToSmarts(rdreactant,fragidx)
                reacmap=reacmap.union(fragmap) #Probably can retrieve directly from template_dict
                reacidx=reacidx.union(fragidx) #Probably can retrieve directly from template_dict
            else:
                reacfrag=Chem.MolFragmentToSmarts(rdreactant,fragidx)
                reacmap=fragmap
                reacidx=fragidx
    
        prodfrag=''
        reacmapiter=copy.copy(reacmap)
        for idx,rdproduct in enumerate(rdrxn.GetProducts()):
            prodidx=set()
            for mapnum in reacmap:
                if Chem.MolToSmarts(mapdict[mapnum][0])==Chem.MolToSmarts(rdproduct) and mapnum in reacmapiter:
                    prodidx.add(mapdict[mapnum][2])
                    reacmapiter.remove(mapnum)
            if prodfrag:
                prodfrag=prodfrag+'.'+Chem.MolFragmentToSmarts(rdproduct,prodidx)
            else:
                prodfrag=Chem.MolFragmentToSmarts(rdproduct,prodidx)
            if not reacmapiter:
                break
    
        template=gen_template(reacfrag,prodfrag)
        template_rxn=rdChemReactions.ReactionFromSmarts(template, useSmiles=True)
        template_sketch=drawReaction(template_rxn)
        rxn.update({'Template': template, 'Template Reaction': template_rxn,'Template Sketch': template_sketch})
        rxn.update({'Query Reactants':[smles[reacid]['Query Compound'] for reacid in rctlist]})
    
    #%% Apply template to generate impurities
    for rxnid,rxn in template_dict.items():
        template_rxn=rxn['Template Reaction']
        query_reactants=rxn['Query Reactants']
        imp_raw=template_rxn.RunReactants([Chem.MolFromSmiles(query_reactant) for query_reactant in query_reactants]) #Raw output of code
        imp_smles=[tuple(Chem.MolToSmiles(imp) for imp in imp_prod) for imp_prod in imp_raw]
        imp_smles=list(set(tuple(sorted(t)) for t in imp_smles))
        imp_mol=[tuple(Chem.MolFromSmiles(impsmles) for impsmles in tup) for tup in imp_smles]
        rxn.update({'Impurity Smiles': imp_smles, 'Impurity Molecules': imp_mol})
        
        #Generating final impurity reaction
        reacfrag='.'.join(query_reactants)
        imp_RSmiles=[gen_template(reacfrag,'.'.join(tup)) for tup in imp_smles]
        imp_RDKit_rxn=[rdChemReactions.ReactionFromSmarts(imp_rxn,useSmiles=True) for imp_rxn in imp_RSmiles]
        imp_sketch=[drawReaction(rxn) for rxn in imp_RDKit_rxn]
        rxn.update({'Impurity Reaction Smiles': imp_RSmiles, 'Impurity RDKit Reaction': imp_RDKit_rxn, 'Impurity Reaction Sketch': imp_sketch})


     #%% Write libraries to pickle file
    direc=os.path.join(os.getcwd(),'Output_Dictionaries\\'+casenum)
    if not os.path.isdir(direc):
        os.makedirs(direc)
    writepickle(rxnlib,os.path.join(direc,'rxnlib'))
    writepickle(smles,os.path.join(direc,'smles'))
    writepickle(analogue_rxns,os.path.join(direc,'analogue_rxns'))
    writepickle(template_dict,os.path.join(direc,'template_dict'))
    
    
    
    #%% Writing images to file (remove comment if case changes)
    if visualizeoutput:
        sep=os.sep
        directory=os.path.join(os.getcwd(),'Images'+sep+casenum+sep+'Sketches_Full')
        if not os.path.isdir(directory):
            os.makedirs(directory)
        for rxnid in rxnlib.keys():
            writetofile(rxnlib[rxnid]['Sketch'],os.path.join(directory,rxnid+'.svg'))
            convSVGtoPNG(os.path.join(directory,rxnid),os.path.join(directory,rxnid))
        
        
        directory=os.path.join(os.getcwd(),'Images'+sep+casenum+sep+'Sketches_ReactantOnly')
        if not os.path.isdir(directory):
            os.makedirs(directory)
        for rxnid in rxnlib.keys():
            writetofile(rxnlib[rxnid]['Sketch2'],os.path.join(directory,rxnid+'.svg'))
            convSVGtoPNG(os.path.join(directory,rxnid),os.path.join(directory,rxnid))
        
        directory=os.path.join(os.getcwd(),'Images'+sep+casenum+sep+'Sketches_Analogue')
        if not os.path.isdir(directory):
            os.makedirs(directory)
        for rxnid in analogue_rxns.keys():
            if 'Balanced Sketch' in analogue_rxns[rxnid].keys():
                writetofile(analogue_rxns[rxnid]['Balanced Sketch'],os.path.join(directory,rxnid+'.svg'))
            else:
                writetofile(analogue_rxns[rxnid]['Sketch'],os.path.join(directory,rxnid+'.svg'))
            convSVGtoPNG(os.path.join(directory,rxnid),os.path.join(directory,rxnid))
        
        directory=os.path.join(os.getcwd(),'Images'+sep+casenum+sep+'Impurities')
        if not os.path.isdir(directory):
            os.makedirs(directory)
        for rxnid in template_dict.keys():
            directory=os.path.join(directory,rxnid)
            if not os.path.isdir(directory):
                os.makedirs(directory)
            writetofile(template_dict[rxnid]['Template Sketch'],os.path.join(directory,'Template.svg'))
            convSVGtoPNG(os.path.join(directory,'Template'),os.path.join(directory,'Template'))
            for idx,sketch in enumerate(template_dict[rxnid]['Impurity Reaction Sketch']):
                writetofile(sketch,os.path.join(directory,'Impurity'+str(idx)+'.svg'))
                convSVGtoPNG(os.path.join(directory,'Impurity'+str(idx)),os.path.join(directory,'Impurity'+str(idx)))   
        

    #%% Write libraries to json file
    # with open('rxnlib.json', 'w') as handle:
    #     json.dump(rxnlib, handle)
    return rxnlib,smles,analogue_rxns,template_dict

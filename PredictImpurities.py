# -*- coding: utf-8 -*-
"""
Created on Fri Oct  9 11:24:30 2020

This is the main script, and uses outputs from DataProcessing.py as well as
main functions defined in MainFunctions.py. Extracts relevant templates
and generates impurities

@author: AADA01
"""
#%%
from MainFunctions import Chem, rdChemReactions,molfromsmiles,openpickle, getfragments, \
maprxn,rdMolDraw2D,drawReaction,parsemap,get_changed_atoms, os, writetofile, \
delcontents,getlist, convSVGtoPNG, balance_stoichiometry,writepickle,copy,gen_template

import itertools
from Screening import screening
from DataProcessing import processdata
from collections import Counter

#%%
def predict_impurities(casenum,datapreprocessingdone=True, screeningdone=False,predictiondone=False):
    if not datapreprocessingdone: #ONLY RUN THIS ON SERVER AND SPECIFY DATAPREPROCESSINGDONE AS FALSE
        processdata(casenum)
    sep=os.sep   
    direc=os.path.join(os.getcwd(),'Output_Screening'+sep+casenum)

    #%% Writing to screening folder
    if not screeningdone:
        rxnlib,smles,analogue_rxns,template_dict,error_dict=screening(casenum)
        if not os.path.isdir(direc):
            os.makedirs(direc)
        if rxnlib:
            writepickle(rxnlib,os.path.join(direc,'rxnlib'))
        if smles:
            writepickle(smles,os.path.join(direc,'smles'))
        if analogue_rxns:
            writepickle(analogue_rxns,os.path.join(direc,'analogue_rxns'))
        if template_dict:
            writepickle(template_dict,os.path.join(direc,'template_dict'))
        if error_dict:
            writepickle(error_dict,os.path.join(direc,'error_dict'))

    #%% Reading from screening folder
    elif not predictiondone:
        rxnlib=openpickle(os.path.join(direc,'rxnlib.pickle'))
        smles=openpickle(os.path.join(direc,'smles.pickle'))
        analogue_rxns=openpickle(os.path.join(direc,'analogue_rxns.pickle'))
        template_dict=openpickle(os.path.join(direc,'template_dict.pickle'))
        error_dict=openpickle(os.path.join(direc,'error_dict.pickle'))
    else:
        direc=os.path.join(os.getcwd(),'Output_Final'+sep+casenum)
        rxnlib=openpickle(os.path.join(direc,'rxnlib.pickle'))
        smles=openpickle(os.path.join(direc,'smles.pickle'))
        analogue_rxns=openpickle(os.path.join(direc,'analogue_rxns.pickle'))
        template_dict=openpickle(os.path.join(direc,'template_dict.pickle'))
        error_dict=openpickle(os.path.join(direc,'error_dict.pickle'))
        return rxnlib,smles,analogue_rxns,template_dict,error_dict         

    #%% Extracting template
    for rxnid,rxn in template_dict.items():
        print(rxnid)
        reacfragloc=rxn['Fragment Location']
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
            rctlist=[rspecies for rspecies in itertools.chain(rxn['Reactants'],rxn['Reagents']) if rspecies]

        oldreacid='' #To keep track of repeat reactants
        reacfrag='' #Store smart strings of fragments
        reacmap=set() #Cumulatively store mapping indices of fragments
        reacidx=set() #Cumulatively store atom indices of fragments
        for idx,rdreactant in enumerate(rdrxn.GetReactants()):
            recur=0
            clean_reactant=clean_rxn.GetReactants()[idx]
            for reacid in rctlist:
                if smles[reacid]['Smiles']==Chem.MolToSmiles(clean_reactant): #Redundant...just to cross-check
                    if rxn.get('Query Reactants'):
                        rxn['Query Reactants'].extend([smles[reacid]['Query Compound']])
                    else:
                        rxn.update({'Query Reactants':[smles[reacid]['Query Compound']]})
                    break
            if not oldreacid:
                oldreacid=reacid
            elif reacid==oldreacid:
                recur+=1
            else:
                oldreacid=reacid

            fragloc=reacfragloc[reacid][recur] #Reactant store of mapping indices, atom indices of all carrier fragments involved in reaction
            fragidx={idx for carrfrag in fragloc.keys() for match in fragloc[carrfrag][1] for idx in match}
            fragmap={mapnum for carrfrag in fragloc.keys() for match in fragloc[carrfrag][0] for mapnum in match}
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

    # return template_dict

    #%% Apply template to generate impurities
    for rxnid,rxn in template_dict.items():
        print(rxnid)
        template_rxn=rxn['Template Reaction']
        query_reactants=rxn['Query Reactants']
        try:
            imp_raw=template_rxn.RunReactants([Chem.MolFromSmiles(query_reactant) for query_reactant in query_reactants]) #Raw output of code
        except Exception:
            errormsg='Likely mismatch in number of template and query reactants. Check template.'
            if error_dict:
                error_dict.update({rxnid:{'Reason for screen out': errormsg}})
            else:
                error_dict={rxnid:{'Reason for screen out': errormsg}}
            continue
        else:
            if not imp_raw:
                errormsg='Template did not work. Likely incompatibility with query reactants.'
                if error_dict:
                    error_dict.update({rxnid:{'Reason for screen out': errormsg}})
                else:
                    error_dict={rxnid:{'Reason for screen out': errormsg}}
                continue
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


     #%% Write final output libraries to pickle file
    sep=os.sep
    direc=os.path.join(os.getcwd(),'Output_Final'+sep+casenum)
    if not os.path.isdir(direc):
        os.makedirs(direc)
    writepickle(rxnlib,os.path.join(direc,'rxnlib'))
    writepickle(smles,os.path.join(direc,'smles'))
    writepickle(analogue_rxns,os.path.join(direc,'analogue_rxns'))
    writepickle(template_dict,os.path.join(direc,'template_dict'))
    writepickle(template_dict,os.path.join(direc,'error_dict'))
    
    return rxnlib,smles,analogue_rxns,template_dict,error_dict

    # with open('rxnlib.json', 'w') as handle: write to json file
    #     json.dump(rxnlib, handle)


def visualizeoutput(casenum,rxnlib,analogue_rxns,template_dict,delandreplace=True,vizfull=True,vizrr=True,vizanalogue=True,vizimpurities=True):
    sep=os.sep
    # %% Full sketches
    if vizfull:
        directory=os.path.join(os.getcwd(),'Images'+sep+casenum+sep+'Sketches_Full')
        if not os.path.isdir(directory):
            os.makedirs(directory)
        elif delandreplace==True:
            delcontents(directory)

        for rxnid in rxnlib.keys():
            writetofile(rxnlib[rxnid]['Sketch'],os.path.join(directory,rxnid+'.svg'))
            convSVGtoPNG(os.path.join(directory,rxnid),os.path.join(directory,rxnid))

    #%% Reactants and reagents only
    if vizrr:
        directory=os.path.join(os.getcwd(),'Images'+sep+casenum+sep+'Sketches_Reactant_Reagent_Only')
        if not os.path.isdir(directory):
            os.makedirs(directory)
        elif delandreplace==True:
            delcontents(directory)
        for rxnid in rxnlib.keys():
            writetofile(rxnlib[rxnid]['Sketch2'],os.path.join(directory,rxnid+'.svg'))
            convSVGtoPNG(os.path.join(directory,rxnid),os.path.join(directory,rxnid))
    #%% Analogue reactions only
    if vizanalogue:
        directory=os.path.join(os.getcwd(),'Images'+sep+casenum+sep+'Sketches_Analogue')
        if not os.path.isdir(directory):
            os.makedirs(directory)
        elif delandreplace==True:
            delcontents(directory)
        for rxnid in analogue_rxns.keys():
            if 'Balanced Sketch' in analogue_rxns[rxnid].keys():
                writetofile(analogue_rxns[rxnid]['Balanced Sketch'],os.path.join(directory,rxnid+'.svg'))
            else:
                writetofile(analogue_rxns[rxnid]['Sketch2'],os.path.join(directory,rxnid+'.svg'))
            convSVGtoPNG(os.path.join(directory,rxnid),os.path.join(directory,rxnid))
    #%% Impurities
    if vizimpurities:
        directory1=os.path.join(os.getcwd(),'Images'+sep+casenum+sep+'Impurities')
        if not os.path.isdir(directory1):
            os.makedirs(directory1)
        elif delandreplace==True:
            delcontents(directory1)
        for rxnid in template_dict.keys():
            if not template_dict[rxnid].get('Template Sketch'):
                continue
            directory=os.path.join(directory1,rxnid)
            if not os.path.isdir(directory):
                os.makedirs(directory)
            writetofile(template_dict[rxnid]['Template Sketch'],os.path.join(directory,'Template.svg'))
            convSVGtoPNG(os.path.join(directory,'Template'),os.path.join(directory,'Template'))
            if not template_dict[rxnid].get('Impurity Reaction Sketch'):
                continue
            for idx,sketch in enumerate(template_dict[rxnid]['Impurity Reaction Sketch']):
                writetofile(sketch,os.path.join(directory,'Impurity'+str(idx)+'.svg'))
                convSVGtoPNG(os.path.join(directory,'Impurity'+str(idx)),os.path.join(directory,'Impurity'+str(idx)))
        #%% Rejected sketches
        long_ids=[key for key,rxn in rxnlib.items() if rxn.get('Balanced RSmiles') and not rxn.get('Balanced Mapping')]
        if long_ids:
            directory=os.path.join(os.getcwd(),'Images'+sep+casenum+sep+'Sketches_Rejected')
            if not os.path.isdir(directory):
                os.makedirs(directory)
            elif delandreplace==True:
                delcontents(directory)
            for long_id in long_ids:
                writetofile(drawReaction(rdChemReactions.ReactionFromSmarts(rxnlib[long_id]['Balanced RSmiles'],useSmiles=True)),os.path.join(directory,long_id+'.svg'))
                # convSVGtoPNG(os.path.join(directory,long_id),os.path.join(directory,long_id))




#%% Post processing



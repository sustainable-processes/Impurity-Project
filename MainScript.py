# -*- coding: utf-8 -*-
"""
Created on Fri Oct  9 11:24:30 2020

This is the main script, and uses outputs from dataprocessing.py as well as
main functions defined in MainFunctions.py

@author: AADA01
"""

from MainFunctions import Chem, rdChemReactions,molfromsmiles,openpickle, getfragments,maprxn,rdMolDraw2D,drawReaction,parsemap,get_changed_atoms, os, writetofile,getlist,convSVGtoPNG, balance_stoichiometry,writepickle

#%% Importing libraries directly from DataProcessing.py

from DataProcessing import rxnlib, smles,analogue_rxns, template_dict  #To reload module just clear namespace and press f5 or runfile, this will cause UMR to autoreload


#%% Extracting template
# for rxnid,rxn in template_dict.items():
    




#%% Impurity Prediction

#%%  OR Reading libraries from file
# directory=os.path.join(os.getcwd(),'Dictionaries')
# rxnlib=openpickle(os.path.join(directory,'rxnlib.pickle'))
# smles=openpickle(os.path.join(directory,'smles.pickle'))



# {k:adict[k] for k in ('key1','key2','key99') if k in adict #Put earlier if want to accelerate code



#%% Writing images to file (remove comment if images change)
# directory=os.path.join(os.getcwd(),'Images','Sketches_ReactantOnly')
# for rxnid in rxnlib.keys():
#     writetofile(rxnlib[rxnid]['Sketch2'],os.path.join(directory,rxnid+'.svg'))
#     convSVGtoPNG(os.path.join(directory,rxnid),os.path.join(directory,rxnid))

# Compare atom list of carrier fragments to which reactants belong with reaction center
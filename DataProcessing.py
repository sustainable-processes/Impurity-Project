# -*- coding: utf-8 -*-
"""
Created on Fri Oct  9 10:39:17 2020

This script processes candidate reactions into a dictionary for easier referencing. Also 
produces an ID list which is used to retrieve SMILES strings and mole files from the server.

@author: AADA01
"""

from MainFunctions import openpickle

#%% CANDIDATE REACTION PROCESSING (0 to 14; refer to excel sheet)
#  Data is in pickle format (need to deserialize). No need pandas or numpy, can 
# directly read.This cell assembles a dictionary based on candidate reactions extracted from
# Reaxys. ONLY RUN ONCE IF candidate_rxns.pkl changes

candidate_rxns=openpickle('candidate_rxns.pkl')
rxn_list=candidate_rxns['rxns'][1:] #Removing header (Contains formatting)
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


analogue_compds=openpickle('analogue_compounds.pkl') #Not needed. For reference only

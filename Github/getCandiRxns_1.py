#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 25 16:39:23 2020

@author: Guo

searching for candidate reactions
"""
def checkRxts(rxts, comp_pools):
    # this function determines if the reaction records should be kept or not
    # rxts: a list of rxtIDs
    # comp_pools: dict of dict, { 'rxt/pro_smi':{'carr_frag_smarts': ['rxtID', ...], ...} , ...}
    #               { 'NCCCOC':{'CCOC': ['5068', ...], ...} , ...}
    # return: (bool, dict)
    #       bool, True if condition 1: each pro can be found in an analogue comp_pool
    #               and condition 2: no pros share the same analogue comp_pool
    #       dict: rxtsAssign = {rxtID: (rxt/pro_smi_fromQueryRxn, carr_frag_smi), ...}
    #             this indicate a rxt in a rxn under which comp_pool
    rxtsAssign = {k: []  for k in rxts}
    pool_records = set([]) # {(rxt/pro_smi, carr_frag_smi), ...} for condition 2
    for rxt in rxts: # rxt: rxtID '1234'
        # iterate through all pools under all compounds from the query rxn
        for k1 in comp_pools.keys(): # k1: rxt/pro_smi, 'NCCCCOC'
            frag_pool = comp_pools[k1] # frag_pool: {'[#7](-[#6](:[#6])': ['5068', ...], ...}
            for k2 in frag_pool.keys():# k2, frag_smi, 'CCOC'
                if rxt in set(frag_pool[k2]):
                    k1k2 = (k1, k2) # (rxt/pro_smi, frag_smi)
                    if k1k2 in pool_records: # condition 2: no compounds share the same pool
                        return (False, rxtsAssign)
                    rxtsAssign[rxt] = rxtsAssign[rxt] + [k1k2]
                    pool_records.update(set([k1k2]))
        if len(rxtsAssign[rxt]) == 0: # condition 1: all rxt or reagent should be in one pool
            return (False, rxtsAssign)
    return (True, rxtsAssign)



"""load analogue compounds pools"""




""" 
get candidate reactions 
searching candi_rxns are too slow at this step, because rows are read one by one
possible solution: using SQL? 
or multiprocessing, read multiple lines batch by batch 
and then process using multiple cores
"""                
keys_Qcomps = list(comp_pools.keys())
keys_carriFrags = [list(comp_pools[_].keys()) for _ in keys_Qcomps]
keys_Qcomps_carriFrags = [[(keys_Qcomps[i],_) for _ in keys_carriFrags[i]] for i in range(len(keys_Qcomps))]
keys_Qcomps_carriFrags = [_ for sublist in keys_Qcomps_carriFrags for _ in sublist]
  
candi_rxns = []              
rxtsAssigns = {k: [] for k in keys_Qcomps_carriFrags}  
# rxtsAssigns: { ('rxt/pro_smi', 'carr_frag_smarts'): ['rxtID', ...], ... }
with open(rxnSource, 'r') as infile:
    infile.readline(); infile.readline()
    while True:
        rxnRecords = infile.readline()
        if rxnRecords == '':
            break
        else:
            rxnRecordslist = rxnRecords.split('\t')
            rxts = rxnRecordslist[1].split(',')
            rxts = rxts + rxnRecordslist[11].split(',') # <<< include reagents also
            flag, rxts_assign = checkRxts(rxts, comp_pools) 
            # flag: bool; rxts_assign: {'15752734': [(rxt/pro_smi_fromQueryRxn, carr_frag_smi), ...], ...}
            if flag:
                candi_rxns = candi_rxns + [rxnRecords]
                for k1 in rxts_assign.keys():
                    for k2 in rxtsAssigns.keys():
                        for assign in rxts_assign[k1]:
                            if assign == k2:
                                rxtsAssigns[k2] = rxtsAssigns[k2] + [k1] 
                
                         









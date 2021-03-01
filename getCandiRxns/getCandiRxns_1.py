#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 25 16:39:23 2020

@author: Guo

searching for candidate reactions
"""

import os
import copy
from collections import namedtuple
from rdkit import Chem
from rdkit import DataStructs
import rdkit.Chem.Descriptors
from rdkit.Chem.Fingerprints import FingerprintMols
from FindFunctionalGroups import identify_functional_groups as IFG
import multiprocessing
from functools import partial
import gc



def multiprocessor(func, lsOrArray, numCores, unlist=True, **kwargs):
    """
    func, list/np.array, int, args to func -> list
    func: function that process list or np.array
    lsOrArray: list or array fed to func, list will be segraged into shorter list, 
                NOTE array will be segragated based on rows, row based multiprocessing
    numCores: number of cores
    **kwargs: keyword args to func
    [sliced_array[:,:3][_,:]]
    """
    core_pool = multiprocessing.Pool(processes=numCores, maxtasksperchild=1000)  # maxtasksperchild=10000
    len_obj = len(lsOrArray)
    sizeChunk = len_obj//numCores+1
    if isinstance(lsOrArray, list):
        chunks = [lsOrArray[_:_+sizeChunk] for _ in range(0, len_obj, sizeChunk)]
    elif isinstance(lsOrArray, np.ndarray):
        idxChunks = []
        for _ in range(0, len_obj, sizeChunk):
            idxChunks = idxChunks + [list(range(len_obj))[ _ : _ + sizeChunk]]
        chunks = []
        for _ in idxChunks:
            chunks = chunks + [lsOrArray[_,:]]  
    # multiprocessing start here
    results = core_pool.map(partial(func, **kwargs), chunks)
#    print('multiprocessor got results')
    core_pool.close()   
#    print('multiprocessor join results')
    core_pool.join()
#    print('multiprocessor join done')
    if unlist:
#        print('multiprocessor running unlisting')
        results = [_ for sublist in results for _ in sublist] # if sublist is not None
#        print('multiprocessor running done')
    gc.collect()
    return results

def getMols(IDs):
    '''
    Retrieves smiles strings from a set of Reaxys IDs using the cambridge server
    Connect via VPN

    Parameters
    ----------
    IDs : List
        List of Reaxys substance IDs

    Returns
    -------
    mols: List
        List of molecule files

    '''
    str_cwd = os.getcwd()
    os.chdir('/home/projects/graph/data/')
    folderNames = [_ for _ in os.listdir('.') if os.path.isdir(_)] # folder name and IDslist file name with .dat are same 
    IDslen = len(IDs)
    mols = [None]*IDslen
    for i in range(IDslen):
        for folderName in folderNames:
            try:
                os.chdir('/home/projects/graph/data/' + folderName)
                mols[i] = Chem.MolFromMolFile(IDs[i])
                break    # if get mol
            except:
                continue
    os.chdir(str_cwd)
    return mols 

def getIDsByFrag_v1(IDs, patt, MwThresh, addHs=False):
    # memo save version
    # list, mol, int -> list
    # get compounds based on pattern of fragment
    # get mols of IDs
    # IDs: list of IDs; patt: pattern to match eg.: Chem.MolFromSmarts('[C;R]-[OH]'), MwThresh: Mw threshold e.g. 1000
    # return: list of qualified IDs
    IDslen = len(IDs)
    compsWithFrag = [None]*IDslen
    for i in range(IDslen):
        mol = getMols([IDs[i]])[0] # may contain None
    # find compounds wtht he patt or fragment   
        if mol is not None:
            if addHs:
                mol = Chem.rdmolops.AddHs(mol)
            if Chem.Descriptors.MolWt(mol) < MwThresh:
                if len(mol.GetSubstructMatch(patt))!=0:
                    compsWithFrag[i] = IDs[i]
    compsWithFrag = [_ for _ in compsWithFrag if _ is not None]
    return  compsWithFrag

def getSimilarComp(compPool, targetComp, thresh, UseCompID=True):
    # list of IDs, str ID, thresh, type of targetComp-> list of IDs
    # get similar compounds of target compound based on tanimoto
    # compPool: list of compounds to compare with targetComp
    # targetComp: str ID of the targetComp
    # thresh: float, threshold to decide if to keep the compound as similar one
    # UseCompID: bool, True -> targetComp is a reaxys ID, False will be smiles str
    simComps = []
    if UseCompID:
        targetMol = getMols([targetComp])[0]
    else: # should be a smiles str
        targetMol = Chem.MolFromSmiles(targetComp)
    target_fp = FingerprintMols.FingerprintMol(targetMol)
    for comp in compPool:
        compMol = getMols([comp])[0]
        comp_fp = FingerprintMols.FingerprintMol(compMol)
        simValue = DataStructs.FingerprintSimilarity(target_fp, comp_fp)
        if simValue >= thresh:
            simComps = simComps + [comp]
    return simComps


def getCarrierFrags(smi, size, resFormat='smarts', addHs=True):
    """
    str (smiles), int -> list (str_smiles/smarts) 
    smi: str, smiles of a compound
    size: size of carrier fragments, number atoms
    resFormat: 'smiles' or 'smarts'
    addHs: bool, if True H will be consider for generation of substructure, recommend True
            otherwise terminal atoms of a molecule are not differentiated with other atoms
    this function return list of strings of smarts representing carrier frags with miniSize of size
    a carrier frag carries a functional group defined by using Ertl's method
    find out the bonds to cut
    define a cutter to cut out the target fragments based on bonds to cut
    1) get a list of function groups using IFG
    2) if IFG list empty meaning no functional groups, directly return compound smiles
    3) evaluate the size of current fragment,  
       if frag_size >= size, directly return the comp smile
       if frag_size < miniSizeFrag, expand the fragment by searching for neighors, only use terminal atoms to search
       repeat expansion till frag_size >= size
    """
    mol = Chem.MolFromSmiles(smi)
    if addHs:
        mol = Chem.AddHs(mol)
    # -- get the list of functional groups FG
    IFG_ls = IFG(mol) # e.g., [IFG(atomIds=(1, 4, 7), atoms='NC=O', type='cNC(C)=O'), IFG(atomIds=(10,), atoms='O', type='cO')]
    # if IFG_ls is empty, directly return this compounds
    if len(IFG_ls) == 0:
        if resFormat=='smiles':
            return smi
        elif resFormat=='smarts':
            return Chem.MolToSmarts(Chem.MolFromSmiles(smi))
    # -- get atomIDs (FGs_atomIDs_expan) and terminalAtomIDs (FGs_terminal_atomIDs) for all frags
    FGs_atomIDs = [_.atomIds for _ in IFG_ls] # e.g., [(1, 4, 7), ...]
    n_FGs = len(FGs_atomIDs)
    FGs_atomIDs_expan = [None]*n_FGs # expan all FGs that < size # e.g., [[1, 4, 7, 8, 9], ...]
    # terminal atoms: atoms on which bonds to cut will be searched
    FGs_terminal_atomIDs = [None]*n_FGs # [[2, 3, 8], ...] or [[], [], ...], note [2, 3, 8] are terminal of comp not frag
    for i in range(n_FGs):
        # initialization before search and expand fragments
        FG_atomIDs = list(FGs_atomIDs[i]) # e.g., [1, 4, 7]
        FG_size = len(FG_atomIDs)
        # find initial terminal atoms
        # terminal atoms: atoms with neis that not belong to FG_atomIDs
        # note terminal atoms here are terminal of fragments not comp
        # C=O has no fragment terminal, but both C and O are comp terminals
        # there are 3 cases that no fragment terminal ie, FG_terminal_atomIDs = [] one case FG_terminal_atomIDs not empty
        # case 1: FG_size < size, and FG_terminal_atomIDs = [], ie., comp too small, no chane to expand e.g., C=O, HCl, HNO3
        # case 2: FG_size < size, and FG_terminal_atomIDs = [2,3,8], 
        #         got chane to expand, but after expand still cannot reach size, e.g., OCC=O
        #         eventually FG_terminal_atomIDs = []
        # case 3: FG_size >= size, and FG_terminal_atomIDs = [], cos the whole big comp is a functional group
        # for all above cases the original comp will be return
        # two case FG_terminal_atomIDs is not empty, cutting will performed based on FG_terminal_atomIDs
        # case 4: FG_size >= size, and FG_terminal_atomIDs = [2,3,8], just cut using FG_terminal_atomIDs
        # case 5: FG_size < size, and FG_terminal_atomIDs = [2,3,8], after expand, size reached and FG_terminal_atomIDs changed
        FG_terminal_atomIDs = [] # e.g., [2, 3, 8] or [], HNO3 or C=O -> []
        for atomID in FG_atomIDs:
            neis_IDs = [_.GetIdx() for _ in mol.GetAtomWithIdx(atomID).GetNeighbors()]
            if len(set(neis_IDs) - set(FG_atomIDs)) != 0:
                FG_terminal_atomIDs = FG_terminal_atomIDs + [atomID]
        # perform search / expand fragments
        if (FG_size >= size) | (len(FG_terminal_atomIDs) == 0): # [],  case 1 cannot expand, although FG_size < size
            # make sure all elements in FGs_atomIDs_expan are lists
            FGs_atomIDs_expan[i] = FG_atomIDs
            FGs_terminal_atomIDs[i] = FG_terminal_atomIDs 
            # case 3 still could be [] even if FG_size >= size, though less likely
        else:
            # init FG to be expand
            FG_expan_atomIDs = copy.deepcopy(FG_atomIDs) # e.g., [1, 4, 7]
            for rep in range(size): # max repeat size times, since repeat size should reach the size
                FG_expan_atomIDs_old = copy.deepcopy(FG_expan_atomIDs)
                # not all atoms need be searched for neis, only ones not in FG_expan_atomIDs for each epoch
                for atomID in FG_terminal_atomIDs:
                    neis_IDs = [_.GetIdx() for _ in mol.GetAtomWithIdx(atomID).GetNeighbors()] 
                    FG_expan_atomIDs = FG_expan_atomIDs + neis_IDs
                FG_expan_atomIDs = list(set(FG_expan_atomIDs))
                FG_terminal_atomIDs = list(set(FG_expan_atomIDs)-set(FG_expan_atomIDs_old))
                if len(FG_terminal_atomIDs)==0: # cannot expand due to small comp size
                    FGs_atomIDs_expan[i] = FG_expan_atomIDs
                    FGs_terminal_atomIDs[i] = FG_terminal_atomIDs # [] case 2
                    break # for expan of next FG
                if len(FG_expan_atomIDs) >= size: # reach the desired size
                    FGs_atomIDs_expan[i] = FG_expan_atomIDs
                    FGs_terminal_atomIDs[i] = FG_terminal_atomIDs
                    break 
    # it seems not possible that [[2, 3, 8], [], ...], if there one [] then all should be [], 
    # either [[]] or [[], [], ...] (i.e., small comp with single or multiple FGs)     
    if len(FGs_terminal_atomIDs[0])==0: #sum([len(_) for _ in FGs_terminal_atomIDs]) == 0:  
        if resFormat=='smiles':
            return smi
        elif resFormat=='smarts':
            return Chem.MolToSmarts(Chem.MolFromSmiles(smi))
    else:
        FGs_strs = [] # smiles or smarts
        for FG_atomIDs in FGs_atomIDs_expan:
            if resFormat=='smiles':
                FGs_str = Chem.MolFragmentToSmiles(mol, FG_atomIDs, canonical=True)
                FGs_strs = FGs_strs + [FGs_str]
            elif resFormat=='smarts':
                FGs_str = Chem.MolFragmentToSmarts(mol, FG_atomIDs, isomericSmarts=False)
                FGs_strs = FGs_strs + [FGs_str]
        return FGs_strs
    
def getCompPool(comp, source, fragSize, smiThresh, MwThresh, numCores, addHs=True):
    """
    str, lst, int, float, str, int -> dict
    comp: str, smiles of a comp
    source: list, reaxys IDs of compounds e.g.,: ['1234', ...]
    fragSize: int, mininum fragment size
    smiThresh: float, 0 to 1, recommend 0.5, similarity threshold
    MwThresh: int/float, threshold for molecular weights
    numCores: int, number of cpus
    addHs: normally true, since carrier fragments from getCarrierFrags normally count on H
    return: CompPool
        dict: qualified compounds, keys: smiles of frag, values: list of compIDs
                e.g.,: {'CCOC': ['5068', ...], ...}
    dependencies:
        getCarrierFrags()
        multiprocessor()
        getSimilarComp()
        getIDsByFrag_v1()
    """       
    # get fragments of each comp
    frags = getCarrierFrags(smi=comp, size=fragSize)
    # initiate the dict to return
    CompPool = {k: None for k in frags}
    for frag in frags:
        # find analogue comps contains frag
        patt_frag = Chem.MolFromSmarts(frag)
        if len(source) >= 500 and numCores!=0: # just to determined if multiprocessing is necessary
            AnaComps = multiprocessor(func=getIDsByFrag_v1, lsOrArray=source, 
                                        numCores=numCores, unlist=True, patt=patt_frag, 
                                        MwThresh=MwThresh, addHs=addHs)
        else:
            AnaComps = getIDsByFrag_v1(source, patt_frag, MwThresh, addHs)
        # keep comps pass smiThresh
        if len(AnaComps) >=500 and numCores!=0:
            AnaComps = multiprocessor(func=getSimilarComp, lsOrArray=AnaComps, 
                                                numCores=numCores, unlist=True, 
                                                targetComp=comp, thresh=smiThresh, UseCompID=False)
        else:
            AnaComps = getSimilarComp(AnaComps, comp, smiThresh, UseCompID=False)
        CompPool[frag] = AnaComps
    return CompPool    



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

""" 
get candidate reactions 
searching candi_rxns are too slow at this step, because rows are read one by one
possible solution: using SQL? 
or multiprocessing, read multiple lines batch by batch 
and then process using multiple cores
"""     
def get_candirxns(comp_pools,inputdir,rxnSource):
    keys_Qcomps = list(comp_pools.keys())
    keys_carriFrags = [list(comp_pools[_].keys()) for _ in keys_Qcomps]
    keys_Qcomps_carriFrags = [[(keys_Qcomps[i],_) for _ in keys_carriFrags[i]] for i in range(len(keys_Qcomps))]
    keys_Qcomps_carriFrags = [_ for sublist in keys_Qcomps_carriFrags for _ in sublist]
    candi_rxns = []              
    rxtsAssigns = {k: set() for k in keys_Qcomps_carriFrags}  
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
                # rxts = rxts + rxnRecordslist[11].split(',') # <<< include reagents also
                # reagents will not be checked if they contain carrierfrag or not but
                # reagents provided by user will be check if they exist in current reaction
                # unfortunately, if no reagents, there is no ''
                allReagsIn = all([_ in reags for _ in ReagIDs])
                # in case important reactions are missed, if number of references higher than certain threshold
                # here >=5, the reaction will be kept
                numRef = rxnRecordslist[3] # number of references for this reaction, higher means more reliable
                if numRef != '': # <<< may possibly got other unexpected string
                    numRef = int(numRef)
                flag, rxts_assign = checkRxts(rxts, comp_pools) 
                # flag: bool; rxts_assign: {'15752734': [(rxt/pro_smi_fromQueryRxn, carr_frag_smi), ...], ...}
                if (flag & allReagsIn) | (flag & numRef >=5):
                    candi_rxns = candi_rxns + [rxnRecords]
                    for k1 in rxts_assign.keys():
                        for k2 in rxtsAssigns.keys():
                            for assign in rxts_assign[k1]:
                                if assign == k2:
                                    rxtsAssigns[k2].update([k1])
    return candi_rxns, rxtsAssigns
                         









# %load ./FunctionsDB.py

import shutil
import itertools
import copy
from MainFunctions import writepickle, openpickle
from collections import Counter, namedtuple

from rdkit import Chem  # Importing RDKit
from rdkit.Chem.Draw import rdMolDraw2D  # Drawing 2D molecules/reactions
from rdkit.Chem import rdChemReactions  # Reaction processing
from rdkit.Chem import Draw  # For drawing molecules/reactions
from rdkit.Chem import AllChem  # Overall support
from FindFunctionalGroups import identify_functional_groups as IFG

import sqlite3
import os
import dask.delayed as delayed
import dask.dataframe as dd
import dask.array as da
import dask.bag as dba
import numpy as np
import pandas as pd

#%% Building basic substance database

def info(molfile):
    mol = Chem.MolFromMolFile(molfile)
    if mol:
        Chem.SanitizeMol(mol)
        mol.UpdatePropertyCache(strict=False)
        smils = Chem.MolToSmiles(mol)
        return mol, smils
    else:
        return molfile


def basic(ID, folder):
    if str.isdecimal(ID):
        molfileaddress = folder+os.sep+ID
        try:
            res = info(molfileaddress)
        except Exception as e:
            error = e
            compaddrs = {'SubstanceID': int(
                ID), 'MolFileAddress': molfileaddress, 'Error': error}
        else:
            if type(res) == tuple:
                smiles = res[1]
                compaddrs = {'SubstanceID': int(
                    ID), 'MolFileAddress': molfileaddress, 'Smiles': smiles}
            else:
                error = 'Valence error'
                compaddrs = {'SubstanceID': int(
                    ID), 'MolFileAddress': molfileaddress, 'Error': error}
    else:
        compaddrs = {}
    return compaddrs


def basicgroup(molfilelist, folder):
    return [basic(ID, folder) for ID in molfilelist]


def substancedblist(folderName, partitions):
    dem = os.listdir(folderName)
    b = db.from_sequence(dem, npartitions=partitions)
    dflist = b.map_partitions(basicgroup, folderName).compute()
    return dflist


#%% Fragment detection
def getCarrierFrags0(smi,expand=1,userinput='smiles',resFormat='smarts',addHs=True):
    """
    str (smiles), int -> list (str_smiles/smarts) 
    smi: str, smiles of a compound
    size: Level to expand functional group by (eg./ 1 will expand to first degree neighbours, 2 will expand to second degree neighbors etc.)
    resFormat: 'smiles' or 'smarts' 
    addHs: bool, if True H will be consider for generation of substructure, recommend True
            otherwise terminal atoms of a molecule are not differentiated with other atoms
    this function return list of strings of smarts representing carrier frags with miniSize of size
    a carrier frag carries a functional group defined by using Ertl's method
    find out the bonds to cut define a cutter to cut out the target fragments based on bonds to cut
    1) get a list of function groups using IFG
    2) if IFG list empty meaning no functional groups, directly return compound smiles
    3) expand to nearest neighbors based on expand value
    """
    if userinput=='smiles':
        mol = Chem.MolFromSmiles(smi)
        Chem.SanitizeMol(mol)
        mol.UpdatePropertyCache(strict=False)
    else:
        mol=smi
    if addHs:
        mol = Chem.AddHs(mol)
    # -- get the list of functional groups FG
    # e.g., [IFG(atomIds=(1, 4, 7), atoms='NC=O', type='cNC(C)=O'), IFG(atomIds=(10,), atoms='O', type='cO')]
    IFG_ls = IFG(mol)
    # if IFG_ls is empty, directly return this compounds
    if len(IFG_ls) == 0:
        if resFormat == 'smiles':
            if userinput!='smiles':
                return Chem.MolToSmiles(mol)
            else:
                return smi
        elif resFormat == 'smarts':
            if userinput!='smiles':
                return Chem.MolToSmarts(mol)
            else:
                return Chem.MolToSmarts(Chem.MolFromSmiles(smi))
    # -- get atomIDs (FGs_atomIDs_expan) and terminalAtomIDs (FGs_terminal_atomIDs) for all frags
    FGs_atomIDs = [_.atomIds for _ in IFG_ls]  # e.g., [(1, 4, 7), ...]
    n_FGs = len(FGs_atomIDs)
    # expan all FGs that < size # e.g., [[1, 4, 7, 8, 9], ...]
    FGs_atomIDs_expan = [None]*n_FGs
    # terminal atoms: atoms on which bonds to cut will be searched
    # [[2, 3, 8], ...] or [[], [], ...], note [2, 3, 8] are terminal of comp not frag
    FGs_terminal_atomIDs = [None]*n_FGs
    for i in range(n_FGs):
#         breakpoint()
        # initialization before search and expand fragments
        FG_atomIDs = list(FGs_atomIDs[i])  # e.g., [1, 4, 7]
        FG_size = len(FG_atomIDs)
        FG_terminal_atomIDs = []  # e.g., [2, 3, 8] or [], HNO3 or C=O -> [] Terminal IDs of fragments
        for atomID in FG_atomIDs:
            neis_IDs = [_.GetIdx() for _ in mol.GetAtomWithIdx(atomID).GetNeighbors()]
            if len(set(neis_IDs) - set(FG_atomIDs)) != 0:
                FG_terminal_atomIDs = FG_terminal_atomIDs + [atomID]
        if len(FG_terminal_atomIDs) == 0:
            # make sure all elements in FGs_atomIDs_expan are lists
            FGs_atomIDs_expan[i] = FG_atomIDs
            FGs_terminal_atomIDs[i] = FG_terminal_atomIDs
            # case 3 still could be [] even if FG_size >= size, though less likely
        else:
            FG_expan_atomIDs = copy.deepcopy(FG_atomIDs)  # e.g., [1, 4, 7]
             # max repeat size times, since repeat size should reach the size
            for rep in range(expand):
                FG_expan_atomIDs_old = copy.deepcopy(FG_expan_atomIDs)
                # not all atoms need be searched for neis, only ones not in FG_expan_atomIDs for each epoch
                for atomID in FG_terminal_atomIDs:
                    neis_IDs = [_.GetIdx() for _ in mol.GetAtomWithIdx(atomID).GetNeighbors()]
                    FG_expan_atomIDs = FG_expan_atomIDs + neis_IDs
                FG_expan_atomIDs = list(set(FG_expan_atomIDs))
                FG_terminal_atomIDs = list(set(FG_expan_atomIDs)-set(FG_expan_atomIDs_old))
                FGs_atomIDs_expan[i] = FG_expan_atomIDs
                FGs_terminal_atomIDs[i] = FG_terminal_atomIDs  # [] case 2
                if len(FG_terminal_atomIDs) == 0:  # cannot expand due to small comp size
                    break  # for expan of next FG
    # it seems not possible that [[2, 3, 8], [], ...], if there one [] then all should be [],
    # either [[]] or [[], [], ...] (i.e., small comp with single or multiple FGs)
    # sum([len(_) for _ in FGs_terminal_atomIDs]) == 0:
#     breakpoint()
    if len(FGs_terminal_atomIDs[0]) == 0:
        if resFormat == 'smiles':
            if userinput!='smiles':
                return Chem.MolToSmiles(mol)
            else:
                return smi
        elif resFormat == 'smarts':
            if userinput!='smiles':
                return Chem.MolToSmarts(mol)
            else:
                return Chem.MolToSmarts(Chem.MolFromSmiles(smi))
    else:
        FGs_strs = []  # smiles or smarts
        for FG_atomIDs in FGs_atomIDs_expan:
            if resFormat == 'smiles':
                FGs_str = Chem.MolFragmentToSmiles(mol, FG_atomIDs, canonical=True)
                FGs_strs = FGs_strs + [FGs_str]
            elif resFormat == 'smarts':
                FGs_str = Chem.MolFragmentToSmarts(mol, FG_atomIDs, isomericSmarts=False)
                FGs_strs = FGs_strs + [FGs_str]
        return FGs_strs
    
    
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
    find out the bonds to cut define a cutter to cut out the target fragments based on bonds to cut
    1) get a list of function groups using IFG
    2) if IFG list empty meaning no functional groups, directly return compound smiles
    3) evaluate the size of current fragment,  
       if frag_size >= size, directly return the comp smile
       if frag_size < miniSizeFrag, expand the fragment by searching for neighors, only use terminal atoms to search
       repeat expansion till frag_size >= size
    """
    mol = Chem.MolFromSmiles(smi)
    Chem.SanitizeMol(mol)
    mol.UpdatePropertyCache(strict=False)
    if addHs:
        mol = Chem.AddHs(mol)
    # -- get the list of functional groups FG
    # e.g., [IFG(atomIds=(1, 4, 7), atoms='NC=O', type='cNC(C)=O'), IFG(atomIds=(10,), atoms='O', type='cO')]
    IFG_ls = IFG(mol)
    # if IFG_ls is empty, directly return this compounds
    if len(IFG_ls) == 0:
        if resFormat == 'smiles':
            return smi
        elif resFormat == 'smarts':
            return Chem.MolToSmarts(Chem.MolFromSmiles(smi))
    # -- get atomIDs (FGs_atomIDs_expan) and terminalAtomIDs (FGs_terminal_atomIDs) for all frags
    FGs_atomIDs = [_.atomIds for _ in IFG_ls]  # e.g., [(1, 4, 7), ...]
    n_FGs = len(FGs_atomIDs)
    # expan all FGs that < size # e.g., [[1, 4, 7, 8, 9], ...]
    FGs_atomIDs_expan = [None]*n_FGs
    # terminal atoms: atoms on which bonds to cut will be searched
    # [[2, 3, 8], ...] or [[], [], ...], note [2, 3, 8] are terminal of comp not frag
    FGs_terminal_atomIDs = [None]*n_FGs
    for i in range(n_FGs):
        # initialization before search and expand fragments
        FG_atomIDs = list(FGs_atomIDs[i])  # e.g., [1, 4, 7]
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
        FG_terminal_atomIDs = []  # e.g., [2, 3, 8] or [], HNO3 or C=O -> []
        for atomID in FG_atomIDs:
            neis_IDs = [_.GetIdx() for _ in mol.GetAtomWithIdx(atomID).GetNeighbors()]
            if len(set(neis_IDs) - set(FG_atomIDs)) != 0:
                FG_terminal_atomIDs = FG_terminal_atomIDs + [atomID]
        # perform search / expand fragments
        # [],  case 1 cannot expand, although FG_size < size
        if (FG_size >= size) | (len(FG_terminal_atomIDs) == 0):
            # make sure all elements in FGs_atomIDs_expan are lists
            FGs_atomIDs_expan[i] = FG_atomIDs
            FGs_terminal_atomIDs[i] = FG_terminal_atomIDs
            # case 3 still could be [] even if FG_size >= size, though less likely
        else:
            # init FG to be expand
            FG_expan_atomIDs = copy.deepcopy(FG_atomIDs)  # e.g., [1, 4, 7]
            # max repeat size times, since repeat size should reach the size
            for rep in range(size):
                FG_expan_atomIDs_old = copy.deepcopy(FG_expan_atomIDs)
                # not all atoms need be searched for neis, only ones not in FG_expan_atomIDs for each epoch
                for atomID in FG_terminal_atomIDs:
                    neis_IDs = [_.GetIdx() for _ in mol.GetAtomWithIdx(atomID).GetNeighbors()]
                    FG_expan_atomIDs = FG_expan_atomIDs + neis_IDs
                FG_expan_atomIDs = list(set(FG_expan_atomIDs))
                FG_terminal_atomIDs = list(set(FG_expan_atomIDs)-set(FG_expan_atomIDs_old))
                if len(FG_terminal_atomIDs) == 0:  # cannot expand due to small comp size
                    FGs_atomIDs_expan[i] = FG_expan_atomIDs
                    FGs_terminal_atomIDs[i] = FG_terminal_atomIDs  # [] case 2
                    break  # for expan of next FG
                if len(FG_expan_atomIDs) >= size:  # reach the desired size
                    FGs_atomIDs_expan[i] = FG_expan_atomIDs
                    FGs_terminal_atomIDs[i] = FG_terminal_atomIDs
                    break
    # it seems not possible that [[2, 3, 8], [], ...], if there one [] then all should be [],
    # either [[]] or [[], [], ...] (i.e., small comp with single or multiple FGs)
    # sum([len(_) for _ in FGs_terminal_atomIDs]) == 0:
    if len(FGs_terminal_atomIDs[0]) == 0:
        if resFormat == 'smiles':
            return smi
        elif resFormat == 'smarts':
            return Chem.MolToSmarts(Chem.MolFromSmiles(smi))
    else:
        FGs_strs = []  # smiles or smarts
        for FG_atomIDs in FGs_atomIDs_expan:
            if resFormat == 'smiles':
                FGs_str = Chem.MolFragmentToSmiles(mol, FG_atomIDs, canonical=True)
                FGs_strs = FGs_strs + [FGs_str]
            elif resFormat == 'smarts':
                FGs_str = Chem.MolFragmentToSmarts(mol, FG_atomIDs, isomericSmarts=False)
                FGs_strs = FGs_strs + [FGs_str]
        return FGs_strs
    
#%% Adding fragment smarts column accounting for mixtures
def getfrags(series,expand=1): #natoms changed to expand
    smiles=series['Smiles']
    if smiles=='Error':
        return 'Error','Error'
    if series['>1 Compound']==True: #This compound is a mixture. Need to split and apply getcarrierfrags to each smiles
        fragsmarts=getmixturefrags(smiles,expand=expand)
        fragsmiles=getmixturefrags(smiles,expand=expand,resFormat='smiles')
    else:
        try:
            fragsmarts=getCarrierFrags0(smiles,expand=expand)
            fragsmiles=getCarrierFrags0(smiles,expand=expand,resFormat='smiles')
        except Exception:
            return 'Error','Error'
    if type(fragsmarts)!=list:
        fragsmarts=[fragsmarts]
        fragsmiles=[fragsmiles]
    
#     fragsmiles=[Chem.MolToSmiles(Chem.MolFromSmarts(fragsmart)) for fragsmart in fragsmarts] #Does not capture aromaticity
    return fragsmiles,fragsmarts
        
    
def getfragpartition(partition,natoms):
    return partition.apply(getfrags,natoms=natoms,axis=1)


#%% Adding fragment smiles column
def getfragsmiles(fragsmarts):
    if type(fragsmarts)==list:
        fragsmiles=[Chem.MolToSmiles(Chem.MolFromSmarts(fragsmart)) for fragsmart in fragsmarts]    
    elif fragsmarts=='Error':
        return 'Error'
    else:
         fragsmiles=Chem.MolToSmiles(Chem.MolFromSmarts(fragsmarts))
    return fragsmiles

def getsmiles(series):
    fragsmarts=series['FragmentSmarts']
    return getfragsmiles(fragsmarts)

def getsmilespartition(partition):
    return partition.apply(getsmiles,axis=1)

#%% Creating mixture column (True if mixture, False if not mixture, Error if smiles not present)

def mixtures(smiles):
    if smiles=='Error':
        return 'Error'
    elif len(smiles.split('.'))>1:
        return True
    else:
        return False

def findMixtures(series):
    smiles=series['Smiles']
    return mixtures(smiles)
    
def findMixturespartition(partition):
       return partition.apply(findMixtures,axis=1)

#%% Changing fragment entries in mixture rows

def getmixturefrags(mixsmiles,expand=1,resFormat='smarts', addHs=True): 
    try:
        res=[]
        for smiles in mixsmiles.split('.'):
            reslist=getCarrierFrags0(smiles,expand=expand,resFormat=resFormat,addHs=addHs)
            if type(reslist)!=list:
                res+=[reslist]
            else:
                res+=reslist
    except Exception:
        return 'Error'
    else:
        return res    

def getMixturefrags(series,expand=1):
    mixsmiles=series['Smiles'] # add .values[0] if column is a multiindex, otherwise droplevel = 1 to remove list
    return getmixturefrags(mixsmiles,expand=expand)  

def getMixturefragspartition(partition,natoms):
    return partition.apply(getMixturefrags,expand=expand,axis=1)

# def collapsepartition(partition): #Doesn't work yet. Once exploded it is extremely time-consuming to collapse the dataframe
#     temp=partition.groupby(['Smiles']).agg([list])
#     return temp

    
#%% Joining columns to a dataframe

def joindf(seriesdf,DB,explodeDB=None):
    if seriesdf.index.name!=DB.index.name or seriesdf.index.names!=DB.index.names:
        if DB.index.name or DB.index.names:
            DB.reset_index(inplace=True)
        if seriesdf.index.name or seriesdf.index.names:
            seriesdf.reset_index(inplace=True)        
    DB=DB.join(seriesdf)
    if explodeDB:
        DB=DB.explode(explodeDB)
    return DB




def buildfragdb(sdb=None,sdbd=None,writesdbd=False, \
                sdbdc=None,fragseries=None,natoms=None, \
                writefragseries=False, fdbm=None, fdb=None, dfdb=None, \
                fragsmiles=False, mixtures=False, mixturefrags=False, \
                index=None, writefdb=False,writedfdb=False):
    
    # Note: cluster and client must be initiated for this function to work. Substance database, either dask or pandas
    # should also be loaded. It is advised to avoid loading dask dataframes due to high memory usage. It is recommended to
    #persist and create dask dataframes outside this function, otherwise overheads may be added. Never reset dask dataframe
    #index or else it will reindex all partitions...always start with pandas, reindex and then convert to dask when large computations
    #need to be done.
    
    # sdb = substance database, sdbd = substance database dask, sdbdc = substance database dask cleaned 
    #(remove smile errors), # fragseries = series containing fragment information, natoms =  size of fragment, 
    # writefragseries = True if write to file false otherwise, fdbm = fragment database master, pandas version of sdbdc,
    # fdb = final fragment database exploded and unindexed pfdb = final pandas fragment dataframe, multiindexed, 
    # dfdb = final dask fragment dataframe
   
    # Step 1: Create dask dataframe from pandas substance database (output: sdbd)
    
    if sdb and not type(sdbd)==dd.core.DataFrame:
        sdbd = dd.from_pandas(sdb, npartitions=16)
        sdbd = client.persist(sdbd)
        print('Dask substance dataframe created and persisted')
        if writesdbd:
            sdbd.to_parquet(writesdbd)
            print('Dask substance dataframe written to file: ' + writesdb)
        if not natoms and not fragseries:
            print("Please specify size of fragments that should be retrieved")
            return sdbd
        
   # Step 2: Clean dask data frame and select only substance ID and smiles column, removing errors (output: sdbdc).

    if sdbd and not type(sdbdc)==dd.core.DataFrame:
        sdbd = client.persist(sdbd)
        sdbdc = sdbd.reset_index()[['SubstanceID', 'Smiles','>1 Compound']]
        sdbdc= sdbdc[sdbdc.Smiles!='Error']
        sdbdc=client.persist(sdbdc)
        print('Cleaned dask substance dataframe created and persisted')
        if not natoms and not type(fragseries)==pd.core.series.Series:
            print("Please specify size of fragments that should be retrieved")
            return sdbdc      
        
    # Step 3: Scrape 16 million compounds, and extract series of active fragments for each (output: fragseries)
    
    if not type(fragseries)==pd.core.series.Series and not type(fdb)==pd.core.frame.DataFrame:
        if not type(sdbdc)==dd.core.DataFrame:
            return "Please include a cleaned dask substance dataframe for fragment retrieval"
        if not natoms:
            return "Please specify size of fragments retrieved"
        sdbdc=client.persist(sdbdc)
#         if natoms==0:
#             name='ActiveFragmentSmarts'
#         else:
#            filename='CarrierFragments'+'(n='+str(natoms)+')'
        name='FragmentSmarts'
        fragseries=sdbdc.map_partitions(getfragpartition,natoms=natoms,meta=(name,'O')).compute()
        print('fragseries retrieved')
    
    if writefragseries:
        if not type(fragseries)==pd.core.series.Series:
            return "Supply fragment series to write to file"
        writepickle(fragseries,writefragseries)
        print('fragseries writted to file: '+ writefragseries)
      
     # Step 4: Prepare cleaned pandas dataframe for fragment series attachment (output: fdbm)  
        
    if not type(fdbm)==pd.core.frame.DataFrame and not type(fdb)==pd.core.frame.DataFrame:
        if not type(sdbdc)==dd.core.DataFrame:
            return "Please include a cleaned dask substance dataframe to which fragment information can be attached"
        fdbm=sdbdc.compute()
#         fdbm.reset_index(inplace=True)
#         fdbm.drop('index',axis=1,inplace=True)
        
    # Step 5: Attaching fragment series, generating an exploded fragment database (output: fdb)
    
    if type(fdbm)==pd.core.frame.DataFrame:
        if not type(fragseries)==pd.core.series.Series:
            return "Please include the fragment series that should be attached"
        fragdf=pd.DataFrame(fragseries)
        fdb=joindf(fragdf,fdbm,explodeDB=fragdf.columns[0])
        print("Unindexed fragment database completed")
    
    # Step 6: Adding additional columns, formatting, indexing and writing fragment database to file

    if fragsmiles:
        if not type(fdb)==pd.core.frame.DataFrame:
            return "Please supply fragment database to analyze"
        if not type(dfdb)==dd.core.DataFrame:
            if fdb.index.name or fdb.index.names:
                fdb.reset_index(inplace=True)
            dfdb=dd.from_pandas(fdb,npartitions=181)
            dfdb=client.persist(dfdb)
        fragseries=dfdb.map_partitions(getsmilespartition,meta=('FragmentSmiles','O')).compute()
        if writefragseries:
            if not type(fragseries)==pd.core.series.Series:
                return "Supply fragment series to write to file"
            writepickle(fragseries,writefragseries)
            print('fragseries writted to file: '+ writefragseries)    
        smilesdf=pd.DataFrame(fragseries)
        fdb=joindf(smilesdf,fdb)
        
        
    if mixtures:
        if not type(fdb)==pd.core.frame.DataFrame:
            return "Please supply fragment database to analyze"
        if not type(dfdb)==dd.core.DataFrame:
            if fdb.index.name or fdb.index.names:
                fdb.reset_index(inplace=True)
            dfdb=dd.from_pandas(fdb,npartitions=181)
            dfdb=client.persist(dfdb)
        fragseries=dfdb.map_partitions(findMixturespartition,meta=('>1 Compound','boolean'))
        if writefragseries:
            if not type(fragseries)==pd.core.series.Series:
                return "Supply fragment series to write to file"
            writepickle(fragseries,writefragseries)
            print('fragseries writted to file: '+ writefragseries) 
        mixturedf=pd.DataFrame(fragseries)
        fdb=joindf(mixturedf,fdb)
    
    if mixturefrags:
        #NOTE: If exploded, fragment database should be aggregated by Smiles/Smarts
        if not type(fdb)==pd.core.frame.DataFrame or '> 1 Compound' not in fdb.columns:
            return "Please supply fragment database to analyze, with mixture indication (specify mixtures=True)"
        if not natoms:
            return "Please specify size of fragments retrieved"
        if not type(dfdb)==dd.core.DataFrame:
            dfdb=dd.from_pandas(fdb,npartitions=16)
            dfdb=client.persist(dfdb)
        fragseries=dfdb.map_partitions(getMixturefragspartition,natoms=natoms,meta=('FragmentSmarts','O'))
        if writefragseries:
            if not type(fragseries)==pd.core.series.Series:
                return "Supply fragment series to write to file"
            writepickle(fragseries,writefragseries)
            print('fragseries writted to file: '+ writefragseries) 
        mixturesmarts=pd.DataFrame(fragseries)
        fdb=joindf(mixturesmarts,fdb)
        
        
    
# Deprecated..Dask dataframes take up too much memory when read from file

#     if dfdb:
#         if fdb.index.name:
#             fdb.reset_index(inplace=True)
#         fdb.set_index('SubstanceID',inplace=True)
#         dfdb=dd.from_pandas(fdb,npartitions=64)
# #         dfdb=dfdb.reset_index()
# #         dfdb=dfdb.set_index(dfdb.columns[dfdb.columns.str.contains('ragment')][0])
#         print('Dask fragment database created')
#     if writedfdb:
#         dfdb.to_parquet(writedfdb)
#         print('Dask fragment database written to file: ' + writedfdb)

        
    #%% Indexing dataframe
    
    if index:
        if not type(fdb)==pd.core.frame.DataFrame:
             return "Please supply fragment database to index"
        if fdb.index.name==index or fdb.index.names==index:
            print("Database is already indexed.") 
        elif fdb.index.name or fdb.index.names:
            fdb.reset_index(inplace=True)
            fdb.set_index([index],inplace=True)
        else:
            fdb.set_index([index],inplace=True)
        print('Pandas fragment database indexed')
    if writefdb:
        writepickle(fdb,writefdb)
        print('Pandas fragment database written to file: ' + writefdb)
    
    return fdb
    

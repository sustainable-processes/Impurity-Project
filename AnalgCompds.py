# %load ./AnalgCompds.py

from rdkit import Chem
from rdkit.Chem.Descriptors import MolWt
from MainFunctions import molfromsmiles,CustomError
from collections import Counter
import os
import sqlite3
import pandas as pd
from FindFunctionalGroups import identify_functional_groups as IFG
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit import DataStructs
import ray
import modin.pandas as mpd
from ReaxysAPIv2 import initray

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
        # initialization before search and expand fragments
        FG_atomIDs = list(FGs_atomIDs[i])  # e.g., [1, 4, 7]
        FG_size = len(FG_atomIDs)
        FG_terminal_atomIDs = []  # e.g., [2, 3, 8] or [], HNO3 or C=O -> [] Terminal IDs of fragments
        for atomID in FG_atomIDs:
            neis_IDs = [_.GetIdx() for _ in mol.GetAtomWithIdx(atomID).GetNeighbors()]
            if len(set(neis_IDs) - set(FG_atomIDs)) != 0:
                FG_terminal_atomIDs = FG_terminal_atomIDs + [atomID]
        if (expand==0) | (len(FG_terminal_atomIDs) == 0):
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
    
    

def processquery(userinput,expand=1,refquery={},debug=False):
    '''
    Takes in user input (reaction), splits into species, and creates dictionary with carrier fragments for each
    
    Expand determines degree of expansion from the identified functional group to generate carrier fragments
    refquery refers to past user inputs that may have already computed results (to avoid computation as this is expensive)
    Specify debug as true for feedback and opportunity to manually specify carrier fragments for species
    '''
    if not type(userinput)==str:
        raise CustomError("Please input a reaction smiles string. Include '>>' even if no products are inputted")
    splitrxn=userinput.split('>>')
    if len(splitrxn)==1: #Only reactants specified
        specs=set(splitrxn[0].split('.'))
    else:
        specs=set(splitrxn[0].split('.')).union(set(splitrxn[1].split('.')))
    inputquery={'smiles':userinput,'species':{}}
    specs={Chem.MolToSmiles(molfromsmiles(spec)) for spec in specs}
    species=inputquery['species']  
    for spec in specs:
        species.update({spec:{}})
        if refquery and spec in refquery['species']:
            carrierfrags=list(refquery['species'][spec].keys())
        else:
            carrierfrags2=Counter(getCarrierFrags0(spec,expand=expand,resFormat='smiles'))
            carrierfrags=list(carrierfrags2.keys())
        ### User input ####
        userinput2='N'
        userexpand=expand
        usermodified=False
        if debug:
            while userinput2=='N':
                userinput2=input("Identified carrier fragments for species "+str(spec)+" are: "+str(carrierfrags)+'. Input Y/N: ')
                if userinput2=='N':
                    choice=input("Do you want to specify expansion number (A) or carrier fragments directly (B): ")
                    if choice=='A':
                        usermodified=False
                        userexpand=input('Specify expansion number:')
                        userexpand=int(userexpand)
                        carrierfrags2=Counter(getCarrierFrags0(spec,expand=userexpand,resFormat='smiles'))
                        carrierfrags=list(carrierfrags2.keys())
                    elif choice=='B':
                        usermodified=True
                        userexpand=expand
                        fragvalid=False
                        iters=0
                        carrierfrags2=Counter(carrierfrags)
                        while not fragvalid: 
                            iters+=1
                            if iters>5:
                                raise CustomError('You have exceeded the maximum tries in inputting valid fragments.')
                            inputstr=input("Please supply valid carrier fragments instead of "+str(carrierfrags)+': ')
                            carrierfrags=inputstr.replace(' ','').split(',')
                            for carrierfrag in carrierfrags:
                                matches=findfragsub(spec,carrierfrag,fragment=False,addHs=True,returnindices=True)
                                if matches:
                                    count=len(matches)
                                    carrierfrags2[carrierfrag]=count
                                    fragvalid=True
                                else:
                                    fragvalid=False
                        
        species[spec].update({carrierfrag:{'parent':spec,'count':carrierfrags2[carrierfrag],'expand':userexpand,'usermodified':usermodified} for carrierfrag in carrierfrags})
    return inputquery

def getanaloguespecies(inputquery,DBsource,SQL=False,refquery={},ncpus=16,fragtable=None,substancedbsource=None,includeisotopes=False):
    '''
    Takes in input query dictionary and fragment database (DBsource) and populates each carrier fragment, returning updated dictionary
    with analogue species pools. Also returns a fragment dictionary (keys are fragments and values are query species).
    
    Pass in an SQL connection under DBsource and put SQL as true if memory is low, or pass in a file or database location address
    refquery refers to past user inputs that may have already computed results (to avoid computation as this is expensive)
    ncpus indicates how many CPUs or cores for parallel execution
    fragtable is a dataframe of possible carrier fragments in the current database, and is needed if the user specifies 
    custom fragments (as results are not prestored)
    
    substancedbsource can either be a substance dataframe or address, mainly for fragments containing only carbon and
    hydrogen (aromatic/aliphatic hydrocarbons)as the fragment identification algorithm does not extract these 
    (Only heteroatoms)
    
    Specify includeisotopes as True if isotope fragments need to be considered
    
    '''
    if type(DBsource)==str:
        if not SQL: #Address/location specified
            DB=pd.read_pickle(DBsource)
        else:
            DB=sqlite3.connect(DBsource)
    elif type(DBsource)==pd.core.frame.DataFrame:
        DB=DBsource
    elif type(DBsource)==sqlite3.Connection:
        DB=DBsource
    inputquery2=copy.deepcopy(inputquery) #Output populated dictionary
    fragdict={} #To document fragments that are completed
    for spec in inputquery['species']:
        for frag in inputquery['species'][spec]:
#             breakpoint()
            fraginfo=inputquery2['species'][spec][frag]
            if frag in fragdict: #If fragment has already been processed previously
                refspec=fragdict[frag][0]
                analoguepool=inputquery2['species'][refspec][frag]['analoguepool']
                fragdict[frag].extend([fraginfo['parent']])
            else:
                if frag not in fragdict:
                    fragdict.update({frag:[fraginfo['parent']]})
                if refquery and spec in refquery['species'] and frag in refquery['species'][spec] and 'analoguepool' in refquery['species'][spec][frag]: #If user provides a reference query
                    analoguepool=refquery['species'][spec][frag]['analoguepool']
                else:
                    #Check if fragment is in database
                    if SQL:
                        sql3='''SELECT Exists(SELECT 1 from FragmentDB1 Where FragmentSmiles="'''+ fragmentsmiles + '''")'''
                        result=pd.read_sql_query(sql3,DB)
                        present=result.iloc[0].values[0]
                    else:
                        if frag in fragdb.index:
                            present=True
                        else:
                            present=False
                            
                    if present:
                        if includeisotopes: #Other fragments that include specified fragment as substructure
                            if fragtable is None:
                                raise CustomError('Please supply a fragment table or dataframe of fragments')
                            fraglist=[]
                            initray(num_cpus=ncpus)
                            freqtabledis=mpd.DataFrame(freqtable)
                            fraglist=freqtabledis.apply(ffsrow,patt=frag,colname='FragSmiles',axis=1,result_type='reduce')
                            fraglist=pd.Series(data=fraglist.values,index=fraglist.index)
                            fraglist=list(fraglist[fraglist.values==True].index)
                            analoguepool=getCompPool(DB,fraglist,SQL=SQL)
                            if not analoguepool.empty:
                                analoguepool=analoguepool.droplevel(0)
                                aggreg=analoguepool.groupby(level='SubstanceID')['count'].sum()
                                analoguepool=analoguepool[~analoguepool.index.duplicated(keep='first')]
                                analoguepool[['count']]=aggreg
                        else:
                            analoguepool=getCompPool(DB,frag,SQL=SQL)
                    else: #User modified
                        if substancedbsource is None:
                            raise CustomError('Please supply a substance dataframe or address for custom fragments')
                        elif type(substancedbsource)==str:
                            substancedb=pd.read_pickle(substancedbsource)
                        else:
                            substancedb=substancedbsource
                        initray(num_cpus=ncpus)
                        substancedbdis=mpd.DataFrame(substancedb)
                        matches=substancedbdis.apply(ffsrow,patt=frag,colname='Smiles',axis=1,result_type='reduce')
                        matches=pd.Series(data=matches.values,index=matches.index)
                        match=list(matches[matches.values==True].index) 
                        analoguepool=substancedb.loc[substancedb.index.isin(match)]
            fraginfo.update({'analoguepool':analoguepool})     
    return inputquery2,fragdict


           
def updatequery(inputquery,fragchoice={},similarity=True,fingerprint='morgan',morganradius=2,addHs=True,
                molwt=True,refquery={},ncpus=16):
    '''
    Takes an input query dictionary with analogue species pools and adds additional information such as 
    fingerprint similarity and molecular weight
    
    Specify fragments under fragchoice if user wants to selectively apply the update (or leave blank if intent is to apply to all)
    Specify similarity as true for calculation of fingerprint similarity (fingerprint user-defined as morgan or topological,
    with user-defined morgan radius, addHs reflects whether hydrogens are considered or not)
    Specify molwt as true if molecular weights need to be calculated
    refquery refers to past user inputs that may have already computed results (to avoid computation as this is expensive)
    ncpus indicates how many CPUs or cores for parallel execution
    
    '''
    species=inputquery['species']
    for spec in species:
        for frag in species[spec]:
            if fragchoice and frag not in fragchoice:
                continue
            error=False
            fraginfo=species[spec][frag]
            if refquery:
                try:
                    analoguepool=refquery['species'][spec][frag]['analoguepool']
                    if similarity:
                        if 'similarity' in analoguepool.dtypes:
                            fraginfo['analoguepool']=analoguepool
                        else:
                            raise Exception
                    if molwt:
                        if 'molwt' in analoguepool.dtypes:
                            fraginfo['analoguepool']=analoguepool
                        else:
                            raise Exception
                except Exception:
                    error=True
                    pass
            elif not refquery or error:
                analoguepool=fraginfo['analoguepool']
                initray(num_cpus=ncpus)
                analoguepooldis=mpd.DataFrame(analoguepool)
                if similarity: 
                    sim=analoguepooldis.apply(getSimilarityrow,queryspec=spec,fingerprint=fingerprint,
                                             morganradius=morganradius,addHs=addHs,axis=1,result_type='reduce')
                    if fingerprint=='morgan':
                        analoguepool['relevance'+'_'+'morgan'+'_'+str(morganradius)]=pd.DataFrame(data=sim.values,index=sim.index,columns=['relevance'+'_'+'morgan'+'_'+str(morganradius)])
                    else:
                        analoguepool['relevance_top']=pd.DataFrame(data=sim.values,index=sim.index,columns=['relevance_top'])
                if molwt:
                    mw=analoguepooldis.apply(MWrow,axis=1,result_type='reduce')
                    analoguepool['molwt']=pd.DataFrame(data=mw.values,index=mw.index,columns=['molwt'])
            fraginfo['analoguepool']=analoguepool
    return inputquery
                    

def getcombinedpool(inputquery,fragchoice={},ST=None,fingerprint='morgan',morganradius=2,MWT=None,
                    nomixtures=True,res_format='list'):
    '''
    Retrieves combined analogue pool across all fragments. Additional filters can be applied. Output is either dataframe with
    extensive information (compound ID, SMILES, count, fragment, and query species) or a list of analogue IDs 
    depending on res_format (df and list respectively)
    
    Leave fragchoice empty if all fragments need to be considered or specify which ones need to be retrieved
    Specify nomixtures as True if mixtures need to be filtered out
    Specify similarity (lower bound, based on fingerprint and morgan radius if morgan) and molecular weight threshold 
    (upper bound) under ST and MWT respectively

    '''
    species=inputquery['species']
    combinedpooldf=[]
    for spec in species:
        for frag in species[spec]:
            if fragchoice and frag not in fragchoice:
                continue
            fraginfo=species[spec][frag]
            analoguepool=fraginfo['analoguepool']
            if nomixtures: #No mixtures desired
                analoguepool=analoguepool.loc[analoguepool['>1 Compound']==False]
            if ST is not None: #User specified similarity theshold
                ST=float(ST)
                if fingerprint=='morgan':
                    colname='relevance'+'_'+'morgan'+'_'+str(morganradius)
                else:
                    colname='relevance_top'
                if colname not in analoguepool.dtypes:
                    raise CustomError('Similarity threshold specified but no similarity information in data')
                analoguepool=analoguepool.loc[(analoguepool[colname]!='Error') & (analoguepool[colname]>=ST)]
            if MWT is not None:
                MWT=float(MWT)
                if 'molwt' not in analoguepool.dtypes:
                    raise CustomError('Molecular weight threshold specified but no molecular weight information in data')
                analoguepool=analoguepool.loc[(analoguepool['molwt']!='Error') & (analoguepool['molwt']<=MWT)]
            if res_format=='df': #User requires a dataframe
                analoguepool['FragmentSmiles']=[frag]*len(analoguepool)
                analoguepool['queryspecies']=[spec]*len(analoguepool)
            combinedpooldf.append(analoguepool)
    combinedpooldf2=pd.concat(combinedpooldf)
    if res_format=='list':
        combinedpool=set(combinedpooldf2.index)
        return combinedpool
    else:
        return combinedpooldf2
    
def getCompPool(DB,fragmentsmiles,SQL=False): #For sql to work, must specify db connection under DB. SQL can be very fast or very slow (no consistency)
    '''
    Retrieve analogue compounds for a given fragment SMILES (or list of fragment SMILES), given a database or SQL connection (DB)
    '''
#     DB=sqlite3.connect("/home/aa2133/Impurity-Project/Reaxys_Data/SQL/Reaxys_Data.db")
#     cursor = DB.cursor()
    if type(fragmentsmiles)==list:
        if not SQL:
            return DB.loc[DB.index.isin(fragmentsmiles,level=0)]     
        else:
            fragmentsmiles=['''"'''+frag+'''"''' for frag in fragmentsmiles]
            sql3='''SELECT FragmentSmiles,SubstanceID,Smiles,">1 Compound",count from FragmentDB1 Where FragmentSmiles In '''+'('+','.join(fragmentsmiles)+')'
            dat=pd.read_sql_query(sql3,DB)
            if dat.empty:
                raise CustomError('Fragment does not exist in database')
            dat.set_index(['FragmentSmiles','SubstanceID'],inplace=True)
            return dat
    else:
        if not SQL:
            return DB.xs(fragmentsmiles)
        else:
            sql3='''SELECT SubstanceID,Smiles,">1 Compound",count from FragmentDB1 Where FragmentSmiles= "'''+ fragmentsmiles + '''"'''
            dat=pd.read_sql_query(sql3,DB)
            if dat.empty:
                raise CustomError('Fragment does not exist in database')
            dat.set_index('SubstanceID',inplace=True)
            return dat 
        
def findfragsub(parent,patt,fragment=False,addHs=True,returnindices=False):
    '''
    Given a parent and pattern (patt) smiles, calculates if there is a substructure match or not and, optionally, returns
    atom index matches. Can also process mixtures, and will return a dictionary containing atom indices. 
    
    Fragments are always passed in as the pattern, but the parent can either be another fragment (fragment = True) 
    or a molecule (fragment = False).
    The addHs option is mainly to account for hydrogens but may not work if a fragment is passed in due to 
    valence/sanitization errors. In these cases, False is returned.
    Specify returnindices as True if a list of matches needs to be retrieved for fragment with atom indices
    
    
    
    '''
#     breakpoint()
    if '.' in parent or '.' in patt: #mixture present
        if '.' in parent:
            parent=parent.split('.')
        else:
            parent=[parent]
        if '.' in patt:
            patt=patt.split('.')
        else:
            patt=[patt]
        patdict={}
        for pat in patt:
            match=False
            for par in parent:
                substructmatch=findfrag(par,pat,fragment=fragment,addHs=addHs,returnindices=returnindices)
                if substructmatch:
                    match=True
                    if returnindices:
                        if pat not in patdict:
                            patdict.update({pat:{'Smiles':[par],'Indexmatch':[substructmatch]}})
                        else:
                            patdict[pat]['Smiles'].extend([par])
                            patdict[pat]['Indexmatch'].extend([substructmatch])
                    else:
                        break
        if not match:
            if returnindices:
                return {}
            else:
                return False
        else:
            if returnindices:
                return patdict
            else:
                return True
    else:
        return findfrag(parent,patt,fragment=fragment,addHs=addHs,returnindices=returnindices)
            
                
    
def findfrag(parent,patt,fragment=False,addHs=True,returnindices=False): # [HH] is not a valid SMARTS but is a valid SMILES
    '''
    Given a parent and pattern (patt) smiles, calculates if there is a substructure match or not and, optionally, returns
    atom index matches
    
    Fragments are always passed in as the pattern, but the parent can either be another fragment (fragment = True) 
    or a molecule (fragment = False).
    The addHs option is mainly to account for hydrogens but may not work if a fragment is passed in due to 
    valence/sanitization errors. In these cases, False is returned.
    Specify returnindices as True if a list of matches needs to be retrieved for fragment with atom indices

    '''
    try:
        if fragment:
            parentmol=Chem.MolFromSmarts(parent)
        else:
            parentmol=molfromsmiles(parent)
            if addHs:
                parentmol=Chem.AddHs(parentmol)
        pattmol=Chem.MolFromSmarts(patt) #Pattern 
        if not returnindices:
            return parentmol.HasSubstructMatch(pattmol)
        else:
            return parentmol.GetSubstructMatches(pattmol)
    except Exception:
        if not returnindices:
            return False 
        else:
            return tuple()

def ffsrow(row,patt,colname='FragSmiles',reverse=False,returnindices=False):
    '''
    Given a fragment dataframe, applies findfragsub across each row (parallelisation)
    colname is either 'FragSmiles' if main species is a fragment or 'Smiles' if main species is whole molecule
    Specify reverse as True if a fragment needs to be searched for a substructure parent match instead of the other way 
    round
    Specify returnindices as True if a list of matches needs to be retrieved for fragment with atom indices
    '''
    if colname=='FragSmiles':
        addHs=False
        fragment=True
        parent=row.name
    elif colname=='Smiles':
        fragment=False
        addHs=True
        parent=row['Smiles']
    if parent=='Error':
        return False
    if not reverse:
        return findfragsub(parent,patt,fragment=fragment,addHs=addHs,returnindices=returnindices)
    else:
        return findfragsub(patt,parent,fragment=fragment,addHs=addHs,returnindices=returnindices)
    

        
def getSimilarity(spec1,spec2,fingerprint='morgan',morganradius=2,addHs=True):
    '''
    Given SMILES of two species, and a fingerprint type (morgan or topological), computes fingerprint similarity.
    For morgan, a radius can be specified (default is 2). Note that fragments do not work as valence/structure errors will be thrown out by RDKit
    '''
    mol1=molfromsmiles(spec1)
    mol2=molfromsmiles(spec2)
    if addHs:
        mol1=Chem.AddHs(mol1)
        mol2=Chem.AddHs(mol2)
    if fingerprint=='morgan':
        fp1 = AllChem.GetMorganFingerprint(mol1,morganradius)
        fp2 = AllChem.GetMorganFingerprint(mol2,morganradius)
        return DataStructs.DiceSimilarity(fp1,fp2)
    elif fingerprint=='topological':
        fp1=FingerprintMols.FingerprintMol(mol1)
        fp2=FingerprintMols.FingerprintMol(mol2)
        return DataStructs.FingerprintSimilarity(fp1,fp2)

def getSimilarityrow(row,queryspec,fingerprint='morgan',morganradius=2,addHs=True):
    '''
    Calculates similarity per row of a dataframe, returns error if there are issues in representing the SMILES
    
    '''
    analoguespec=row['Smiles']
    try:
        return getSimilarity(analoguespec,queryspec,fingerprint=fingerprint,morganradius=morganradius,addHs=addHs)
    except Exception:
        return 'Error'
    
def MWspec(spec,smiles=True,mol=False):
    '''
    Calculates molecular weight of species, given in
    SMILES or mol form
    '''
    if smiles:
        molwt=MolWt(molfromsmiles(spec))
    elif mol:
        molwt=MolWt(spec)
    return molwt

def MWrow(row):
    '''
    Calculates molecular weight of species per row in a
    dataframe, returns error if there are issues in representing the SMILES
    
    '''
    smiles=row['Smiles']
    try:
        return MWspec(smiles)
    except Exception:
        return 'Error'

def updatecombinedpool(combinedpool,exemptionlist=[],catalyst=[],DBsource=None,SQL=True):
    '''
    Updates combined analogue pool with other species exempted from filtering, and returns new pool
    
    exemptionlist refers to a list of compounds that are exempted from the search (eg. catalysts);
    catalyst is user-specified and is used to find matches in the database, which are also exempted;
    DBsource (optional) refers to the substance database, used if a catalyst is specified (either file address or dataframe)
    Pass in an SQL connection under DBsource and put SQL as true if memory is low
    
    '''
    if exemptionlist:
        combinedpoolex=combinedpool.union(set(exemptionlist))
    if catalyst:
        if DBsource is None:
            raise CustomError('Catalyst specified but no reference database indicated under DBsource for checking')
        elif type(DBsource)==str:
            if not SQL: #Address/location specified
                DB=pd.read_pickle(DBsource)
            else:
                DB=sqlite3.connect(DBsource)
        elif type(DBsource)==pd.core.frame.DataFrame:
            DB=DBsource
        elif type(DBsource)==sqlite3.Connection:
            DB=DBsource
        if type(catalyst)!=list:
            catalyst=[catalyst]
        catidx=[]
        for cat in catalyst:
            try:
                catsmiles=Chem.MolToSmiles(Chem.MolFromSmiles(cat))
            except Exception:
                continue
            sql3='''SELECT SubstanceID from SubstanceDB Where Smiles= "'''+ catsmiles + '''"'''
            dat=pd.read_sql_query(sql3,db)
            dat.set_index('SubstanceID',inplace=True)
            catidx+=list(dat.index)
        combinedpoolex=combinedpoolex.union(set(catidx))
    return combinedpoolex
    

# %load ./AnalgCompds.py

#%% This section is to define classes (In progress --for now focusing on getting the workflow to work..Note none of the
# classes defined auto-update attributes when any are changed. The instance must be reinitialized. Refer to 
# https://gist.github.com/jessicald/2861038 for a comparison of dicts and classes

from rdkit import Chem
from rdkit.Chem.Descriptors import MolWt
from MainFunctions import molfromsmiles,CustomError,mol_with_atom_index
import sqlite3
import pandas as pd
from FunctionsDB import getmixturefrags,getCarrierFrags,mixtures,joindf,getfragsmiles
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit import DataStructs
import ray
import modin.pandas as mpd
import dask.dataframe as dd

#%% Parallel functions

def getCompPool(DB,fragmentsmiles,SQL=False): #For sql to work, must specify db connection under DB. SQL can be very fast or very slow (no consistency)
    if not SQL:
        return DB.xs(fragmentsmiles)
    else:
        sql3='''SELECT SubstanceID,Smiles,">1 Compound",FragmentSmarts from FragmentDB Where FragmentSmiles= "'''+ fragmentsmiles + '''"'''
        dat=pd.read_sql_query(sql3,DB)
        dat.set_index('SubstanceID',inplace=True)
        return dat 
def MWCompd(series,rows=True, Mol=False,Smiles=False):
    if rows:
        try:
            return MolWt(molfromsmiles(series['Smiles']))
        except Exception:
            return 'Error'
    elif Smiles:
        return MolWt(molfromsmiles(series))
    elif Mol:
        return MolWt(series)
        
def MW(partition):
    return partition.apply(MWCompd,axis=1)

def getSimilarity(series,target,rows=True, CompdMol=False, CompdSmiles=False, CompdFingerprint=False, 
                  TargetFingerprint=True, TargetSmiles=False,
                  TargetMol=False, UseCompID=False,DB=None):
    # list of IDs, str ID, thresh, type of targetComp-> list of IDs
    # get similar compounds of target compound based on tanimoto
    # compPool: list of compounds to compare with targetComp
    # targetComp: str ID of the targetComp
    # UseCompID: bool, True -> targetComp is a reaxys ID, False will be smiles str
    if not TargetFingerprint:
        if UseCompID:
            if DB is None:
                return "Please include a fragment database"
            targetMol = molfromsmiles(DB.xs(target,level=1)['Smiles'].unique()[0])
        elif TargetSmiles: # should be a smiles str
            targetMol = molfromsmiles(target)
        else:
            targetMol=target
        target_fp = FingerprintMols.FingerprintMol(targetMol)
    else:
        target_fp=target
    
    if rows:
        try:
            return DataStructs.FingerprintSimilarity(target_fp,FingerprintMols.FingerprintMol(molfromsmiles(series['Smiles'])))
        except Exception:
            return 'Error'
    elif CompdSmiles:
        return DataStructs.FingerprintSimilarity(target_fp,FingerprintMols.FingerprintMol(molfromsmiles(series)))
    elif CompdMol:
        return DataStructs.FingerprintSimilarity(target_fp,FingerprintMols.FingerprintMol(series))
    elif CompdFingerprint:
        return DataStructs.FingerprintSimilarity(target_fp,series)
    

def getSim(partition,target):
    return partition.apply(getSimilarity,target=target,axis=1)
        
        
        
#%% Classes
        
class queryrxn:
    def __init__(self,inputsmiles=None,species=None,fragsize=6,client=None,FragDB=None,SQL=False,modinray=True,MWThresh=300,SimThresh=0.5):
        '''
        input smiles: Smiles of reaction
        species: All species involved in the reaction (should be either compound or querycompd class)
        fragsize: Only fragsize=6 works at the moment as FragDB has fragments of this size
        client: Distributed client for parallel procesing. If defined, analogue compound lists will be further filtered
        based on molecular weight and similarity
        FragDB: Fragment database for grouping analogue compounds containing the same fragment
        MWThresh: Molecular weight threshold for filtering analogue compound lists
        SimThresh: Similarity threshold for filtering analogue compound lists

        Parameters
        ----------
        inputsmiles : TYPE, optional
            DESCRIPTION. The default is None.
        species : TYPE, optional
            DESCRIPTION. The default is None.
        fragsize : TYPE, optional
            DESCRIPTION. The default is 6.
        client : TYPE, optional
            DESCRIPTION. The default is None.
        FragDB : TYPE, optional
            DESCRIPTION. The default is None.
        SQL : TYPE, optional
            DESCRIPTION. The default is False.
        modinray : TYPE, optional
            DESCRIPTION. The default is True.
        MWThresh : TYPE, optional
            DESCRIPTION. The default is 300.
        SimThresh : TYPE, optional
            DESCRIPTION. The default is 0.5.

        Raises
        ------
        CustomError
            DESCRIPTION.

        Returns
        -------
        None.

        '''

        if inputsmiles is not None:
            self.smiles=inputsmiles
            if not type(inputsmiles)==str:
                raise CustomError("Please input a reaction smiles string. Include '>>' even if no products are inputted")
            self.species=self.SpeciesList(fragsize=fragsize,client=client,FragDB=FragDB,SQL=SQL,modinray=modinray,MWThresh=MWThresh,SimThresh=SimThresh)
        elif species is not None:
            self.species=species
        self.speciesdict=self.SpeciesDict(self.species)
            
    def SpeciesList(self,fragsize=6,client=None,FragDB=None,SQL=False,modinray=True,MWThresh=300,SimThresh=0.5):
        specieslist=[]
        splitrxn=self.smiles.split('>>')
        if len(splitrxn)==1: #Only reactants specified
            rcts=splitrxn[0].split('.')
            prods=[]
        else:
            rcts=set(splitrxn[0].split('.'))
            prods=set(splitrxn[1].split('.'))
        specieslist+=[querycompd(Type='rct',smiles=rct,verifysmiles=True,fragsize=fragsize,client=client,FragDB=FragDB,SQL=SQL,
                    modinray=modinray,MWThresh=MWThresh,SimThresh=SimThresh) for rct in rcts]
        specieslist+=[querycompd(Type='prod',smiles=prod,verifysmiles=True,fragsize=fragsize,client=client,FragDB=FragDB,SQL=SQL,
                    modinray=modinray,MWThresh=MWThresh,SimThresh=SimThresh) for prod in prods]
        return specieslist

    def SpeciesDict(self,specieslist):
        return {species.smiles:species for species in specieslist}
        

class compound:
    '''
    General compound class to store important information such as smiles, smarts, mol, fingerprint, type
    (reactant, reagent or product)
    
    '''
    def __init__(self,Type=None,smiles=None,smarts=None,verifysmiles=False,verifysmarts=False, fingerprint=None,
                formula=None):
        if Type:
            self.type=Type
        if smiles:
            if verifysmiles:
                self.smiles=Chem.MolToSmiles(molfromsmiles(smiles))
            else:
                self.smiles=smiles
            self.mol=self.MolFromSmiles()
            self.smarts=Chem.MolToSmarts(self.mol)
        if smarts:
            if verifysmarts:
                self.smarts=Chem.MolToSmarts(Chem.MolFromSmarts(smarts))
            else:
                self.smarts=smarts
            self.mol=self.MolFromSmarts()
            self.smiles=Chem.MolToSmiles(self.mol)

        if fingerprint:
            self.fingerprint=fingerprint
        else:
            self.fingerprint=self.set_fingerprint()
        
        if formula:
            self.formula=formula
        else:
            self.formula=self.Formula()
    
    def MolFromSmiles(self):
        '''
        Converts a smiles string into a mol object for RDKit use. Sanitizes and makes sure 
        molecule is valid. This is not stored as an attribute as it will take up excessive
        memory especially given there are 16 million compounds stored in Reaxys.

        '''
        return molfromsmiles(self.smiles)
    
    def MolFromSmarts(self):
        '''
        Converts smarts string into a mol object for RDKit use. Sanitizes and makes sure 
        molecule is valid.
        
        '''
        mol=Chem.MolFromSmarts(self.smarts)
        try:
            Chem.SanitizeMol(mol)
        except Exception:
            pass
        mol.UpdatePropertyCache(strict=False)        
        return mol
    
    def atom_index(self):
        '''
        Draw a molecule with atom index numbers (set by RDKit)

        '''
        return mol_with_atom_index(self.mol)
    
    def set_fingerprint(self):
        
        return FingerprintMols.FingerprintMol(self.mol)
    
    def loc_record(self,FragDB):
        '''
        Locates and returns record of compound in the Reaxys database (stored on the server)
        No point in retrieving record for one compound (takes around 4 seconds). Might as well
        recompute everything in less than a second
        '''
        if FragDB.index.name=='Smiles':
            return FragDB.loc[self.smiles]
        else:
            return FragDB.loc[FragDB['Smiles']==self.smiles]
        
    def Mixtures(self):
        '''
        Returns True if the compound is a mixture and False if not
        '''
        return mixtures(self.smiles)
    
    def MolWt(self):
        return MWCompd(self.mol,rows=False,Mol=True)
    
    def Formula(self):
        return Chem.rdMolDescriptors.CalcMolFormula(self.mol)
    

class querycompd(compound):
    '''
    Query compound class built from user input. Has all attributes of the compound class, and also contains relevant
    fragment information.
    
    '''
    def __init__(self,count=1,fragments=None,fragsize=6,client=None,FragDB=None,SQL=False,modinray=True,MWThresh=300,SimThresh=0.5,**kwargs):
        super(querycompd,self).__init__(**kwargs)
        self.count=count
        if fragments is not None:
            self.frags=fragments
        else:
            self.frags=self.FragsList(fragsize,FragDB=FragDB,SQL=SQL,client=client,modinray=modinray,MWThresh=MWThresh,SimThresh=SimThresh)
        self.fragsdict=self.FragsDict(self.frags)
            
    def Frags(self, fragsize, resFormat='smarts', addHs=True):
        '''
        Returns fragments of a specified size
        '''
        if self.Mixtures():
            return getmixturefrags(self.smiles,fragsize,resFormat=resFormat,addHs=addHs)
        else:
            return getCarrierFrags(self.smiles,fragsize,resFormat=resFormat,addHs=addHs)
    
    def FragsList(self,fragsize,FragDB=None,SQL=False,client=None,modinray=True,MWThresh=300,SimThresh=0.5):
        '''
        Builds a list of fragment objects. If there are repeated fragments, keeps a count
        
        '''
        fragdict={}
        frags=self.Frags(fragsize)
        if type(frags)!=list:
            frags=[frags]
        for frag in frags:
            fragsmiles=getfragsmiles(frag)
            if fragsmiles in fragdict.keys():
                fragdict[fragsmiles][1]+=1
            else:
                fragdict.update({fragsmiles:[frag,1]})
        return [fragment(parent_instance=self,count=fragdict[fragsmiles][1],smarts=fragdict[fragsmiles][0],FragDB=FragDB,SQL=SQL,
                         client=client,modinray=modinray,MWThresh=MWThresh,SimThresh=SimThresh) for fragsmiles in fragdict]
        
    def FragsDict(self,fragslist):
        return {frag.smiles:frag for frag in fragslist}

class fragment(compound):
    '''
    Fragment class. Objects/instances are fragments containing all attributes of the compound class. Important
    attributes include CompPool, which reflects all compounds in Reaxys containing the fragment.
    
    __init__:
    
    FragDB: Fragment database, which needs to be correct with relation to size of fragment belonging to query compound.
    
    parent_instance: Query compound (querycompd class) to which the fragment originally belongs
    
    parent_smiles: Query compound smiles. Either above or this can be specified, thus allowing
    the user to decide how the fragment instance should be initialized (In the latter case, the smiles should be valid)
    
    MWThresh:
    
    SimThresh:
    
    client: If the user specifies a client (cluster), parallel computing can be used. Otherwise sequential will be used
    (not advisable as this is much slower). Refer to dask section to specify.
    
    **kwargs: 
    
    Attributes:
    
    '''
    def __init__(self,count=1,parent_instance=None,parent_smiles=None,FragDB=None,SQL=False,client=None,modinray=True,
                 MWThresh=300,SimThresh=0.5,**kwargs):
        
        super(fragment,self).__init__(**kwargs)
        self.count=count
        if parent_smiles and parent_instance is None:
            parent_instance=querycompd(smiles=parent_smiles)
        if parent_instance is not None:
            self.parent=parent_instance
        if FragDB is not None:
            self.comppool=self.CompPool(FragDB,SQL=SQL,parent_instance=parent_instance,client=client,modinray=modinray)
            self.filteredpoolMW=self.FilterPool(MWThresh=MWThresh)
            if parent_instance is not None:
                self.filteredpoolSim=self.FilterPool(SimThresh=SimThresh)
                self.filteredpool=self.FilterPool(MWThresh=MWThresh,SimThresh=SimThresh)
        
        
    def CompPool(self,FragDB,SQL=False,parent_instance=None,client=None,modinray=True):
        '''
        Retrieves general analogue compound pool containing all compounds that contain the fragment
        
        '''

        comppool=getCompPool(FragDB,self.smiles,SQL=SQL) #Database needs to be multiindexed at fragment smiles and substance ID
        comppool['Count']=comppool.groupby([comppool.index])['Smiles'].transform('count')
        comppool=comppool[~comppool.index.duplicated(keep='first')]
        if modinray: #This means modin+ ray is chosen
            pool=mpd.DataFrame(comppool)
            molwt=pool.apply(MWCompd,axis=1,result_type='reduce')
            comppool=joindf(pd.DataFrame(data=molwt.values,index=molwt.index,columns=['MolWt']),comppool)
            if parent_instance:
                fp_similarity=pool.apply(getSimilarity,target=parent_instance.fingerprint,axis=1,result_type='reduce')
                comppool=joindf(pd.DataFrame(data=fp_similarity.values,index=fp_similarity.index,columns=['Similarity']),comppool)
            comppool=comppool[comppool.MolWt!='Error']
        elif client is not None: #This means parallel processing dask is chosen
            DaskPool=client.persist(dd.from_pandas(comppool,npartitions=16))
            molwt=DaskPool.map_partitions(MW,meta=('MolWt','int64')).compute()
            comppool=joindf(pd.DataFrame(molwt),comppool)
            if parent_instance:
                fp_similarity=DaskPool.map_partitions(getSim,target=parent_instance.fingerprint,
                                                      meta=('Similarity','int64')).compute()
                comppool=joindf(pd.DataFrame(fp_similarity),comppool)
            comppool=comppool[comppool.MolWt!='Error']
            del DaskPool

        return comppool
                
            
    
    def FilterPool(self,MWThresh=None,SimThresh=None):
        if MWThresh and not SimThresh:
            return self.comppool.loc[self.comppool['MolWt']<=MWThresh]
        elif SimThresh and not MWThresh:
            return self.comppool.loc[self.comppool['Similarity']>=SimThresh]
        elif MWThresh and SimThresh:
            return self.comppool.loc[(self.comppool['Similarity']>=SimThresh) & (self.comppool['MolWt']<=MWThresh)] 
            

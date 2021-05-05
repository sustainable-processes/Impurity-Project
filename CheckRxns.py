#%load ./CheckRnxs.py

from collections import OrderedDict #Preserve order of keys relative to reaction from left to right
from rdkit import Chem
from rdkit.Chem import rdChemReactions #Reaction processing
import copy
from FunctionsDB import getCarrierFrags

def checkrxnrow(row):
    return checkrxn(mappedrxn=row['mapped_rxn'],Rdata=row['LHSdata'],Pdata=row['RHSdata'])


def checkrxn(mappedrxn,LHSdata=OrderedDict({}),RHSdata=OrderedDict({}),Rdata=None,Pdata=None): #Assume same rxn smiles stored next to each other
    LHSdata=OrderedDict({})
    RHSdata=OrderedDict({})
    rdrxn=rdChemReactions.ReactionFromSmarts(mappedrxn,useSmiles=True)
    cleanrxn=copy.copy(rdrxn)
    rdChemReactions.RemoveMappingNumbersFromReactions(cleanrxn)
    rmixtures=[ID for ID in Rdata if len(Rdata[ID]['smiles'].split('.'))>1]
    pmixtures=[ID for ID in Pdata if len(Pdata[ID]['smiles'].split('.'))>1]
    msg=''
    if rmixtures:
        msg+='LHS species '+', '.join([str(ID) for ID in rmixtures])+' is a mixture'
    if pmixtures:
        if msg:
            msg=msg+', '+'RHS species '+', '.join([str(ID) for ID in pmixtures])+' is a mixture'
        else:
            msg+='RHS species '+', '.join([str(ID) for ID in pmixtures])+' is a mixture'
    if msg:
        return LHSdata,RHSdata,msg
    else:
        msg+='Valid'

    for ID,rct in enumerate(cleanrxn.GetReactants()):
        rctsmiles=Chem.MolToSmiles(rct)
        for ID0 in Rdata:
            if Rdata[ID0]['smiles']==rctsmiles:
                if ID0 in LHSdata.keys():
                    LHSdata[ID0]['mappedsmiles'].extend([Chem.MolToSmiles(rdrxn.GetReactants()[ID])])
                    LHSdata[ID0]['cleanmol'].extend([rct])
                else:
                    LHSdata.update({ID0:Rdata[ID0]})
                    LHSdata[ID0].update({'mappedsmiles':[Chem.MolToSmiles(rdrxn.GetReactants()[ID])]})
                    LHSdata[ID0].update({'cleanmol':[rct]})
                break
    for ID,prod in enumerate(cleanrxn.GetProducts()):
        prodsmiles=Chem.MolToSmiles(prod)
        for ID0 in Pdata:
            if Pdata[ID0]['smiles']==prodsmiles:
                if ID0 in RHSdata.keys():
                    RHSdata[ID0]['mappedsmiles'].extend([Chem.MolToSmiles(rdrxn.GetProducts()[ID])])
                    RHSdata[ID0]['cleanmol'].extend([prod])
                else:
                    RHSdata.update({ID0:Pdata[ID0]})
                    RHSdata[ID0].update({'mappedsmiles':[Chem.MolToSmiles(rdrxn.GetProducts()[ID])]})
                    RHSdata[ID0].update({'cleanmol':[prod]})
                break
    return LHSdata,RHSdata,msg


def assignfragsrow(row,querydict,natoms):
    LHSdata=row.LHSdata
    return assignfrags(LHSdata,querydict,natoms)
def assignfrags(LHSdata,querydict,natoms):
    fragdata=copy.deepcopy(LHSdata)
    nonanalogue=[]
    for rctid in LHSdata:
        analgfrags=getCarrierFrags(fragdata[rctid]['smiles'],natoms)
        if type(analgfrags)==str:
            analgfrags=[analgfrags]
        commonfrags=set([Chem.MolToSmiles(Chem.MolFromSmarts(fragsmarts)) for fragsmarts in analgfrags]).intersection(set(querydict.keys()))
        if not commonfrags:    
            nonanalogue+=[rctid]
            continue
        fragdata[rctid]['querycompds']={}
        fragdata[rctid]['querycompds'].update({fragsmiles:querydict[fragsmiles] for fragsmiles in commonfrags})
    if nonanalogue:
        msg='Species '+', '.join([str(rctid) for rctid in nonanalogue])+' not analogue'
    else:
        msg='Valid'
    return fragdata,msg
    
        

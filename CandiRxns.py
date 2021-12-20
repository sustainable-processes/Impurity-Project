# %load ./CandiRxns.py 

from decimal import Decimal, ROUND_HALF_UP
import math
import copy
from MainFunctions import getcompdict

def addfilters(queryrxn,ST=None,MWT=None,nomixtures=True):
    '''
    Adds additional filters on filtered pool. Results stored as filtered pool
    Specify nomixtures as True if mixtures need to be filtered out
    Specify similarity and molecular weight threshold udner ST and MWT respectively
    
    '''
    if ST is not None:
        for species in queryrxn.species:
            for frag in species.frags:
                frag.filteredpoolSim=frag.FilterPool(SimThresh=ST)
                if nomixtures:
                    frag.filteredpoolSim=frag.filteredpoolSim.loc[frag.filteredpoolSim['>1 Compound']==False]
                    
    if MWT is not None:
        for species in queryrxn.species:
            for frag in species.frags:
                frag.filteredpoolMW=frag.FilterPool(MWThresh=MWT)
                if nomixtures:
                    frag.filteredpoolMW=frag.filteredpoolMW.loc[frag.filteredpoolMW['>1 Compound']==False]
               
                
    if ST is not None and MWT is not None:
        for species in queryrxn.species:
            for frag in species.frags:
                frag.filteredpool=frag.FilterPool(SimThresh=ST,MWThresh=MWT)
                if nomixtures:
                    frag.filteredpool=frag.filteredpool.loc[frag.filteredpool['>1 Compound']==False]
    if ST is None and MWT is None and nomixtures:
        for species in queryrxn.species:
            for frag in species.frags:
                frag.filteredpoolSim=frag.filteredpoolSim.loc[frag.filteredpoolSim['>1 Compound']==False]
                frag.filteredpoolMW=frag.filteredpoolMW.loc[frag.filteredpoolMW['>1 Compound']==False]
                frag.filteredpool=frag.filteredpool.loc[frag.filteredpool['>1 Compound']==False]
                
def checkunresolved(row,unresolvedids,exemptionlist=[]):
    '''
    Given a row, checks if reactants, reagents and products have valid/representable smiles
    based on a set of unresolved IDs
    
    If an exemption list is given (misclassified catalysts as reagents) it will be omitted from checking for reagents
    
    
    '''
    
    rct=set(row['ReactantID'])
    rgt=set(row['ReagentID'])
    prod=set(row['ProductID'])
    rctmatches=[val in unresolvedids for val in rct]
    prodmatches=[val in unresolvedids for val in prod]
    if exemptionlist:
        rgtmatches=[val in unresolvedids and val not in exemptionlist for val in rgt]
    else:
        rgtmatches=[val in unresolvedids for val in rgt]
    if any(rgtmatches) or any(rgtmatches) or any(prodmatches):
        return True
    else:
        return False
    
                

def checkrxts(row,combinedpool,combinedpoolex=None,relevance_s='loose',relevance_m='loose'):
    '''
    Given a row, checks if either rcts, or rcts and reagents are analogue;
    Include additional exemptions (eg. misclassified catalysts) under combinedpoolex. Records with these reagent IDs 
    even though not analogue will be kept;
    relevance_s is for single reference records, relevance_m is for multiple reference records
    For either option, specify relevance as 'loosest' if only reactants need to be analogue (combinedpool);
    as 'loose' if any reagent analogue is acceptable, taking into account exemption list (combinedpoolex);
    as 'strict' if all reagents need to be analogue, taking into account exemption list (combinedpoolex);
    as 'strictest' if all record reagents need to be analogue, not taking into account exemption list (combinedpool);
    
    ''' 
    if combinedpoolex is None:
        combinedpoolex=combinedpool
    rct=set()
    rct=set(row['ReactantID'])
    numrefs=row['NumRefs']
    if numrefs==1: #Single reference record
        relevance=relevance_s
    else: #Multiple reference record
        relevance=relevance_m
    if relevance=='loose':
        if row['ReagentID']!='NaN' and row['ReagentID']:
            matches=[val in combinedpoolex for val in set(row['ReagentID'])]
            if any(matches) and rct.issubset(combinedpool):
                return True
            else:
                return False
    if relevance=='strict':
        if row['ReagentID']!='NaN' and row['ReagentID']:
            matches=[val in combinedpoolex for val in set(row['ReagentID'])]
            if all(matches) and rct.issubset(combinedpool):
                return True
            else:
                return False
    if relevance=='strictest':
        if row['ReagentID']!='NaN' and row['ReagentID']:
            matches=[val in combinedpool for val in set(row['ReagentID'])]
            if all(matches) and rct.issubset(combinedpool):
                return True
            else:
                return False 
    if rct.issubset(combinedpool):
        return True
    else:
        return False
                       
#     else: #Multiple reference record
        
#     if not rctonly:
#         if row['ReagentID']!='NaN' and row['ReagentID']:
#             if relevance=='strict':
#                 matches=[val in combinedpoolex for val in set(row['ReagentID'])]
#                 if all(matches) and rct.issubset(combinedpool):
#                     return True
#                 elif row['NumRefs']>1 and any(matches) and rct.issubset(combinedpool): #Keep more than one reference
#                     return True
#                 else:
#                     return False
#             elif relevance=='strictest':
#                 matches=[val in combinedpoolex for val in set(row['ReagentID'])]
#                 if all(matches) and rct.issubset(combinedpool):
#                     return True
#                 else:
#                     return False
#             else:
#                 matches=[val in combinedpoolex for val in set(row['ReagentID'])]
#                 if any(matches) and rct.issubset(combinedpool):
#                     return True 
#                 else:
#                     return False
#     if rct.issubset(combinedpool):
#         return True
#     else:
#         return False
                
                
# def rxp(partition,db):
#     partition.apply(rx,db=db,axis=1)

def rx(row,db,Rdata=True,Pdata=False,Rgtdata=False,Solvdata=False,database=False):
#     Rdata={ID:getcompdict(ID=ID,FragDB=FragDB)[ID] for ID in copy.copy(row.ReactantID)} #Stores atom count, type, smiles for each LHS species keyed to ID
#     breakpoint()
    if Rdata:
        col='ReactantID'
    elif Pdata:
        col='ProductID'
    elif Rgtdata:
        col='ReagentID'
    elif Solvdata:
        col='SolventID'
    dat={}
    if row[col]=='NaN' or row[col] is None or not row[col]:
        return {}
    for ID in set(copy.copy(row[col])):
        try:
            if database:
                smiles=db.loc[ID].Smiles
            else:
                smiles=pd.read_sql_query('''SELECT Smiles from SubstanceDB Where SubstanceID=  "'''+ str(ID) + '''"''',db).Smiles[0]
            dat.update(getcompdict(ID=ID,smiles=smiles))  
        except Exception:
            continue
    return dat
                
#     try:
# #         smiles=[pd.read_sql_query('''SELECT Smiles from SubstanceDB Where SubstanceID=  "'''+ str(ID) + '''"''',db).Smiles[0] for ID in copy.copy(row['ReactantID'])]
#         dat={ID:getcompdict(ID=ID,smiles=pd.read_sql_query('''SELECT Smiles from SubstanceDB Where SubstanceID=  "'''+ str(ID) + '''"''',db).Smiles[0])[ID] for ID in set(copy.copy(row[col]))}
# #         Pdata={ID:getcompdict(ID=ID,smiles=pd.read_sql_query('''SELECT Smiles from SubstanceDB Where SubstanceID=  "'''+ str(ID) + '''"''',db).Smiles[0])[ID] for ID in copy.copy(row['ProductID'])}
#     except Exception as e:
#         return str(e)
#     else:
#         return dat                
                
                
                
# def checkrxts(row,combinedpool,rctonly=True,relevance='loose'):
#     '''
#     Given a row, checks if either rcts, or rcts and reagents are analogue
#     Specify rcctonly as True if only reactant IDs need to be checked, otherwise reagent IDs will also be checked
#     Specify relevance as 'loose' if any reagent analogue is acceptable, including more references;
#     as 'strict' if all reagents need to be analogue, and for more than one reference at least (1/numrefs * reagent number)
#     *rounded up* reagents should be analogue;
#     as 'balanced' if all reagents need to be analogue, and for more than one reference at least 1 reagent should be analogue
#     as 'strictest' if all record reagents need to be analogue, including morerefs; 
    
#     '''
#     rct=set()
#     rct=set(row['ReactantID'])    
#     if not rctonly:
#         if row['ReagentID']!='NaN':
#             matches=[val in combinedpool for val in set(row['ReagentID'])]
#             if relevance=='balanced':
#                 if all(matches) and rct.issubset(combinedpool):
#                     return True
#                 elif row['NumRefs']>1 and any(matches) and rct.issubset(combinedpool): #Keep more than one reference
#                     return True
#                 else:
#                     return False
# #                 idealanaloguenum=Decimal((1/row['NumRefs'])*(len(row['ReagentID']))).to_integral_value(rounding=ROUND_HALF_UP)
# #                 if sum(matches)>=idealanaloguenum and rct.issubset(combinedpool):
# #                     return True
# #                 else:
# #                     return False
#             elif relevance=='strict':
#                 if all(matches) and rct.issubset(combinedpool):
#                     return True
#                 elif row['NumRefs']>1:
#                     idealanaloguenum=math.ceil((1/row['NumRefs'])*(len(row['ReagentID'])))
#                     if sum(matches)>=idealanaloguenum and rct.issubset(combinedpool):
#                         return True
#                     else:
#                         return False
            
#             elif relevance=='strictest': 
#                 if all(matches) and rct.issubset(combinedpool):
#                     return True
#                 else:
#                     return False
#             else:
#                 if any(matches) and rct.issubset(combinedpool):
#                     return True
#                 elif relevance=='loose':
#                     if row['NumRefs']==1: #Line was added as reactions with more than one reference are removed if reagents not in analogue compd pool (even if reagents aren't involved)
#                                        # This adds more candirxns 15230 to 25347
#                         return False 
#                 else:
#                     return False
# #             if row['NumRefs']>1:
# #                 rct=rct.union(set(row['ReagentID'][0]))
# #             else:
# #                 rct=rct.union(set(row['ReagentID']))
#     if rct.issubset(combinedpool):
#         return True
#     else:
#         return False

# %load ./CandiRxns.py 

from decimal import Decimal, ROUND_HALF_UP
import math

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
                
def checkunresolved(row,unresolvedids):
    '''
    Given a row, checks if reactants, reagents and products have valid/representable smiles
    Based on a set of unresolved IDs
    
    '''
    rct=set(row['ReactantID'])
    rgt=set(row['ReagentID'])
    if rct.issubset(unresolvedids) or rgt.issubset(unresolvedids):
        return False
    else:
        return True
    
                

def checkrxts(row,combinedpool,rctonly=True,relevance='loose',combinedpoolex=None,morerefs=True):
    '''
    Given a row, checks if either rcts, or rcts and reagents are analogue
    Specify rcctonly as True if only reactant IDs need to be checked, otherwise reagent IDs will also be checked;
    Include additional exemptions (eg. misclassified catalysts) under combinedpoolex. Records with these reagent IDs even though not analogue
    will be kept;
    Specify morerefs as True if records with more than 1 ref (combined analogues) is kept as long as 1 analogue is present
    Specify relevance as 'loose' if any reagent analogue is acceptable, taking into account exemption list;
    as 'strict' if all reagents need to be analogue, taking into account exemption list;
    as 'strictest' if all record reagents need to be analogue, not taking into account exemption list
    
    ''' 
    if combinedpoolex is None:
        combinedpoolex=combinedpool
    rct=set()
    rct=set(row['ReactantID'])    
    if not rctonly:
        if row['ReagentID']!='NaN' and row['ReagentID']:
            if relevance=='strict':
                matches=[val in combinedpoolex for val in set(row['ReagentID'])]
                if all(matches) and rct.issubset(combinedpool):
                    return True
                elif morerefs:
                    if row['NumRefs']>1 and any(matches) and rct.issubset(combinedpool): #Keep more than one reference
                        return True
                    else:
                        return False
                else:
                    return False
            elif relevance=='strictest':
                matches=[val in combinedpool for val in set(row['ReagentID'])]
                if all(matches) and rct.issubset(combinedpool):
                    return True
                else:
                    return False
            else:
                matches=[val in combinedpoolex for val in set(row['ReagentID'])]
                if any(matches) and rct.issubset(combinedpool):
                    return True 
                else:
                    return False
    if rct.issubset(combinedpool):
        return True
    else:
        return False
                
                
def rxp(partition,db):
    partition.apply(rx,db=db,axis=1)
def rx(row,db,Rdata=True,Pdata=False,Rgtdata=False):
#     Rdata={ID:getcompdict(ID=ID,FragDB=FragDB)[ID] for ID in copy.copy(row.ReactantID)} #Stores atom count, type, smiles for each LHS species keyed to ID
    if Rdata:
        col='ReactantID'
    elif Pdata:
        col='ProductID'
    elif Rgtdata:
        dat={}
        col='ReagentID'
        if row[col]=='NaN' or row[col] is None:
            return {}
        for ID in set(copy.copy(row[col])):
            try:
                dat.update(getcompdict(ID=ID,smiles=pd.read_sql_query('''SELECT Smiles from SubstanceDB Where SubstanceID=  "'''+ str(ID) + '''"''',db).Smiles[0]))  
            except Exception:
                continue
        return dat
                
    try:
#         smiles=[pd.read_sql_query('''SELECT Smiles from SubstanceDB Where SubstanceID=  "'''+ str(ID) + '''"''',db).Smiles[0] for ID in copy.copy(row['ReactantID'])]
        dat={ID:getcompdict(ID=ID,smiles=pd.read_sql_query('''SELECT Smiles from SubstanceDB Where SubstanceID=  "'''+ str(ID) + '''"''',db).Smiles[0])[ID] for ID in set(copy.copy(row[col]))}
#         Pdata={ID:getcompdict(ID=ID,smiles=pd.read_sql_query('''SELECT Smiles from SubstanceDB Where SubstanceID=  "'''+ str(ID) + '''"''',db).Smiles[0])[ID] for ID in copy.copy(row['ProductID'])}
    except Exception as e:
        return str(e)
    else:
        return dat                
                
                
                
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

# %load ./CandiRxns.py 

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
                
        
                
def checkrxts(row,combinedpool,rctonly=True,relevance='loose'):
    '''
    Given a row, checks if either rcts, or rcts and reagents are analogue
    Specify rcctonly as True if only reactant IDs need to be checked, otherwise reagent IDs will also be checked
    Specify relevance as 'loose' if any reagent analogue is acceptable, and more refs included by reactant basis;
    as 'strictest' if all record reagents need to be analogue, including morerefs; 
    as 'balanced' if algorithm is used to check if sufficient analogues are present
    
    '''
    rct=set()
    rct=set(row['ReactantID'])    
    if not rctonly:
        if row['ReagentID']!='NaN':
            matches=[val in combinedpool for val in set(row['ReagentID'])]
            if relevance=='balanced':
                idealanaloguenum=round((1/row['NumRefs'])*(len(row['ReagentID'])))
                if sum(matches)>=idealanaloguenum and rct.issubset(combinedpool):
                    return True
                else:
                    return False
            elif relevance=='strictest': 
                if all(matches) and rct.issubset(combinedpool):
                    return True
                else:
                    return False
            else:
                if any(matches) and rct.issubset(combinedpool):
                    return True
                elif relevance=='loose':
                    if row['NumRefs']==1: #Line was added as reactions with more than one reference are removed if reagents not in analogue compd pool (even if reagents aren't involved)
                                       # This adds more candirxns 15230 to 25347
                        return False 
                else:
                    return False
#             if row['NumRefs']>1:
#                 rct=rct.union(set(row['ReagentID'][0]))
#             else:
#                 rct=rct.union(set(row['ReagentID']))
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

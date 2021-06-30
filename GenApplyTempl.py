# %load ./GenApplyTempl.py
import copy,itertools
from rdkit import Chem
from MainFunctions import molfromsmiles
from rdkit.Chem import rdChemReactions
from rdkit.Chem import Descriptors

def gen_template_row(row):
    LHSdata=copy.deepcopy(row.LHSdata)
    RHSdata=copy.deepcopy(row.RHSdata)
    res=copy.deepcopy(row.mapdict)
    nofg=row.nofg
    return gen_template(LHSdata,RHSdata,res,nofg=nofg)

def gen_template(LHSdata,RHSdata,res,nofg=[]): #Don't add invalid reactions here
    reacfrag=[]
    prodfrag=[]
    farfg=[]
    unusedprod=[]
    for analoguecompd in LHSdata:
        if analoguecompd in nofg and not LHSdata[analoguecompd]['reacfrag'] : #Hydrogen will not be mapped and won't show up in reaction center 
            reacfrag+=[Chem.MolToSmarts(molfromsmiles(mappedsmiles)) for mappedsmiles in LHSdata[analoguecompd]['mappedsmiles']]
            continue
        for inst,fraginf in LHSdata[analoguecompd]['reacfrag'].items():
            fragatomidx=set()
            for frag,matchidx in fraginf.items():
                fragatomidx=fragatomidx.union({atomidx for match in matchidx for atomidx in list(LHSdata[analoguecompd]['fragloc'][inst][frag]['corrmatches'])[match]})
            reacfragcurr=[Chem.MolFragmentToSmarts(molfromsmiles(LHSdata[analoguecompd]['mappedsmiles'][inst]),fragatomidx)]
            if len(reacfragcurr[0].split('.'))>1 and analoguecompd not in farfg: # Fragments involved in the reaction are too far away, and are decomposed into more than 1 species. Invalid template.
                farfg+=[analoguecompd]
            reacfrag+=reacfragcurr
            for atomidx in fragatomidx:
#                 breakpoint()
                assigned=False
                for val in res.values():
                    if val[0]==analoguecompd and val[1]==inst and val[2]==atomidx:
#                         breakpoint()
                        assigned=True
                        prodid=val[3]
                        prodinst=val[4]
                        prodatomidx=val[5]
                        if 'rxnatomidx' not in RHSdata[prodid].keys():
                            RHSdata[prodid]['rxnatomidx']={prodinst:{prodatomidx}}
                            assigned=True
                            break
                        elif prodinst in RHSdata[prodid]['rxnatomidx']:
                            RHSdata[prodid]['rxnatomidx'][prodinst].add(prodatomidx)
                            assigned=True
                            break
                        else:
                            RHSdata[prodid]['rxnatomidx'].update({prodinst:{prodatomidx}})
                            assigned=True
                            break
        
    for prodid in RHSdata:
        if 'rxnatomidx' not in RHSdata[prodid].keys() and RHSdata[prodid]['formula']=='H2': #H2 product
            prodfrag+=[Chem.MolToSmarts(molfromsmiles(mappedsmiles)) for mappedsmiles in RHSdata[prodid]['mappedsmiles']]
            continue
        elif 'rxnatomidx' not in RHSdata[prodid].keys(): #SUSPECT
            unusedprod+=[prodid]
            continue
        for inst,prodatomidx in RHSdata[prodid]['rxnatomidx'].items():
            prodfrag+=[Chem.MolFragmentToSmarts(molfromsmiles(RHSdata[prodid]['mappedsmiles'][inst]),prodatomidx)]
    addstr=[]
    if farfg:
        addstr.append('Reacting functional groups are too far away in '+'species '+', '.join([str(rctid) for rctid in farfg]))
    if unusedprod:
        addstr.append('Product '+', '.join([str(ID) for ID in unusedprod])+'is not produced from reacting fragments')
    if not addstr:
        msg4='Valid'
    else:
        msg4=', '.join(addstr)
    return '>>'.join(['.'.join(reacfrag),'.'.join(prodfrag)]),LHSdata,RHSdata,msg4,farfg,unusedprod

def apply_template_row(row,inputquery):
    LHSdata=copy.deepcopy(row.LHSdata)
    templt=copy.deepcopy(row.template)
    nofg=row.nofg
    return apply_template(LHSdata,templt,inputquery=inputquery,nofg=nofg) 


def apply_template(LHSdata,templt,inputquery,nofg=[]):
    combs=[]
    imp_raw=[]
    imp_smles=[]
    imp_mol=[]
    imp_Rsmiles=[]
    msg5=''
    
    template_rxn=rdChemReactions.ReactionFromSmarts(templt,useSmiles=False)
    querycompdbin=[]
    for analoguecompd in LHSdata:
        if analoguecompd in nofg and not LHSdata[analoguecompd]['reacfrag'] : #Hydrogen will not be mapped and won't show up in reaction center 
            querycompdbin+=[{LHSdata[analoguecompd]['smiles']} for i in range(LHSdata[analoguecompd]['count'])]
            continue
        for inst,fraginf in LHSdata[analoguecompd]['reacfrag'].items():
#             breakpoint()
            querycompdset=set()
            for frag,matchidx in fraginf.items():
                for querycompd in LHSdata[analoguecompd]['querycompds'][frag]:
                    if inputquery.speciesdict[querycompd].fragsdict[frag].count<len(matchidx):
                        return [],[],[],'Species '+str(analoguecompd)+' has too many reacting functional groups' #Invalid templates as one/more reactants have more functional groups reacting relative to query
                if not querycompdset:
                    querycompdset=set(LHSdata[analoguecompd]['querycompds'][frag])
                else:
                    querycompdset=querycompdset.intersection(set(LHSdata[analoguecompd]['querycompds'][frag])) #There may be more than one query compound with match
            if querycompdset:
                querycompdbin+=[querycompdset]
        if not querycompdbin:
            return [],[],[],'Species '+str(analoguecompd)+' reacts at fragments from different query compounds'
    combs=list(itertools.product(*querycompdbin))
    try:
        imp_raw=[template_rxn.RunReactants([Chem.MolFromSmiles(query_reactant) for query_reactant in comb]) for comb in combs]
    except Exception as e:
        msg5='Template cannot be applied to reactants'
        return combs,'Error','Error',msg5
    if not imp_raw[0]:
        msg5='Template cannot be applied to reactants'
        return combs,'Error','Error',msg5
#     breakpoint()
    imp_smles=[tuple(tuple(Chem.MolToSmiles(imp) for imp in imp_prod) for imp_prod in comb) for comb in imp_raw]
    imp_smles=[{imp_prod for imp_prod in comb} for comb in imp_smles] # Remove duplicates
    combs2=[] #Check chemical validity--only add chemically valid combinations
    comb2=set()
    imp_smles2=[]
    for idx,comb in enumerate(imp_smles):
        comb2=set()
        for imp_prod in comb:
            try:
                imp_mol={tuple(molfromsmiles(imp) for imp in imp_prod)}
                comb2.add(imp_prod)
            except Exception as e:
                continue
        if comb2:
            imp_smles2+=[comb2]
            combs2+=[combs[idx]]
    if not imp_smles2:
        msg5='Impurities not chemically valid'
        return combs,imp_smles,'Error',msg5
#     try:
#         imp_mol=[{tuple(molfromsmiles(impsmles) for impsmles in impset) for impset in comb} for comb in imp_smles]
#     except Exception as e:
#         msg5='Impurities not chemically valid'
#         return combs,imp_smles,'Error',msg5
    imp_Rsmiles=[set('>>'.join(['.'.join(combs2[idx]),'.'.join(impprod)]) for impprod in imp_smles2[idx]) for idx in range(len(imp_smles2))]
    msg5='Valid'
    return combs2,imp_smles2,imp_Rsmiles,msg5

def removeduplicates(row):
    rejidx=[]
    querycompds=copy.deepcopy(row.querycompds)
    impurities=copy.deepcopy(row.impurities)
    impurityrxns=copy.deepcopy(row.impurityrxn)
    collset=copy.deepcopy(set(tuple(sorted(t)) for t in querycompds))
    for idx,comb in enumerate(querycompds):
        norm=tuple(sorted(comb))
        if norm in collset:
            collset=collset-{norm}   
        else:
            rejidx+=[idx]
    if not rejidx: #No duplicates
        return querycompds,impurities,impurityrxns
    rejidx2=[]
    for idx in rejidx:
        compareset=set(tuple(sorted(t)) for t in impurities[idx])
        otheridx=[i for i in range(len(impurities)) if i!=idx]
        otherset=set(tuple(sorted(t)) for oidx in otheridx for t in impurities[oidx])
        if compareset.issubset(otherset): #Duplicate entry is present in other query compound combinations
            rejidx2+=[idx]
    querycompds=[comb for idx,comb in enumerate(querycompds) if idx not in rejidx2]
    impurities=[comb for idx,comb in enumerate(impurities) if idx not in rejidx2]
    impurityrxns=[comb for idx,comb in enumerate(impurityrxns) if idx not in rejidx2]

    return querycompds,impurities,impurityrxns

def addconditions(row,db):
    ID=row.name
    dat=pd.read_sql_query('''SELECT ReactionID,ReagentID,Temperature,Pressure,ResidenceTime,SolventID,CatalystID from ReactionDB Where ReactionID=  "'''+ str(ID) + '''"''',db)
    dat.set_index('ReactionID',inplace=True)
    return dat


# Removing unrealistic impurities (atoms and query compounds)
def validimpurities(row,candirxnsimpfinal,conditions,hc_prod,querycompds):
    impset=copy.deepcopy(set(row.impurities))
    ID=row.name
    if candirxnsimpfinal.loc[ID].hcprod:
        impset=impset-set([hc_prod[hcid]['smiles'] for hcid in candirxnsimpfinal.loc[ID].hcprod])
    if all([imp in row.querycompds for imp in impset]):
        return 'No transformation of interest'
    elif (conditions.loc[ID].ReagentID.values[0]=='NaN' or conditions.loc[ID].ReagentID.values[0] is None) and (conditions.loc[ID].SolventID.values[0]=='NaN' or conditions.loc[ID].SolventID.values[0] is None) and (conditions.loc[ID].CatalystID.values[0]=='NaN' or conditions.loc[ID].CatalystID.values[0] is None) and conditions.loc[ID].Temperature.values[0] is None and conditions.loc[ID].Pressure.values[0] is None and conditions.loc[ID].ResidenceTime.values[0] is None:
        return 'No reaction conditions detected. Check reaction record to verify plausibility'
    elif all([imp in querycompds for imp in impset]):
        return 'Query compounds suggested as impurities'
    elif any([Descriptors.NumRadicalElectrons(Chem.MolFromSmiles(imp))>0 for imp in impset]):
#     elif any(['[' in imp[0] and ']' in imp[-1] for imp in impset]):
        return 'Impurities are atoms/radicals'
    elif (len(set(row.querycompds))==1 and (conditions.loc[ID].ReagentID.values[0]=='NaN' or conditions.loc[ID].ReagentID.values[0] is None)) or (len(set(candirxnsimpfinal.loc[ID].LHS))==1 and (conditions.loc[ID].ReagentID.values[0]=='NaN' or conditions.loc[ID].ReagentID.values[0] is None)):
        return 'Self-reaction/single reactant detected with no reagents. Check reaction record to verify plausibility.'
    else:
        return 'Valid'

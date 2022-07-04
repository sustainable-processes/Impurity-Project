# %load ./ApplyTempl.py

def applytemplate(analoguerxnstemplfilt,inputquery,ncpus=16,restart=True):
#     breakpoint()
    if ncpus>1:
        if restart:
            initray(num_cpus=ncpus)
        analoguerxnstemplfiltdis=mpd.DataFrame(analoguerxnstemplfilt)
    else:
        analoguerxnstemplfiltdis=analoguerxnstemplfilt
    impurities=analoguerxnstemplfiltdis.apply(apply_template_row,inputquery=inputquery,axis=1,result_type='reduce')
    impurities=pd.Series(data=impurities.values,index=impurities.index) #Optional convert modin back to pandas
    impuritiesdf=pd.DataFrame(data=impurities.tolist(),index=impurities.index,columns=['querycompds','impurities','impurityrxn','msg5'])
    analoguerxnsimp=copy.deepcopy(analoguerxnstemplfilt)
    analoguerxnsimp[['querycompds','impurities','impurityrxn','msg5']]=impuritiesdf
    return analoguerxnsimp


def apply_template_row(row,inputquery):
    LHSdata=copy.deepcopy(row.LHSdata)
    templt=copy.deepcopy(row.template)
    unusedanalogue=row.unusedanalogue
    return apply_template(LHSdata,templt,inputquery,unusedanalogue=unusedanalogue) 


def apply_template(LHSdata,templt,inputquery,unusedanalogue=[]):
#     breakpoint()
    combs=[]
    imp_raw=[]
    imp_smles=[]
    imp_mol=[]
    imp_Rsmiles=[]
    diffreac=set()
    msg5=''
#     breakpoint()
    template_rxn=rdChemReactions.ReactionFromSmarts(templt,useSmiles=False)
    querycompdbin=[]
    for analoguecompd in LHSdata:
        LHSdata_=LHSdata[analoguecompd]
        if 'reacfrag' not in LHSdata_ or not LHSdata_['reacfrag']: #Hydrogen will not be mapped and won't show up in reaction center 
            if analoguecompd not in unusedanalogue:
                querycompdbin+=[{LHSdata_['smiles']} for i in range(LHSdata_['count'])]
            continue
        for inst,fraginf in LHSdata_['reacfrag'].items():
#             breakpoint()
            querycompdset=set()
            for frag,matchidx in fraginf.items():
                for querycompd in LHSdata_['querycompds'][frag]:
                    if inputquery['species'][querycompd][frag]['count']<len(matchidx):
                        return [],[],[],'Species '+str(analoguecompd)+' has too many reacting functional groups' #Invalid templates as one/more reactants have more functional groups reacting relative to query
                if not querycompdset:
                    querycompdset=set(LHSdata_['querycompds'][frag])
                else:
                    querycompdset2=querycompdset.intersection(set(LHSdata[analoguecompd]['querycompds'][frag])) #There may be more than one query compound with match
                    if not querycompdset2:
                        diffreac.add(analoguecompd)
                        querycompdset=querycompdset.union(set(LHSdata[analoguecompd]['querycompds'][frag]))
                    else:
                        querycompdset=querycompdset2
            querycompdbin+=[querycompdset]            
#             if querycompdset:
#                 querycompdbin+=[querycompdset]
#         if not querycompdbin:
#             return [],[],[],'Species '+str(analoguecompd)+' reacts at fragments from different query compounds'
#     breakpoint()
    if diffreac:
        msg5='Species reacts at fragments from different query compounds: '+','.join([str(rct) for rct in diffreac])
    combs=list(itertools.product(*querycompdbin)) #Buggy if too many (keep limit at 6)
    try:
        #Addition of Chem.AddHs significantly impacts results
        imp_raw=[]
        for comb in combs:
            querymols=[]
            for query_reactant in comb:
                querymol=molfromsmiles(query_reactant)
                clear_isotope(querymol)
                querymol=Chem.RemoveHs(querymol)
                querymols.append(querymol)
            imp_raw.append(template_rxn.RunReactants(querymols))     
#         imp_raw=[template_rxn.RunReactants([Chem.AddHs(molfromsmiles(query_reactant)) for query_reactant in comb]) for comb in combs]
    except Exception as e:
        if msg5:
            msg5=msg5+', '+'Template cannot be applied to reactants'
        else:
            msg5='Template cannot be applied to reactants'
        return combs,'Error','Error',msg5
    if not imp_raw[0]:
        if msg5:
            msg5=msg5+', '+'Template cannot be applied to reactants'
        else:
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
        if msg5:
            msg5=msg5+', '+'Impurities not chemically valid'
        else:
            msg5='Impurities not chemically valid'
        return combs,imp_smles,'Error',msg5
    imp_Rsmiles=[set('>>'.join(['.'.join(combs2[idx]),'.'.join(impprod)]) for impprod in imp_smles2[idx]) for idx in range(len(imp_smles2))]
    if not msg5:
        msg5='Valid'
    return combs2,imp_smles2,imp_Rsmiles,msg5

def removeduplicates(analoguerxnsimpfinal,ncpus=16,restart=False):
    if ncpus>1:
        if restart:
            initray(num_cpus=ncpus)
        analoguerxnsimpfinaldis=mpd.DataFrame(analoguerxnsimpfinal)
    else:
        analoguerxnsimpfinaldis=analoguerxnsimpfinal
    duplicatesrem=analoguerxnsimpfinaldis.apply(removeduplicates_,axis=1,result_type='reduce') #remove duplicates
    duplicatesrem=pd.Series(data=duplicatesrem.values,index=duplicatesrem.index)
    duplicatesremdf=pd.DataFrame(data=duplicatesrem.tolist(),index=duplicatesrem.index,columns=['querycompds','impurities','impurityrxn'])
    analoguerxnsimpfinal[['querycompds','impurities','impurityrxn']]=duplicatesremdf
    return analoguerxnsimpfinal

def removeduplicates_(row):
#     breakpoint()
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


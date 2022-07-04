# %load ./ImpurityCleaning.py

def cleanimpurities(analoguerxnsimpfinal,analoguerxnsfilt,analoguerxns_updated,inputquery,includesolv=True,reaxys_update=True,hc_prod={}):
    impfinal=analoguerxnsimpfinal[['querycompds','impurities','impurityrxn','rxnsmiles0','balrxnsmiles','LHS','RHS','hcrct',
                                   'hcprod','LHSdata','RHSdata','NumRefs','NumSteps','NumStages','rxncentermapnum','specmap','rnbmap',
                                   'mapped_rxn','msg1','template','msg4']].copy()
    if includesolv:
        impfinal=updatecolumns(analoguerxns_updated,impfinal,cols=['Rgtdata','Solvdata'])
    else:
        impfinal=updatecolumns(analoguerxns_updated,impfinal,cols=['Rgtdata'])
    if reaxys_update:
        if includesolv:
            impfinal=updatecolumns(analoguerxnsfilt,impfinal,cols=['ReagentID','Temperature','Pressure','ReactionTime',
              'SolventID','MissingSolvent','CatalystID','MissingCatalyst','NameDict','ConditionNotes',
              'ReactionType','YearPublished','Yield'],config=['querycompds','impurities','impurityrxn',
                                                            'template','msg4','rxncentermapnum','specmap','rnbmap',
                                                            'mapped_rxn','msg1','NumRefs','NumSteps','NumStages',
                                                            'rxnsmiles0','balrxnsmiles','LHS','hcrct','LHSdata','RHS','hcprod','RHSdata',
                                                            'ReagentID','Rgtdata','CatalystID','MissingCatalyst',
                                                            'SolventID','MissingSolvent','NameDict','Temperature',
                                                            'Pressure','ReactionTime','ConditionNotes','ReactionType',
                                                            'Yield','YearPublished'])
        else:
            impfinal=updatecolumns(analoguerxnsfilt,impfinal,cols=['ReagentID','Temperature','Pressure','ReactionTime',
             'CatalystID','MissingCatalyst','NameDict','ConditionNotes','ReactionType','YearPublished','Yield'],
                                   config=['querycompds','impurities','impurityrxn',
                                                            'template','msg4','rxncentermapnum','specmap','rnbmap',
                                                            'mapped_rxn','msg1','NumRefs','NumSteps','NumStages',
                                                            'rxnsmiles0','balrxnsmiles','LHS','hcrct','LHSdata','RHS','hcprod','RHSdata',
                                                            'ReagentID','Rgtdata','CatalystID','MissingCatalyst',
                                                            'NameDict','Temperature','Pressure','ReactionTime',
                                                            'ConditionNotes','ReactionType','Yield','YearPublished'])
            
    else:
        if includesolv:
            impfinal=updatecolumns(analoguerxnsfilt,impfinal,cols=['ReagentID','SolventID','Temperature','Pressure','ReactionTime',
                                                                'CatalystID','ReactionType'],
                                                                config=['querycompds','impurities',
                                                                        'impurityrxn', 'template','msg4','rxncentermapnum',
                                                                        'specmap','rnbmap','mapped_rxn','msg1','NumRefs',
                                                                        'NumSteps','NumStages','rxnsmiles0','balrxnsmiles',
                                                                        'LHS','hcrct', 'LHSdata','RHS','hcprod','RHSdata','ReagentID',
                                                                        'Rgtdata','CatalystID','SolventID','Temperature',
                                                                        'Pressure','ReactionTime','ReactionType'])
        else:
            impfinal=updatecolumns(analoguerxnsfilt,impfinal,cols=['ReagentID','Temperature','Pressure','ReactionTime',
                                                                'CatalystID','ReactionType'],
                                                                config=['querycompds','impurities',
                                                                        'impurityrxn', 'template','msg4','rxncentermapnum',
                                                                        'specmap','rnbmap','mapped_rxn','msg1','NumRefs',
                                                                        'NumSteps','NumStages','rxnsmiles0','balrxnsmiles',
                                                                        'LHS','hcrct','LHSdata','RHS','hcprod','RHSdata','ReagentID',
                                                                        'Rgtdata','CatalystID','Temperature',
                                                                        'Pressure','ReactionTime','ReactionType'])
            

    impfinal=impfinal.explode('querycompds')
    # query
    imp1=analoguerxnsimpfinal[['impurities']].copy()
    impquery=imp1.explode('impurities') #Explode impurity set list per query compound set
    # impquery
    imprxn1=analoguerxnsimpfinal[['impurityrxn']].copy()
    imprxnquery=imprxn1.explode('impurityrxn')
    # imprxnquery
    impfinal['impurities']=impquery
    impfinal=impfinal.explode('impurities')# Explode again to get each impurity set
    impfinal['impurityrxn']=imprxnquery.explode('impurityrxn') # Explode again to get each impurity reaction
    impfinal=checkimpurities(impfinal,inputquery,hc_prod,reaxys_update=reaxys_update)
    return impfinal


def addconditions(row,db):
    ID=row.name
    dat=pd.read_sql_query('''SELECT ReactionID,ReagentID,Temperature,Pressure,ResidenceTime,SolventID,CatalystID from ReactionDB Where ReactionID=  "'''+ str(ID) + '''"''',db)
    dat.set_index('ReactionID',inplace=True)
    return dat

def checkimpurities(impfinal,inputquery,hc_prod={},reaxys_update=True,ncpus=16,restart=False):
    querycompds=list(inputquery['species'].keys())
    if ncpus>1:
        if restart:
            initray(num_cpus=ncpus)
        impfinaldis=mpd.DataFrame(impfinal)
    else:
        analoguerxnsimpfinaldis=analoguerxnsimpfinal
    validimp=impfinaldis.apply(valid_impurities,hc_prod=hc_prod,querycompds=querycompds,reaxys_update=reaxys_update,axis=1,result_type='reduce')
    validimp=pd.DataFrame(data=validimp,index=validimp.index,columns=['msg6'])
    impfinal['msg6']=validimp
    return impfinal

# Removing unrealistic impurities (atoms and query compounds)
def valid_impurities(row,hc_prod,querycompds,reaxys_update=True):
    impset=copy.deepcopy(set(row.impurities))
    if row.hcprod:
        impset=impset-set([hc_prod[hcid]['smiles'] for hcid in row.hcprod])
    if all([imp in row.querycompds for imp in impset]):
        return 'No transformation of interest'
    if any([Descriptors.NumRadicalElectrons(Chem.MolFromSmiles(imp))>0 for imp in impset]):
        return 'Impurities are atoms/radicals'
    if reaxys_update:
        if (row.ReagentID=='NaN' or row.ReagentID is None or not row.ReagentID) and ((row.SolventID=='NaN' or row.SolventID is None or not row.SolventID) and not row.MissingSolvent) and ((row.CatalystID=='NaN' or row.CatalystID is None or not row.CatalystID) and not row.MissingCatalyst) and (row.Temperature is None or not row.Temperature) and (row.Pressure is None or not row.Pressure) and (row.ReactionTime is None or not row.ReactionTime):
            return 'No reaction conditions detected. Check reaction record to verify plausibility'
    else:
        if (row.ReagentID=='NaN' or row.ReagentID is None or not row.ReagentID) and (row.SolventID=='NaN' or row.SolventID is None or not row.SolventID) and (row.CatalystID=='NaN' or row.CatalystID is None or not row.CatalystID) and (row.Temperature is None or not row.Temperature) and (row.Pressure is None or not row.Pressure) and (row.ReactionTime is None or not row.ReactionTime):
            return 'No reaction conditions detected. Check reaction record to verify plausibility'
    if reaxys_update:
        if (len(set(row.querycompds))==1 or len(set(row.LHS))==1) and ((row.ReagentID=='NaN' or row.ReagentID is None or not row.ReagentID) and (row.CatalystID=='NaN' or row.CatalystID is None or not row.CatalystID) and not row.MissingCatalyst):
            return 'Self-reaction/single reactant detected with no reagents and catalysts. Check reaction record to verify plausibility.'
    else:
        if (len(set(row.querycompds))==1 or len(set(row.LHS))==1) and ((row.ReagentID=='NaN' or row.ReagentID is None or not row.ReagentID) and (row.CatalystID=='NaN' or row.CatalystID is None or not row.CatalystID)):
            return 'Self-reaction/single reactant detected with no reagents and catalysts. Check reaction record to verify plausibility.'
    if all([imp in querycompds for imp in impset]):
        return 'Query compounds suggested as impurities'    
    return 'Valid'




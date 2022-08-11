# %load ./MainScript.py

import ray
import modin.pandas as mpd
from AnalgCompds import processquery
from MainFunctions import CustomError, drawReaction, openpickle, writepickle
from ReaxysAPIv2 import initray
from helpCompound import hc_Dict,hc_rct
import matplotlib.pyplot as plt
import pandas as pd
import os



# Main Script
masterdbreadpath='/home/aa2133/Impurity-Project/Reaxys_Data/' # Master data source path
datasources={'fraggroupssource':masterdbreadpath+'fragfreq.pickle',
                'fragdbsource':masterdbreadpath+'fragdb.pickle',
                'rxnsource':masterdbreadpath+'ReactionDB.pickle',
                'exemptiondir':masterdbreadpath+'PotCatList.pickle',
                'unresolveddir':masterdbreadpath+'UnresolvedIDs.pickle',
                'substancesource': masterdbreadpath+'SubstanceSmiles.pickle',
                'SQLdatsource': masterdbreadpath+'SQL/Reaxys_Data.db'}
masterparams={'folderwritepath':'/home/aa2133/Impurity-Project/Input/',
              'ncpus':16,
             'showresults':False,
             'writetofile':True,
             'reaxys_update':True}

step2params={'debug':False,
            'iqfilename':'inputquery'}

step3params={'inputquery':None,
            'inputquery_analg':None,
            'SQL':False,
            'refquery':None,
            'includefragparents':False,
            'onlyisotopes':True,
             'iqanalgfilename':'inputquery_analg',
             'fragdict':None,
             'fragdictfilename':'fragdict',
            'fragchoice':{},
            'SimThresh':None,
            'MWThresh':None,
             'inputquery_analg_updated':None,
             'iqanalgupdatedfilename':'inputquery_analg_updated',
             'nomixtures':True,
             'combinedpoolfilename':'combinedpool',
             'combinedpooldffilename':'combinedpooldf'
            }

step4params={'combinedpool':None,
             'combinedpoolex':None,
             'exemptionlist':None,
            'combinedpoolexfilename':'combinedpoolex',
             'workflow':'strict',
             'returnall':True,
             'refanaloguerxns':None,
            'analgidsdictfilename':'analgidsdict',
             'analgrxnsdictfilename':'analgrxnsdict',
             'analoguerxnsrawfilename':'analoguerxnsraw',
             'analoguerxnsfilename':'analoguerxns',
             'analoguerxns':None
            }

step5params={'unresolvedids':None,
             'analoguerxnsfilt':None,
             'analoguerxnsfiltfilename':'analoguerxnsfilt',
             'refanaloguerxns_updated':None,
             'includesolv':True,
             'hc_Dict':hc_Dict,
             'hc_rct':hc_rct,
             'analoguerxns_updatedfilename':'analoguerxns_updated',
             'analoguerxns_updated':None
            }
step6params={'refbalrxns':None,
            'helpprod':True,
            'helpreact':False,
            'addrctonly':False,
            'ignoreH':False,
            'analoguerxnsbal':None,
             'analoguerxnsbalfilename':'analoguerxnsbal',
             'balrxnsfilename':'balrxns',
             'removeLHSinsufficient':True,
             'removeRHSinsufficient':True,
             'removehcarriers':True,
             'removechargeimbalance':True,
             'analoguerxnsbalfilt':None,
             'analoguerxnsbalfiltfilename':'analoguerxnsbalfilt'
            }
step7params={'analoguerxnsmapped':None,
             'refmappedrxns':None,
            'analoguerxnsmappedfilename':'analoguerxnsmapped',
            'analoguerxnsparsed':None,
             'refparsedrxns':None,
            'analoguerxnsparsedfilename':'analoguerxnsparsed',
             'analoguerxnsparsedfilt':None,
            'analoguerxnsparsedfiltfilename':'analoguerxnsparsedfilt',
            'analoguerxnsassigned':None,
            'analoguerxnsassignedfilename':'analoguerxnsassigned',
            'analoguerxnsassignedfilt':None,
            'analoguerxnsassignedfiltfilename':'analoguerxnsassignedfilt',
             'removemappingerror':True,
            'removesmilesdiscrepancy':True,
            'removemandatoryunmapped':True}
step8params={'analoguerxnscent':None,
            'analoguerxnscentfilename':'analoguerxnscent',
            'analoguerxnscentfilt':None,
            'analoguerxnscentfiltfilename':'analoguerxnscentfilt',
             'analoguerxnsvalid':None,
             'analoguerxnsvalidfilename':'analoguerxnsvalid',
             'analoguerxnsfinal':None,
             'analoguerxnsfinalfilename':'analoguerxnsfinal',
            'removeoutfrags':True}
step9params={'analoguerxnstempl':None,
            'analoguerxnstemplfilename':'analoguerxnstempl',
            'analoguerxnstemplfilt':None,
            'analoguerxnstemplfiltfilename':'analoguerxnstemplfilt',
            'onlyvalidtempl':True,
             'removefarfg':True,
             'removeunusedprod':True,
             'specificity':'loose'
            }
step10params={'analoguerxnsimp':None,
              'analoguerxnsimpfilename':'analoguerxnsimp',
             'analoguerxnsimpfilt':None,
             'analoguerxnsimpfiltfilename':'analoguerxnsimpfilt'}
step11params={'impfinal':None,
             'impfinalfilename':'impfinal',
             'impfinalfilt':None,
             'impfinalfiltfilename':'impfinalfilt'}
step12params={}

step13params={'impfinalfilt5':None,
               'includefraginfo':True,
              'impfinalfilt2':None,
              'impfinalfilt3':None,
              'impfinalfilt4':None,
             'impfinalfilt5filename':'impfinalfilt5',
             'impfinalfilt2filename':'impfinalfilt2',
             'impfinalfilt3filename':'impfinalfilt3',
             'impfinalfilt4filename':'impfinalfilt4'
             }

step14params={'summary3':None,
             'summaryfilename':'summary',
             'summary2filename':'summary2',
             'summary3filename':'summary3',
              'relquantile':0.9,     
              'uppertquantile':0.95, 
              'lowertquantile':0.05
             }

step15params={'analoguedisplay':2,
             'missedanaloguedisplay':2
             }

inputparams={**datasources,**masterparams,**step2params,**step3params,**step4params,**step5params,**step6params,
            **step7params,**step8params,**step9params,**step10params,**step11params,**step12params,**step13params,
            **step14params,**step15params}

    
def main(casename,steps,userinput='',catalyst=[],Trange=None,conditions=[],IP=inputparams,**kwargs):
    '''
    casename: case folder
    steps: All steps to execute
    userinput: Input reaction SMILES
    catalyst: Specify catalyst SMILES or text (recommended to be as general as possible) *optional*
    Trange: Temperature range ([lowerbound, upperbound]) of reaction *optional*
    IP: Dictionary of default input parameters and stored results
    **kwargs: Optional user parameters that will replace default parameters in IP
    '''
#     breakpoint()
    if kwargs:
        for key,val in kwargs.items():
            if key in IP:
                IP[key]=val
#     breakpoint()
    # Creating directories to store output (Edit if any changes)
    stages=['DataMining','DataProcessing','ImpurityPrediction','ImpurityRanking']
    casedir=IP['folderwritepath']+casename+'/'
    dmdir=casedir+stages[0]+'/'
    IP['dmdir']=dmdir
    dpdir=casedir+stages[1]+'/'
    IP['dpdir']=dpdir
    ipdir=casedir+stages[2]+'/'
    IP['ipdir']=ipdir
    irdir=casedir+stages[3]+'/'
    IP['irdir']=irdir
    for fdir in [dmdir,dpdir,ipdir,irdir]:
        if not os.path.exists(fdir):
            os.makedirs(fdir)
    import gc
    gc.collect()
    
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Main workflow %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Step 1 & 2  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Step 1 involves inputting query species, here defined as user-supplied reactants, reagents, solvents, and main products
# in the form of SMILES. These are already passed in under userinput

# Step 2 involves processing the input SMILES (userinput), and identifying functional group (FG) fragments in each 
# user-supplied query species. This information is stored in an easily queried dictionary. 
    if 2 in steps: 
        IP=step2(userinput,IP)
        import gc
        gc.collect()
            
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Step 3  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          
# Step 3 involves finding analogue species in Reaxys that contain any of the identified carrier fragments

    if 3 in steps:
        IP=step3(IP)
        import gc
        gc.collect()
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Step 4  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          
# In step 4, analogue reactions that only contain analogue or query species are extracted from the Reaxys data.
    if 4 in steps:
        IP=step4(IP)
        import gc
        gc.collect()        
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Step 5  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Step 5 is the pre-processing, cleaning and filtering the obtained analogue reactions.
    if 5 in steps:
        IP=step5(IP)
        import gc #Garbage collection to free up memory
        gc.collect()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Step 6  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Step 6 resolves the issue of incomplete reaction structures that can be commonly found in Reaxys. Balances reactions
    if 6 in steps:
        IP=step6(IP)
        import gc
        gc.collect()        
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Step 7  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Step 7 involves atom-to-atom mapping to establish a one-to-one correspondence between reactant and product atoms
    if 7 in steps:
        IP=step7(IP)
        import gc
        gc.collect()        
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Step 8  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# In step 8, reaction centres of the remaining analogue reactions are obtained to ensure that they are within designated carrier fragments           
    if 8 in steps:
        IP=step8(IP)
        import gc
        gc.collect()
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Step 9  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Step 9 refers to template generation.
    if 9 in steps:
        IP=step9(IP)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Step 10  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Step 10 refers to template application.        
    if 10 in steps:
        IP=step10(IP)
        import gc
        gc.collect()
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Step 11  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Step 11 refers to impurity cleaning and filtering.  
    if 11 in steps:
        IP=step11(IP)
        import gc
        gc.collect()
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Step 12  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Step 12 refers to morgan fingerprint similarity calculations.      
    if 12 in steps:
        IP=step12(IP)
        import gc
        gc.collect()
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Step 13  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Step 13 refers to reaction condition filtering.    
    if 13 in steps:
        IP=step13(IP,conditions=conditions,catalyst=catalyst)
        import gc
        gc.collect()
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Step 14  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Step 14 refers to impurity ranking.      
    if 14 in steps:
        IP=step14(IP,Trange=Trange)
        import gc
        gc.collect()
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Step 15  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Step 15 refers to result visualization.  
    if 15 in steps:
#             breakpoint()
        impfinalfilt=IP['impfinalfilt']
        if impfinalfilt is None:
            impfinalfilt=pd.read_pickle(IP['irdir']+IP['impfinalfiltfilename']+'.pickle')
        inputquery=IP['inputquery']
        if inputquery is None:
            inputquery=openpickle(IP['dmdir']+IP['iqfilename']+'.pickle')
        summary3=IP['summary3']
        if summary3 is None:
            summary3=pd.read_pickle(IP['irdir']+IP['summary3filename']+'.pickle')
        impfinalfilt2=IP['impfinalfilt2']
        if impfinalfilt2 is None:
            impfinalfilt2=pd.read_pickle(IP['irdir']+IP['impfinalfilt2filename']+'.pickle')
        impfinalfilt3=IP['impfinalfilt3']
        if impfinalfilt3 is None:
            impfinalfilt3=pd.read_pickle(IP['irdir']+IP['impfinalfilt3filename']+'.pickle')
        impfinalfilt4=IP['impfinalfilt4']
        if impfinalfilt4 is None:
            impfinalfilt4=pd.read_pickle(IP['irdir']+IP['impfinalfilt4filename']+'.pickle')
        impfinalfilt5=IP['impfinalfilt5']
        if impfinalfilt5 is None:
            impfinalfilt5=pd.read_pickle(IP['irdir']+IP['impfinalfilt5filename']+'.pickle')
        impfinalfilt['Relevance_morgan']=impfinalfilt['Relevance_morgan'].apply(lambda x: round(x,2))
        print('Query reaction (user input):')
        display(drawReaction(rdChemReactions.ReactionFromSmarts(inputquery['smiles'],useSmiles=True)))
        hc_Dict=IP['hc_Dict']
        for i,rxn in enumerate(summary3.rxn):
            if i==0:
                if not Trange:
                    Trange=summary3.loc[summary3.rxn==rxn].t_range.values[0]
                print('Temperature range (C) of '+str(min(Trange))+' to '+str(max(Trange)))
                mainprods=summary3.loc[summary3.rxn==rxn].products.values[0]
            Trangei=summary3.loc[summary3.rxn==rxn].t_range.values[0]
            if min(Trangei)>max(Trange) or max(Trangei)<min(Trange):
                continue
            impprods=summary3.loc[summary3.rxn==rxn].products.values[0]
            if i>0 and any([impprod in mainprods and impprod not in [hc_Dict[j]['smiles'] for j in hc_Dict] for impprod in impprods]):
                continue
            display(drawReaction(rdChemReactions.ReactionFromSmarts(rxn,useSmiles=True)))
            max_relevance=summary3.loc[summary3.rxn==rxn]['Max relevance_tfiltered'].values[0]
            print('Max relevance: '+str(max_relevance))
            print('Smarts: '+rxn)
            print('Number of hits: '+str(summary3.loc[summary3.rxn==rxn]['Hits_tfiltered'].values[0]))
            print('Temperature range (10th to 90th %tile): '+'-'.join([str(round(t,1)) for t in Trangei]))
            frame4=summary3.loc[summary3.rxn==rxn]['Frame4'].values[0]
            frame2=summary3.loc[summary3.rxn==rxn]['Frame3'].values[0]
            frame=summary3.loc[summary3.rxn==rxn]['Frame'].values[0]
            impframe_old=impfinalfilt.loc[impfinalfilt.impurityrxn==rxn].sort_values(by='Relevance_morgan',ascending=False)
            for i in range(IP['analoguedisplay']):
                print('Analogue reactions with max relevance:')
                print('Reaction ID: '+str(frame4.iloc[i].name[0]))
                print('Instance: '+str(frame4.iloc[i].name[1]))
                print('Relevance: '+str(frame4.iloc[i].Relevance_morgan))
                display(drawReaction(rdChemReactions.ReactionFromSmarts(frame4.iloc[i].mapped_rxn,useSmiles=True)))
                idx=frame4.iloc[i].name
                rgts=set(frame4.iloc[i].ReagentID)-set(frame4.iloc[i].LHS)
                cats=frame4.iloc[i].CatalystID2
                solvs=frame4.iloc[i].SolventID
                print('Reagents: '+', '.join([frame4.iloc[i].NameDict[ID] for ID in rgts if ID not in cats]))
                print('Catalyst: '+', '.join([frame4.iloc[i].NameDict[ID] for ID in cats])+', '.join(frame4.iloc[i].MissingCatalyst))
                print('Solvents: '+', '.join([frame4.iloc[i].NameDict[ID] for ID in solvs])+', '.join(frame4.iloc[i].MissingSolvent))
                tempr=impframe_old.loc[idx].Temperature
                if tempr:
                    if type(tempr)==list:   
                        tempr=tempr[0]
                    print('Temperature (C): '+str(tempr))
                print('Condition notes: '+frame4.iloc[i].ConditionNotes)


            errorcodes=['Invalid condition','Suspect self-reaction (without catalyst or reagent)','Invalid catalyst',
                       'Missing reactants/products','Missing temperature','Temperature outside indicated range']
            for i in range(IP['missedanaloguedisplay']):
                try:
                    maxrejrel=impframe_old.iloc[i].Relevance_morgan
                    maxrejid=impframe_old.iloc[i].name
                except Exception:
                    break
                if maxrejrel>max_relevance:
                    for k,refdb in enumerate([impfinalfilt2,impfinalfilt3,impfinalfilt4,impfinalfilt5,frame2,frame4]):
                        if maxrejid in refdb.index:
                            continue
                        else:
                            errorcode=errorcodes[k]
                            break

                    print('Removed analogue reactions with max relevance:')
                    print('Reaction ID: '+str(impframe_old.iloc[i].name[0]))
                    print('Instance: '+str(frame4.iloc[i].name[1]))
                    print('Relevance: '+str(impframe_old.iloc[i].Relevance_morgan))
                    print('Reaction removed due to: '+errorcode)
                    display(drawReaction(rdChemReactions.ReactionFromSmarts(impframe_old.iloc[i].mapped_rxn,useSmiles=True)))
                    rgts=set(impframe_old.iloc[i].ReagentID)-set(impframe_old.iloc[i].LHS)
                    cats=impframe_old.iloc[i].CatalystID2
                    solvs=impframe_old.iloc[i].SolventID
                    print('Reagents: '+', '.join([impframe_old.iloc[i].NameDict[ID] for ID in rgts if ID not in cats]))
                    print('Catalyst: '+', '.join([impframe_old.iloc[i].NameDict[ID] for ID in cats])+', '.join(impframe_old.iloc[i].MissingCatalyst))
                    print('Solvents: '+', '.join([impframe_old.iloc[i].NameDict[ID] for ID in solvs])+', '.join(impframe_old.iloc[i].MissingSolvent))
                    tempr=impframe_old.iloc[i].Temperature
                    if tempr:
                        if type(tempr)==list:   
                            tempr=tempr[0]
                        print('Temperature (C): '+str(tempr))
                    print('Condition notes: '+impframe_old.iloc[i].ConditionNotes) 
    if IP['showresults']:
        return IP

           
def step2(userinput,IP):
    if userinput:#User has supplied a reaction smiles string
        if not type(userinput)==str:
            raise CustomError('Please provide a query reaction smiles')
        inputsmiles=userinput
        inputquery=processquery(inputsmiles,debug=IP['debug'])
        if IP['writetofile']:
            writepickle(inputquery,IP['dmdir'],IP['iqfilename'])
        IP['inputquery']=inputquery
        print('Step 2 complete')
    return IP

def step3(IP):
    inputquery=IP['inputquery']
    inputquery_analg=IP['inputquery_analg']
    inputquery_analg_updated=IP['inputquery_analg_updated']
    if inputquery_analg is None:
        if inputquery is None:
            inputquery=openpickle(IP['dmdir']+IP['iqfilename']+'.pickle')
            IP['inputquery']=inputquery
        inputquery_analg,fragdict=getanaloguespecies(inputquery,IP['fragdbsource'],SQL=IP['SQL'],
                                                    refquery=IP['refquery'],ncpus=IP['ncpus'],
                                                    fragtable=IP['fraggroupssource'],substancedbsource=IP['substancesource'],
                                                    includefragparents=IP['includefragparents'],onlyisotopes=IP['onlyisotopes'])
        if IP['writetofile']:
            writepickle(inputquery_analg,IP['dmdir'],IP['iqanalgfilename'])
            writepickle(fragdict,IP['dmdir'],IP['fragdictfilename'])
        if IP['showresults']:
            IP['fragdict']=fragdict
    if IP['SimThresh'] is not None or IP['MWThresh'] is not None:
        if inputquery_analg_updated is None:
            similarity=False
            molwt=False
            if IP['SimThresh'] is not None:
                similarity=True
            if IP['MWThresh'] is not None:
                molwt=True
            inputquery_analg_updated=updatequery(inputquery_analg,fragchoice=IP['fragchoice'],similarity=similarity,
                                                 fingerprint='morgan',morganradius=2,addHs=True,molwt=molwt,ncpus=IP['ncpus'])
            if writetofile:
                writepickle(inputquery_analg_updated,IP['dmdir'],IP['iqanalgupdatedfilename'])
            if IP['showresults']:
                IP['inputquery_analg_updated']=inputquery_analg_updated
    else:
        if IP['showresults']:
            IP['inputquery_analg']=inputquery_analg
        inputquery_analg_updated=inputquery_analg
        
    combinedpooldf=getcombinedpool(inputquery_analg_updated,fragchoice=IP['fragchoice'],ST=IP['SimThresh'],
                                   fingerprint='morgan',morganradius=2,MWT=IP['MWThresh'],nomixtures=IP['nomixtures'],
                                   res_format='df')
    combinedpool=getcombinedpool(inputquery_analg_updated,fragchoice=IP['fragchoice'],ST=IP['SimThresh'],
                                 fingerprint='morgan', morganradius=2,MWT=IP['MWThresh'],nomixtures=IP['nomixtures'],
                                 res_format='list')
    if writetofile:
        pd.to_pickle(combinedpooldf,IP['dmdir']+IP['combinedpooldffilename']+'.pickle')
        writepickle(combinedpool,IP['dmdir'],IP['combinedpoolfilename'])
    IP['combinedpool']=combinedpool
    print('Step 3 complete')
    return IP

def step4(IP): #MEMORY INTENSIVE
    combinedpool=IP['combinedpool']
    combinedpoolex=IP['combinedpoolex']
    if combinedpoolex is None:
        if combinedpool is None:
            combinedpool=openpickle(IP['dmdir']+IP['combinedpoolfilename']+'.pickle')
        exemptionlist=openpickle(IP['exemptiondir']) #Catalysts misclassified as reagents (uncertain)
        IP['exemptionlist']=exemptionlist
        combinedpoolex=updatecombinedpool(combinedpool,exemptionlist=exemptionlist)
        if IP['writetofile']:
            writepickle(combinedpoolex,IP['dmdir'],IP['combinedpoolexfilename'])
        IP['combinedpoolex']=combinedpoolex
    analgidsdict,analgrxnsdict,analoguerxnsraw=getanaloguerxns(IP['rxnsource'],combinedpool,combinedpoolex=combinedpoolex, workflow=IP['workflow'],
                              returnall=IP['returnall'],reaxys_update=IP['reaxys_update'],refanaloguerxns=IP['refanaloguerxns'],
                            ncpus=IP['ncpus'])
    if IP['writetofile']:
        writepickle(analgidsdict,IP['dmdir'],IP['analgidsdictfilename'])
        writepickle(analgrxnsdict,IP['dmdir'],IP['analgrxnsdictfilename'])
        if not analoguerxnsraw.empty:
            pd.to_pickle(analoguerxnsraw,IP['dmdir']+IP['analoguerxnsrawfilename']+'.pickle')
        pd.to_pickle(analgrxnsdict[IP['workflow']],IP['dmdir']+IP['analoguerxnsfilename']+'.pickle')
    if IP['showresults']:
        IP['analoguerxns']=analgrxnsdict[IP['workflow']]
    print('Step 4 complete')
    return IP

def step5(IP): #MEMORY INTENSIVE
    analoguerxns=IP['analoguerxns']
    exemptionlist=IP['exemptionlist']
    analoguerxnsfilt=IP['analoguerxnsfilt']
    if analoguerxnsfilt is None:
        if analoguerxns is None:
            analoguerxns=pd.read_pickle(IP['dmdir']+IP['analoguerxnsfilename']+'.pickle')
        if exemptionlist is None:
            exemptionlist=openpickle(IP['exemptiondir'])
        unresolvedids=openpickle(IP['unresolveddir'])
        analoguerxnsfilt=filteranaloguerxns(analoguerxns,unresolvedids,reaxys_update=IP['reaxys_update'],exemptionlist=exemptionlist)
        if IP['writetofile']:
            pd.to_pickle(analoguerxnsfilt,IP['dpdir']+IP['analoguerxnsfiltfilename']+'.pickle')
        if IP['showresults']:
            IP['analoguerxnsfilt']=analoguerxnsfilt
    analoguerxns_updated=addspeciesdata(analoguerxnsfilt,IP['substancesource'],includesolv=IP['includesolv'],
                                       ncpus=IP['ncpus'],SQL=IP['SQL'],reaxys_update=IP['reaxys_update'],hc_Dict=IP['hc_Dict'],
                                       hc_rct=IP['hc_rct'],refanaloguerxns=IP['refanaloguerxns_updated'])
    if IP['writetofile']:
          pd.to_pickle(analoguerxns_updated,IP['dpdir']+IP['analoguerxns_updatedfilename']+'.pickle')
    if IP['showresults']:
          IP['analoguerxns_updated']=analoguerxns_updated
    print('Step 5 complete')
    return IP
            
def step6(IP):
    analoguerxns_updated=IP['analoguerxns_updated']
    analoguerxnsbal=IP['analoguerxnsbal']
    if analoguerxnsbal is None:
        if analoguerxns_updated is None:
            analoguerxns_updated=pd.read_pickle(IP['dpdir']+IP['analoguerxns_updatedfilename']+'.pickle')
        balrxns,analoguerxnsbal=balance_analogue_(analoguerxns_updated,refbalrxns=IP['refbalrxns'],reaxys_update=IP['reaxys_update'],
                                                 includesolv=IP['includesolv'],helpprod=IP['helpprod'],helpreact=IP['helpreact'],
                                                 addrctonly=IP['addrctonly'],ignoreH=IP['ignoreH'],ncpus=IP['ncpus'])
        if IP['writetofile']:
            writepickle(balrxns,IP['dpdir'],IP['balrxnsfilename'])
            pd.to_pickle(analoguerxnsbal,IP['dpdir']+IP['analoguerxnsbalfilename']+'.pickle')
        if IP['showresults']:
            IP['analoguerxnsbal']=analoguerxnsbal
    query='Invalid'
    if IP['removeLHSinsufficient']:
        query+='|LHS species insufficient'
    if IP['removesmilesdiscrepancy']:
        query+='|discrepancy'
    if IP['removemappingerror']:
        query+='|mapping error'
    if IP['removechargeimbalance']:
        query+='|Charge imbalance'
    if IP['removehcarriers']:
        analoguerxnsbalfilt=analoguerxnsbal.loc[(~analoguerxnsbal.msg.str.contains(query,case=False,na=False)) & (~analoguerxnsbal.msg.str.contains('Hydrogen carriers',na=False))]
    else:
        analoguerxnsbalfilt=analoguerxnsbal.loc[(~analoguerxnsbal.msg.str.contains(query,case=False,na=False))]
    if IP['writetofile']:
        pd.to_pickle(analoguerxnsbalfilt,IP['dpdir']+IP['analoguerxnsbalfiltfilename']+'.pickle')
    if IP['showresults']:
        IP['analoguerxnsbalfilt']=analoguernxsbalfilt
    print('Step 6 complete')
    return IP

def step7(IP):
#     breakpoint()
    analoguerxnsbalfilt=IP['analoguerxnsbal']
    analoguerxnsmapped=IP['analoguerxnsmapped']
    analoguerxnsparsed=IP['analoguerxnsparsed']
    analoguerxnsparsedfilt=IP['analoguerxnsparsedfilt']
    analoguerxnsassigned=IP['analoguerxnsassigned']
    if analoguerxnsassigned is None:
        if analoguerxnsparsedfilt is None:
            if analoguerxnsparsed is None:
                if analoguerxnsmapped is None:
                    if analoguerxnsbalfilt is None:
                        analoguerxnsbalfilt=pd.read_pickle(IP['dpdir']+IP['analoguerxnsbalfiltfilename']+'.pickle')
                    analoguerxnsmapped=map_rxns(analoguerxnsbalfilt,refmappedrxns=IP['refmappedrxns'],reaxys_update=IP['reaxys_update'],
                                               ncpus=IP['ncpus'])
                    if IP['writetofile']:
                        pd.to_pickle(analoguerxnsmapped,IP['dpdir']+IP['analoguerxnsmappedfilename']+'.pickle')
                    IP['analoguerxnsmapped']=analoguerxnsmapped
                analoguerxnsmappedfilt=analoguerxnsmapped.loc[analoguerxnsmapped.mapped_rxn!='Error']
                analoguerxnsparsed=checkrxns(analoguerxnsmappedfilt,refparsedrxns=IP['refparsedrxns'],reaxys_update=IP['reaxys_update'],
                                            ncpus=IP['ncpus'])
                if IP['writetofile']:
                    pd.to_pickle(analoguerxnsparsed,IP['dpdir']+IP['analoguerxnsparsedfilename']+'.pickle')
                if IP['showresults']:
                    IP['analoguerxnsparsed']=analoguerxnsparsed
            if IP['removesmilesdiscrepancy']:
                analoguerxnsparsedfilt=analoguerxnsparsed.loc[~analoguerxnsparsed.msg1.str.contains('discrepancy',case=False,na=False)]
            else:
                analoguerxnsparsedfilt=analoguerxnsparsed
            changedrxns=analoguerxnsparsedfilt.loc[(analoguerxnsparsedfilt.msg1.str.contains('Unmapped')) | (analoguerxnsparsedfilt.msg1.str.contains('unmapped'))]
            analoguerxns_updated=IP['analoguerxns_updated']
            if analoguerxns_updated is None:
                analoguerxns_updated=pd.read_pickle(IP['dpdir']+IP['analoguerxns_updatedfilename']+'.pickle')
            if IP['helpprod']:
                hc_prod=IP['hc_Dict']
            else:
                hc_prod={}
            changedrxns=updaterxns(changedrxns,hc_prod=hc_prod,analoguerxns=analoguerxns_updated,ncpus=IP['ncpus'])
            analoguerxnsparsedfilt.update(changedrxns[['mapped_rxn','confidence','balrxnsmiles','msg','LHS','RHS','hcrct','hcprod','LHSdata','RHSdata','msg1']])
            if IP['removemappingerror']:
                analoguerxnsparsedfilt=analoguerxnsparsedfilt.loc[analoguerxnsparsedfilt.mapped_rxn!='Error']
            if IP['removesmilesdiscrepancy']:
                analoguerxnsparsedfilt=analoguerxnsparsedfilt.loc[~analoguerxnsparsedfilt.msg1.str.contains('discrepancy',case=False,na=False)]
            if IP['removeRHSinsufficient']: #WILL REMOVE MOST REACTIONS, limitation of current approach (help compounds insufficient)
                analoguerxnsparsedfilt=analoguerxnsparsedfilt.loc[~analoguerxnsparsedfilt.msg.str.contains('RHS species insufficient')]
            if IP['nomixtures']:
                analoguerxnsparsedfilt=analoguerxnsparsedfilt.loc[~analoguerxnsparsedfilt.msg1.str.contains('Mixture detected for LHS',case=False,na=False)]
            if IP['removemandatoryunmapped']:
                analoguerxnsparsedfilt=analoguerxnsparsedfilt.loc[~analoguerxnsparsedfilt.msg1.str.contains('Mandatory',case=False,na=False)]
            if IP['writetofile']:
                pd.to_pickle(analoguerxnsparsedfilt,IP['dpdir']+IP['analoguerxnsparsedfiltfilename']+'.pickle')
            if IP['showresults']:
                IP['analoguerxnsparsedfilt']=analoguerxnsparsedfilt
        fragdict=IP['fragdict']
        if fragdict is None:
            fragdict=openpickle(IP['dmdir']+IP['fragdictfilename']+'.pickle')
        analoguerxnsassigned=assignfrags(analoguerxnsparsedfilt,fragdict,ncpus=IP['ncpus'])
        if IP['writetofile']:
            pd.to_pickle(analoguerxnsassigned,IP['dpdir']+IP['analoguerxnsassignedfilename']+'.pickle')
        if IP['showresults']:
            IP['analoguerxnsassigned']=analoguerxnsassigned
    analoguerxnsassignedfilt=analoguerxnsassigned.loc[~(analoguerxnsassigned.msg2.str.contains('not analogue',case=False,na=False))]
    if IP['writetofile']:
        pd.to_pickle(analoguerxnsassignedfilt,IP['dpdir']+IP['analoguerxnsassignedfiltfilename']+'.pickle')
    if IP['showresults']:
        IP['analoguerxnsassignedfilt']=analoguerxnsassigned
    print('Step 7 complete')
    return IP
           
def step8(IP):
    analoguerxnsassignedfilt=IP['analoguerxnsassignedfilt']
    analoguerxnscent=IP['analoguerxnscent']
    analoguerxnscentfilt=IP['analoguerxnscentfilt']
    analoguerxnsvalid=IP['analoguerxnsvalid']
    if analoguerxnsvalid is None:
        if analoguerxnscentfilt is None:
            if analoguerxnscent is None:
                if analoguerxnsassignedfilt is None:
                    analoguerxnsassignedfilt=pd.read_pickle(IP['dpdir']+IP['analoguerxnsassignedfiltfilename']+'.pickle')
                analoguerxnscent=reactioncenter(analoguerxnsassignedfilt,ncpus=IP['ncpus'])
                if IP['writetofile']:
                    pd.to_pickle(analoguerxnscent,IP['dpdir']+IP['analoguerxnscentfilename']+'.pickle')
                if IP['showresults']:
                    IP['analoguerxnscent']=analoguerxnscent
            analoguerxnscentfilt=analoguerxnscent.loc[analoguerxnscent.rxncenter==True]
            if IP['writetofile']:
                pd.to_pickle(analoguerxnscentfilt,IP['dpdir']+IP['analoguerxnscentfiltfilename']+'.pickle')
            if IP['showresults']:
                IP['analoguerxnscentfilt']=analoguerxnscentfilt
        analoguerxnsvalid=validreactioncenter(analoguerxnscentfilt,ncpus=IP['ncpus'])
        if IP['writetofile']:
            pd.to_pickle(analoguerxnsvalid,IP['dpdir']+IP['analoguerxnsvalidfilename']+'.pickle')
        if IP['showresults']:
            IP['analoguerxnsvalid']=analoguerxnsvalid
    if IP['removeoutfrags']:
        analoguerxnsfinal=analoguerxnsvalid.loc[~(analoguerxnsvalid['outfrag'].astype(bool))]
    else:
        analoguerxnsfinal=analoguerxnsvalid
    if IP['writetofile']:
        pd.to_pickle(analoguerxnsfinal,IP['dpdir']+IP['analoguerxnsfinalfilename']+'.pickle')
    if IP['showresults']:
        IP['analoguerxnsfinal']=analoguerxnsfinal
    print('Step 8 complete')
    return IP
   
def step9(IP):
    analoguerxnsfinal=IP['analoguerxnsfinal']
    analoguerxnstempl=IP['analoguerxnstempl']
    if analoguerxnstempl is None:
        if analoguerxnsfinal is None:
            analoguerxnsfinal=pd.read_pickle(IP['dpdir']+IP['analoguerxnsfinalfilename']+'.pickle')
        analoguerxnstempl=gentemplate(analoguerxnsfinal,ncpus=IP['ncpus'],specificity=IP['specificity'])
        if IP['writetofile']:
            pd.to_pickle(analoguerxnstempl,IP['ipdir']+IP['analoguerxnstemplfilename']+'.pickle')
        if IP['showresults']:
            IP['analoguerxnstempl']=analoguerxnstempl
    analoguerxnstemplfilt=analoguerxnstempl.loc[analoguerxnstempl.template!='Error']
    if IP['onlyvalidtempl']:
        analoguerxnstemplfilt=analoguerxnstempl.loc[analoguerxnstempl.msg4=='Valid']
    else:
        if IP['removefarfg']:
            analoguerxnstemplfilt=analoguerxnstemplfilt.loc[~analoguerxnstemplfilt.msg4.str.contains('far away',case=False,na=False)]
        if IP['removeunusedprod']:
            analoguerxnstemplfilt=analoguerxnstemplfilt.loc[~analoguerxnstemplfilt.msg4.str.contains('not produced',case=False,na=False)]      
    if IP['writetofile']:
        pd.to_pickle(analoguerxnstemplfilt,IP['ipdir']+IP['analoguerxnstemplfiltfilename']+'.pickle')
    if IP['showresults']:
        IP['analoguerxnstemplfilt']=analoguerxnstemplfilt
    print('Step 9 complete')
    return IP

def step10(IP):
    analoguerxnstemplfilt=IP['analoguerxnstemplfilt']
    inputquery=IP['inputquery']
    analoguerxnsimp=IP['analoguerxnsimp']
    if analoguerxnsimp is None:
        if analoguerxnstemplfilt is None:
            analoguerxnstemplfilt=pd.read_pickle(IP['ipdir']+IP['analoguerxnstemplfiltfilename']+'.pickle')
        if inputquery is None:
            inputquery=openpickle(IP['dmdir']+IP['iqfilename']+'.pickle')
        analoguerxnsimp=applytemplate(analoguerxnstemplfilt,inputquery,ncpus=IP['ncpus'])
        if IP['writetofile']:
            pd.to_pickle(analoguerxnsimp,IP['ipdir']+IP['analoguerxnsimpfilename']+'.pickle')
        if IP['showresults']:
            IP['analoguerxnsimp']=analoguerxnsimp
    analoguerxnsimpfilt=analoguerxnsimp.loc[analoguerxnsimp.msg5=='Valid']
    analoguerxnsimpfilt=removeduplicates(analoguerxnsimpfilt,ncpus=IP['ncpus'],restart=False)
    if IP['writetofile']:
        pd.to_pickle(analoguerxnsimpfilt,IP['ipdir']+IP['analoguerxnsimpfiltfilename']+'.pickle')
    if IP['showresults']:
        IP['analoguerxnsimpfilt']=analoguerxnsimpfilt
    print('Step 10 complete')
    return IP
    
        
def step11(IP):
    analoguerxnsimpfilt=IP['analoguerxnsimpfilt']
    analoguerxnsfilt=IP['analoguerxnsfilt']
    analoguerxns_updated=IP['analoguerxns_updated']
    inputquery=IP['inputquery']
    impfinal=IP['impfinal']
    if impfinal is None:
        if analoguerxnsimpfilt is None:
            analoguerxnsimpfilt=pd.read_pickle(IP['ipdir']+IP['analoguerxnsimpfiltfilename']+'.pickle')
        if analoguerxnsfilt is None:
            analoguerxnsfilt=pd.read_pickle(IP['dpdir']+IP['analoguerxnsfiltfilename']+'.pickle')
        if analoguerxns_updated is None:
            analoguerxns_updated=pd.read_pickle(IP['dpdir']+IP['analoguerxns_updatedfilename']+'.pickle')
        if inputquery is None:
            inputquery=openpickle(IP['dmdir']+IP['iqfilename']+'.pickle')
        impfinal=cleanimpurities(analoguerxnsimpfilt,analoguerxnsfilt,analoguerxns_updated,inputquery,
                                 includesolv=IP['includesolv'],reaxys_update=IP['reaxys_update'],hc_prod=IP['hc_Dict'])
        if IP['writetofile']:
            pd.to_pickle(impfinal,IP['ipdir']+IP['impfinalfilename']+'.pickle')
        if IP['showresults']:
            IP['impfinal']=impfinal
    impfinalfilt=impfinal.loc[(impfinal.msg6=='Valid')|(impfinal.msg6.str.contains('Query'))|(impfinal.msg6.str.contains('Self'))]
    if IP['writetofile']:
        pd.to_pickle(impfinalfilt,IP['ipdir']+IP['impfinalfiltfilename']+'.pickle')
    if IP['showresults']:
        IP['impfinalfilt']=impfinalfilt
    print('Step 11 complete')
    return IP

def step12(IP):
    impfinalfilt=IP['impfinalfilt']
    combinedpool=IP['combinedpool']
    combinedpoolex=IP['combinedpoolex']
    fragdict=IP['fragdict']
    if impfinalfilt is None:
        impfinalfilt=pd.read_pickle(IP['ipdir']+IP['impfinalfiltfilename']+'.pickle')
    if combinedpool is None:
        combinedpool=openpickle(IP['dmdir']+IP['combinedpoolfilename']+'.pickle')
    if combinedpoolex is None:
        combinedpoolex=openpickle(IP['dmdir']+IP['combinedpoolexfilename']+'.pickle')
    if fragdict is None:
        fragdict=openpickle(IP['dmdir']+IP['fragdictfilename']+'.pickle')
    catpool=combinedpoolex-combinedpool
    impfinalfilt=updatecatalyst(impfinalfilt,catpool,ncpus=IP['ncpus'])
    impfinalfilt0=copy.deepcopy(impfinalfilt)
    impfinalfilt=relevance_score_morgan(impfinalfilt,fragdict,includereagents=True,ncpus=IP['ncpus']) # setting up morgan fingerprints
    impfinalfilt_rctonly=relevance_score_morgan(impfinalfilt0,fragdict,includereagents=False,ncpus=IP['ncpus'])
    impfinalfilt['Relevance_morgan_rctonly']=impfinalfilt_rctonly['Relevance_morgan']
    impfinalfilt=standardize(impfinalfilt,reaxys_update=IP['reaxys_update'],ncpus=IP['ncpus'],restart=True)
    if IP['writetofile']:
        pd.to_pickle(impfinalfilt,IP['irdir']+IP['impfinalfiltfilename']+'.pickle')
    if IP['showresults']:
        IP['impfinalfilt']=impfinalfilt
    print('Step 12 complete')
    return IP
   
    
def step13(IP,conditions=[],catalyst=[]):
    impfinalfilt=IP['impfinalfilt']
    impfinalfilt5=IP['impfinalfilt5']
    if impfinalfilt5 is None:
        if impfinalfilt is None:
            impfinalfilt=pd.read_pickle(IP['irdir']+IP['impfinalfiltfilename']+'.pickle')
        impfinalfilt2=conditionfilter(impfinalfilt,conditions=conditions) #Remove invalid conditions
        impfinalfilt3=impfinalfilt2.loc[~impfinalfilt2.msg6.str.contains('Self')] #Remove suspect self-reactions
        impfinalfilt4=catfilter(impfinalfilt3,IP['substancesource'],IP['unresolveddir'],catalyst=catalyst,useray=False)
        impfinalfilt5=copy.deepcopy(impfinalfilt4)
        impfinalfilt5['Temperature']=impfinalfilt5.apply(updatetemp,axis=1)
        if IP['writetofile']:
            pd.to_pickle(impfinalfilt2,IP['irdir']+IP['impfinalfilt2filename']+'.pickle')
            pd.to_pickle(impfinalfilt3,IP['irdir']+IP['impfinalfilt3filename']+'.pickle')
            pd.to_pickle(impfinalfilt4,IP['irdir']+IP['impfinalfilt4filename']+'.pickle')
            pd.to_pickle(impfinalfilt5,IP['irdir']+IP['impfinalfilt5filename']+'.pickle')
        if IP['showresults']:
            IP['impfinalfilt5']=impfinalfilt5
            IP['impfinalfilt2']=impfinalfilt2
            IP['impfinalfilt3']=impfinalfilt3
            IP['impfinalfilt4']=impfinalfilt4
    print('Step 13 complete')
    return IP
        
        
def step14(IP,Trange=None):
    impfinalfilt5=IP['impfinalfilt5']
    if impfinalfilt5 is None:
        impfinalfilt5=pd.read_pickle(IP['irdir']+IP['impfinalfilt5filename']+'.pickle')
        
    if IP['includefraginfo']: #Streamlining needed, slow code
#         breakpoint()
        inputquery=IP['inputquery']
        if inputquery is None:
            inputquery=openpickle(IP['dmdir']+IP['iqfilename']+'.pickle')
        # Adding fragment information
        initray(num_cpus=IP['ncpus'])
        impfinalfilt5dis=mpd.DataFrame(impfinalfilt5)
        fraginfo=impfinalfilt5dis.apply(popfragmentsrow,axis=1,result_type='reduce')
        fraginfo=pd.DataFrame(data=fraginfo.values,index=fraginfo.index,columns=['Fraginfo'])
        impfinalfilt5['Fraginfo']=fraginfo 
        queryfraginfoa=popfragments(inputquery['species'])
        initray(num_cpus=IP['ncpus'])
        impfinalfilt5dis=mpd.DataFrame(impfinalfilt5)
        fragcompare=impfinalfilt5dis.apply(fragcomparerow,queryfraginfoa=queryfraginfoa,axis=1,result_type='reduce')
        fragcompare=pd.DataFrame(data=fragcompare.values,index=fragcompare.index,columns=['Fragcompare'])
        impfinalfilt5['Fragcompare']=fragcompare
        impfinalfilt5dis=mpd.DataFrame(impfinalfilt5)
        fraglists=impfinalfilt5dis.apply(summarizefrags,axis=1,result_type='reduce')
        fraglists=pd.DataFrame(data=fraglists.tolist(),index=fraglists.index,columns=['commonfrags1','uniquefragsa1',
                                                                                  'uniquefragsq1','commonfrags2','uniquefragsa2',
                                                                                  'uniquefragsq2','commonfrags3',
                                                                                  'uniquefragsa3','uniquefragsq3'])

        impfinalfilt5[['commonfrags1','uniquefragsa1','uniquefragsq1','commonfrags2','uniquefragsa2','uniquefragsq2',
                   'commonfrags3','uniquefragsa3','uniquefragsq3']]=fraglists
        if IP['writetofile']:
            pd.to_pickle(impfinalfilt5,IP['irdir']+IP['impfinalfilt5filename']+'.pickle')
        if IP['showresults']:
            IP['impfinalfilt5']=impfinalfilt5
            
    summary=impfinalfilt5.groupby(['impurityrxn','querycompds','impurities'])['impurityrxn'].count()
    summary=pd.DataFrame(data=summary.values,index=summary.index,columns=['Hits'])
#     summary['Max_relevance_top']=impfinalfilt5.groupby(['impurityrxn','querycompds','impurities']).Relevance.max().round(2)
    summary['Max_relevance_morg']=impfinalfilt5.groupby(['impurityrxn','querycompds','impurities']).Relevance_morgan.max().round(2)
    summary['Max_relevance_morg_rctonly']=impfinalfilt5.groupby(['impurityrxn','querycompds','impurities']).Relevance_morgan_rctonly.max().round(2)
    summary=summary.sort_values(by=['Max_relevance_morg','Hits'],ascending=False)
    mainrxn=summary.index.get_level_values(0)[0]

    summary2=rank_impuritiest(summary,mainrxn,impfinalfilt5,trange=Trange,relquantile=IP['relquantile'],
                             uppertquantile=IP['uppertquantile'],lowertquantile=IP['lowertquantile'])
    summary2=pd.DataFrame(summary2)
    summary3=summary2.loc[summary2.Hits_tfiltered>1]
    summary3=summary3.sort_values(by=['Max relevance_tfiltered','Hits_tfiltered','Max relevance_missing','Hits_missing'],ascending=False)
    if IP['writetofile']:
        pd.to_pickle(summary,IP['irdir']+IP['summaryfilename']+'.pickle')
        pd.to_pickle(summary2,IP['irdir']+IP['summary2filename']+'.pickle')
        pd.to_pickle(summary3,IP['irdir']+IP['summary3filename']+'.pickle')
    if IP['showresults']:
        IP['summary']=impfinalfilt5
        IP['summary2']=impfinalfilt2
        IP['summary3']=impfinalfilt3
    print('Step 14 complete')
    return IP

    
    
    
    
    
    
    
    
# impfinalfilt=impfinal.loc[(impfinal.msg6=='Valid')|(impfinal.msg6.str.contains('Query'))|(impfinal.msg6.str.contains('Self'))]
#         impfinalfilt=impfinalfilt.loc[~impfinalfilt.hcrct.astype(bool)] #Remove help reactants
#         print(str(round(len(impfinalfilt.index)/len(impfinal.index)*100,2))+'% reactions remaining with valid, sensible impurities')
#         catpool=combinedpoolex-combinedpool
#         impfinalfilt=updatecatalyst(impfinalfilt,catpool)
#         #Topological fingerprints
#         impfinalfilt=relevance(impfinalfilt,inputquery)
#         impfinalfilt0=copy.deepcopy(impfinalfilt)
#         # Morgan fingerprints
#         impfinalfilt=relevance_score_morgan(impfinalfilt,querydict,includereagents=True) # setting up morgan fingerprints
#         impfinalfilt_rctonly=relevance_score_morgan(impfinalfilt0,querydict,includereagents=False)
#         impfinalfilt['Relevance_morgan_rctonly']=impfinalfilt_rctonly['Relevance_morgan']
#         impfinalfilt=standardize(impfinalfilt,reaxys_update=reaxys_update)
#         pd.to_pickle(impfinalfilt,rxncenterscreendir+'impurities_final_filt.pickle')
#         print('Step 15 done')    
    
#     pd.to_pickle(analoguerxnstempl,ipdir+'analoguerxnstempl'+'_'+str(chunknum)+'.pickle')

#     analoguerxnsbal2=analoguerxnsbal5.loc[~(analoguerxnsbal5.msg.str.contains('LHS species insufficient',na=False)) & ~(analoguerxnsbal5.msg.str.contains('Invalid',na=False)) & ~(analoguerxnsbal5.msg.str.contains('mapping error',case=False,na=False)) & ~(analoguerxnsbal5.msg.str.contains('discrepancy',case=False,na=False))] #& ~(analoguerxnsbal5.msg.str.contains('Mandatory'))
#     pd.to_pickle(analoguerxnsbal5,dpdir+'analoguerxnsbal'+'_'+str(chunknum)+'.pickle')
#     pd.to_pickle(analoguerxnsbal2,dpdir+'analoguerxnsbalfilt'+'_'+str(chunknum)+'.pickle')

#     analoguerxnsmapped=map_rxns(analoguerxnsbal2)
#     pd.to_pickle(analoguerxnsmapped,dpdir+'analoguerxnsmapped'+'_'+str(chunknum)+'.pickle')
#     analoguerxnsmappedfilt=analoguerxnsmapped[analoguerxnsmapped.mapped_rxn!='Error']

#     analoguerxnsparsed=checkrxns(analoguerxnsmappedfilt)
#     pd.to_pickle(analoguerxnsparsed,dpdir+'analoguerxnsparsed'+'_'+str(chunknum)+'.pickle')
#     #     analoguerxnsparsed=pd.read_pickle(dpdir+'analoguerxnsparsed'+'_'+str(chunknum)+'.pickle')
#     analoguerxnsparsedfilt=analoguerxnsparsed.loc[~analoguerxnsparsed.msg1.str.contains('discrepancy',case=False,na=False)]
#     changedrxns=analoguerxnsparsedfilt.loc[(analoguerxnsparsedfilt.msg1.str.contains('Unmapped')) | (analoguerxnsparsedfilt.msg1.str.contains('unmapped'))]
#     changedrxns=updaterxns(changedrxns,hc_prod={},analoguerxns=analoguerxns3)
#     analoguerxnsparsedfilt.update(changedrxns[['mapped_rxn','confidence','balrxnsmiles','msg','LHS','RHS','hcrct','hcprod','LHSdata','RHSdata','msg1']])
#     analoguerxnsparsedfilt=analoguerxnsparsedfilt.loc[(analoguerxnsparsedfilt.mapped_rxn!='Error') & ~(analoguerxnsparsedfilt.msg1.str.contains('discrepancy',case=False,na=False))]

#     fragdict=openpickle(dmdir+'fragdict.pickle')
#     analoguerxnsassigned=assignfrags(analoguerxnsparsedfilt,fragdict,ncpus=8)
#     pd.to_pickle(analoguerxnsassigned,dpdir+'analoguerxnsassigned'+'_'+str(chunknum)+'.pickle')
#     analoguerxnsassignedfilt=analoguerxnsassigned.loc[~(analoguerxnsassigned.msg2.str.contains('not analogue',case=False,na=False))]

#     analoguerxnscent=reactioncenter(analoguerxnsassignedfilt)
#     pd.to_pickle(analoguerxnscent,dpdir+'analoguerxnscent'+'_'+str(chunknum)+'.pickle')
#     analoguerxnscentfilt=analoguerxnscent.loc[analoguerxnscent.rxncenter==True]

#     analoguerxnsvalid=validreactioncenter(analoguerxnscentfilt,ncpus=16)
#     pd.to_pickle(analoguerxnsvalid,dpdir+'analoguerxnsvalid'+'_'+str(chunknum)+'.pickle')
#     analoguerxnsfinal=analoguerxnsvalid.loc[~(analoguerxnsvalid['outfrag'].astype(bool))]

#     analoguerxnstempl=gentemplate(analoguerxnsfinal)
#     pd.to_pickle(analoguerxnstempl,ipdir+'analoguerxnstempl'+'_'+str(chunknum)+'.pickle')

#     reactionfilter=analoguerxnstempl.apply(matchtemplaterow,samplerxn='c[H].Cl>>cCl',axis=1,result_type='reduce')
#     reactionfilterlist=reactionfilter[reactionfilter.values==True].index
#     reaction2=analoguerxnstempl[analoguerxnstempl.index.isin(reactionfilterlist)]
#     pd.to_pickle(reaction2,ipdir+'chlorination'+'_'+str(chunknum)+'.pickle')
#     lb=ub
    
        

# %load ./BalanceRxns.py
import rdkit
from MainFunctions import CustomError,getfragments, getcompdict, molfromsmiles
from chempy import balance_stoichiometry
import copy
from collections import Counter
from decimal import Decimal, ROUND_HALF_UP
from rdkit import Chem #Importing RDKit
from rdkit.Chem import rdChemReactions #Reaction processing
# from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from MapRxns import maprxn


def balance_analogue_(analoguerxns,refbalrxns=None,coefflim=6,reaxys_update=True,includesolv=False,
                     usemapper=True,helpprod=True,helpreact=True,ncpus=16):
    '''
    Applies balance_analogue across each row of a given dataframe
    
    '''
#     breakpoint()
    index=analoguerxns.index.name
    if index and index!='ReactionID':
        analoguerxns.reset_index(inplace=True).set_index('ReactionID',inplace=True)
    if refbalrxns is not None:
        if index:
            commonids=set(refbalrxns.ReactionID).intersection(set(analoguerxns.index))
            commondf=refbalrxns[refbalrxns.ReactionID.isin(commonids)]
            analoguerxns=analoguerxns[~analoguerxns.index.isin(commonids)]
        else:
            commonids=set(refbalrxns.ReactionID).intersection(set(analoguerxns.ReactionID))
            commondf=refbalrxns[refbalrxns.ReactionID.isin(commonids)]
            analoguerxns=analoguerxns[~analoguerxns.ReactionID.isin(commonids)]
        if analoguerxns.empty: #All reactions found
            return commondf
    initray(num_cpus=ncpus)
    analoguerxnsdis=mpd.DataFrame(analoguerxns)
    balrxns=analoguerxnsdis.apply(balance_analogue,coefflim=6,includesolv=includesolv,usemapper=usemapper,
                               helpprod=helpprod,helpreact=helpreact,axis=1,result_type='reduce')
    balrxns=pd.Series(data=balrxns.values,index=balrxns.index) #Optional convert modin back to pandas
    analoguerxnsbal=pd.DataFrame(balrxns,columns=['rxnsmiles'])
    analoguerxnsbal[['rxnsmiles0', 'balrxnsmiles','msg','LHS','RHS','hcrct','hcprod','LHSdata','RHSdata']] = pd.DataFrame(analoguerxnsbal['rxnsmiles'].tolist(), index=analoguerxnsbal.index)
    if reaxys_update:
        analoguerxnsbal[['ReactionID','NumRefs','NumSteps','NumStages']]=analoguerxns[['ReactionID','NumRefs','NumSteps','NumStages']]
        cols=['ReactionID','NumRefs','NumSteps','NumStages','rxnsmiles0','balrxnsmiles','msg','LHS','RHS','hcrct','hcprod','LHSdata','RHSdata']
        analoguerxnsbal=analoguerxnsbal[cols]
    else:
        analoguerxnsbal[['NumRefs','NumSteps']]=analoguerxns[['NumRefs','NumSteps']]
        cols=['NumRefs','NumSteps','rxnsmiles0','balrxnsmiles','msg','LHS','RHS','hcrct','hcprod','LHSdata','RHSdata']
        analoguerxnsbal=analoguerxnsbal[cols]  
    return balrxns,analoguerxnsbal

def balance_analogue(row,basic=True,balance=True,coefflim=6,includesolv=False,
                     usemapper=True,helpprod=True,helpreact=True): #More reliable
    '''
    Applies balancerxn function across a given dataframe row
    
    
    '''
    Rdata={}
    Pdata={}
    Rgtdata={}
    Solvdata={}
    Rgtid=[]
    hc_prod={}
    hc_react={}
    Rdata=copy.deepcopy(row['Rdata'])
    Pdata=copy.deepcopy(row['Pdata'])
    Rgtdata=copy.deepcopy(row['Rgtdata'])
    if row['ReagentID']!='NaN':
        Rgtid=copy.deepcopy(row['ReagentID'])
    if helpprod:
        hc_prod=copy.deepcopy(row['hc_prod'])
    if helpreact:
        hc_react=copy.deepcopy(row['hc_react'])
    if includesolv:
        Solvdata=copy.deepcopy(row['Solvdata'])
    if type(Rdata)!=dict: #Error or compound invalid
        if basic and balance:
            return 'Error','Error',Rdata,'NaN','NaN',[],[],'NaN','NaN'
        else:
            return 'Error',Rdata,'NaN','NaN',[],[],'NaN','NaN' 
    elif type(Pdata)!=dict: #Error or compound invalid
        if basic and balance:
            return 'Error','Error',Pdata,'NaN','NaN',[],[],'NaN','NaN'
        else:
            return 'Error',Pdata,'NaN','NaN',[],[],'NaN','NaN' 
    rxnsmiles0='>>'.join([getfragments([Rdata[r]['smiles'] for r in Rdata],smiles=True),
                          getfragments([Pdata[p]['smiles'] for p in Pdata],smiles=True)])
    if basic and not balance:
        return rxnsmiles0,list(Rdata.keys()),list(Pdata.keys())
    else:
#         return False
#         breakpoint()
        return balancerxn(Rdata,Pdata,Rgtdata=Rgtdata,Solvdata=Solvdata,rxnsmiles0=rxnsmiles0,usemapper=usemapper,coefflim=coefflim,hc_prod=hc_prod,hc_react=hc_react)

def balancerxn(Rdata,Pdata,Rgtdata={},Solvdata={},rxnsmiles0=None,first=True,usemapper=False,
               addedspecies=[],addedhc=[],hc_prod={},hc_react={},coefflim=6,msg=''):
    '''
    Balances reactions given reactant species information (Rdata) and product species information (Pdata)
    
    Rgtdata is optional and refers to reagent species information
    Solvdata is optional and refers to solvent species information
    rxnsmiles0 refers to the reaction SMILES string as represented in Reaxys
    first is True/False depending on whether it is the first time running through the code
    usemapper is True if IBM RXN mapper is used to decide between possible reactants on LHS
    addedspecies refers to new species added to the LHS
    addedhc refers to small species (help reactants) added to the LHS
    hc_prod is optional and is a dictionary of small species (help compounds) to use for balancing the RHS
    hc_react is optional and is a dictionary of small species (help reactants) to use for balancing the LHS
    coefflim refers to the maximum allowed stoichiometric coefficient after balancing
    msg involves any warning messages or updates from the code
    
    '''
    
    if Solvdata: #Add sovents
        Rgtdata={**Rgtdata,**Solvdata} #Newly added
        
    #%% Initialize added species/addedhc (help compound) and add small species
    if first:
        addedspecies=[]
        addedhc=[]
        if Rgtdata:
#             smallspecies=[spec for spec in Rgtdata if Rgtdata[spec]['formula'] in ['H2','O2','H2O','OH2','HOH'] if spec not in Solvdata]
            smallspecies=[spec for spec in Rgtdata if Rgtdata[spec]['formula'] in ['H2','O2','H2O','OH2','HOH']]
            if smallspecies:
                Rdata.update({spec:Rgtdata[spec] for spec in smallspecies})
                addedspecies+=smallspecies
                addedsmallspec=True
    Rcount=sum([Counter(Rdata[ID]['atomdict']) for ID in Rdata for _ in range(Rdata[ID]['count'])],start=Counter()) #Sum of atom counts/types on LHS
    Rcharge=sum([Rdata[ID]['charge'] for ID in Rdata for _ in range(Rdata[ID]['count'])])
    Pcount=sum([Counter(Pdata[ID]['atomdict']) for ID in Pdata for _ in range(Pdata[ID]['count'])],start=Counter()) #Sum of atom counts/types on RHS
    Pcharge=sum([Pdata[ID]['charge'] for ID in Pdata for _ in range(Pdata[ID]['count'])])

    #%% String output handling
    addedstr=''
    if addedspecies:
        addedstr=' with species: '+(','.join([str(species) for species in addedspecies]))
    if addedhc:
        addedstr2=' with help reactant(s): '+(','.join([hc_react[species]['formula'] for species in addedhc]))
        if addedstr:
            addedstr=addedstr+', '+addedstr2
        else:
            addedstr=addedstr2
            
    #%% If reaction is balanced already
    if Rcount==Pcount and Rcharge==Pcharge:
        print('Reaction is fully balanced')
        if first:
            msg='Already balanced'
        elif msg:
            msg=msg+', '+'Balanced'
        else:
            msg='Balanced'
        if addedstr:
            msg+=addedstr
        return update_rxn(Rdata,Pdata,hcrct=addedhc,rxnsmiles0=rxnsmiles0,msg=msg)
        
#%% Otherwise take difference between atom type/count RHS and LHS

    rem=Counter() # Rem contains difference between product and reactant counters
    rem.update(Pcount)  #If atoms not balanced and same charge, attempt to balance. Can try if different charge but more tricky
    rem.subtract(Rcount) #Subtracting reactant atom type index from product atom type index
    postype={key:rem[key] for key in rem.keys() if rem[key]>0} #Finding only positive keys. Note that if counter is positive this means extra molecules need to be added on LHS (eg. reagent).
    negtype={key:abs(rem[key]) for key in rem.keys() if rem[key]<0} #Finding only negative keys. If counter is negative this means extra molecules need to be added on RHS (eg. help compounds)

#     breakpoint()
    if 'hydrogen carriers' in msg:
        return update_rxn(Rdata,Pdata,hc_prod=hc_prod,hcrct=addedhc,rxnsmiles0=rxnsmiles0,msg=msg)
    
    
    elif postype: #Reactants, Reagents may be needed
        #%% Initializing variables
        candirxt=[]
        candirgt=[]
        candihc=[]
        matches=[]
        addedspecies_=[]
        addedhc_=[]
        #%% Get reactant match first
        candirxt=[rctid for rctid in Rdata if set(postype.keys()).issubset(set(Rdata[rctid]['atomdict'].keys()))]
        #%% Get reagent match
        candirgt=[rgtid for rgtid in Rgtdata if set(postype.keys()).issubset(set(Rgtdata[rgtid]['atomdict'].keys())) and rgtid not in candirxt]                        
        #%% Get help compound match
        if hc_react:
            candihc=[hcid for hcid in hc_react if set(postype.keys()).issubset(set(hc_react[hcid]['atomdict'].keys()))]
#         breakpoint() 
        if candirxt and not candirgt: #Only reactant matches
            Rdata,addedspecies_,msg_=resolvecandidates(postype,Rdata,Rdata,candirxt,validate=False)
        elif candirgt and not candirxt:
            Rdata,addedspecies_,msg_=resolvecandidates(postype,Rdata,Rgtdata,candirgt,validate=False)
        elif candirxt and candirgt:
            Rdata,addedspecies_,msg_=resolvecandidates(postype,Rdata,{**Rdata,**Rgtdata},candirxt+candirgt,validate=False)
        elif not candirxt and not candirgt:
            Rdata,addedspecies_,msg_=resolvecandidates(postype,Rdata,{**Rdata,**Rgtdata},list(Rdata.keys())+list(Rgtdata.keys()))
        elif candihc and not candirxt and not candirgt:
            Rdata,addedhc_,msg_=resolvecandidates(postype,Rdata,hc_react,candihc,validate=False)
        if 'hydrogen carriers' in msg_: #Multiple candidates for hydrogen carriers (mapper won't help)
            if msg:
                msg=msg+', '+msg_
            else:
                msg=msg_
            
        elif msg_!='Valid': #Atom surplus on RHS cannot be met by any LHS species
            if msg:
                msg=msg+', '+msg_+addedstr
            else:
                msg=msg_+addedstr
            if rxnsmiles0 is not None:
                return rxnsmiles0,'Error',msg,list(Rdata.keys()),list(Pdata.keys()),[],[],Rdata,Pdata
            else:
                return 'Error',msg,list(Rdata.keys()),list(Pdata.keys()),[],[],Rdata,Pdata 
        else:
            addedspecies+=addedspecies_
        return balancerxn(Rdata,Pdata,Rgtdata=Rgtdata,rxnsmiles0=rxnsmiles0,first=False,usemapper=usemapper,
                         addedspecies=addedspecies,addedhc=addedhc+addedhc_,hc_prod=hc_prod,
                         hc_react=hc_react,coefflim=coefflim,msg=msg) 
    
    elif negtype:
        if usemapper and len(set(addedspecies))>1: #More than one choice or added small species, let the mapper decide
            addedspecies_=[]
            addedhc_=[]
            rxnsmiles=buildrxn(Rdata,Pdata)
            mapped_rxn=maprxn([rxnsmiles])[0]
            if mapped_rxn=='Error':
                if msg:
                    msg=msg+', '+'Mapping error'+addedstr
                else:
                    msg='Mapping error'+addedstr
                if rxnsmiles0 is not None:
                    return update_rxn(Rdata,Pdata,hcrct=addedhc,rxnsmiles0=rxnsmiles0,msg=msg)
                else:
                    return update_rxn(Rdata,Pdata,hcrct=addedhc,msg=msg)
            else:
                mappedrxn=mapped_rxn.get('mapped_rxn')
                conf=mapped_rxn.get('confidence')
                rdrxn=rdChemReactions.ReactionFromSmarts(mappedrxn,useSmiles=True)
                cleanrxn=copy.copy(rdrxn)
                rdChemReactions.RemoveMappingNumbersFromReactions(cleanrxn)
#                 breakpoint()
                rmixtures={}
                LHSdata={}
                for ID,rct in enumerate(cleanrxn.GetReactants()):
                    foundmatch=False
                    mappedmol=rdrxn.GetReactants()[ID]
                    formula=rdkit.Chem.rdMolDescriptors.CalcMolFormula(rct)
                    if any([atom.HasProp('molAtomMapNumber') for atom in mappedmol.GetAtoms()]) or formula=='H2': #Confirmed, mapped reactant
                        rctsmiles=Chem.MolToSmiles(molfromsmiles(Chem.MolToSmiles(rct))) #Ensuring RDKit smiles
                        for ID0 in Rdata:
                            if '.' in Rdata[ID0]['smiles'] and Rdata[ID0]['smiles']!=rctsmiles and rctsmiles in Rdata[ID0]['smiles'].split('.'): #Mixture detected
                                foundmatch=True
                                if ID0 not in rmixtures:
                                    rmixtures.update({ID0:{rctsmiles:1}})
                                elif rctsmiles not in rmixtures[ID0]:
                                    rmixtures[ID0].update({rctsmiles:1})
                                else:
                                    rmixtures[ID0][rctsmiles]+=1
                                break
                            elif Rdata[ID0]['smiles']==rctsmiles:
                                foundmatch=True
                                if ID0 not in LHSdata.keys():
                                    LHSdata.update({ID0:copy.deepcopy(Rdata[ID0])})
                                    LHSdata[ID0]['count']=1
                                else:
                                    LHSdata[ID0]['count']+=1
                                break
                        if not foundmatch:
                            LHSdata=Rdata
                            if msg:
                                msg=msg+', '+'Smiles discrepancy for species'
                            else:
                                msg='Smiles discrepancy for species'
                            break
                if msg:
                    msg=msg+', '+'Mapper used'
                else:
                    msg='Mapper used'
                if rmixtures: #Mixture parsing
                    msg=msg+', '+'Mixture detected for species: '+','.join([str(ID) for ID in rmixtures])
                    for ID0 in rmixtures:
                        count=max([rmixtures[ID0][rctsmiles] for rctsmiles in rmixtures[ID0]])
                        LHSdata.update({ID0:Rdata[ID0]})
                        LHSdata[ID0]['count']=count
            
            if addedhc:
                refhc=Counter(addedhc)
                addedhc=[addedspec for addedspec in LHSdata for _ in range(min([LHSdata[addedspec]['count'],refhc[addedspec]])) if addedspec in addedhc]  
            refspec=Counter(addedspecies)
            addedspecies=[addedspec for addedspec in LHSdata for _ in range(min([LHSdata[addedspec]['count'],refspec[addedspec]])) if addedspec in addedspecies if addedspec not in addedhc]            
            return balancerxn(LHSdata,Pdata,Rgtdata=Rgtdata,rxnsmiles0=rxnsmiles0,first=False,usemapper=False,
                             addedspecies=addedspecies,addedhc=addedhc,hc_prod=hc_prod,hc_react=hc_react,coefflim=coefflim,
                             msg=msg)
        else:
            # Find match with Rdata first
            candidates=[spec for spec in Rdata if Rdata[spec]['atomdict']==negtype if spec in addedspecies]
            if not candidates:
                candidates=[spec for spec in Rdata if set(Rdata[spec]['atomdict'].keys())==set(negtype.keys()) if spec in addedspecies]
            if candidates:
                specremoved=False
                for candi in candidates:
                    matcharray=list(set(Counter({k:negtype[k]/Rdata[candi]['atomdict'][k] for k in negtype}).values()))
                    if all([mult.is_integer() for mult in matcharray]):
                        mult=min(matcharray)
                        newcount=int(Rdata[candi]['count']-mult)
                        if newcount<0:
                            newcount=0
                        if newcount==0:
                            del Rdata[candi]
                        else:
                            Rdata[candi]['count']=newcount
                        specremoved=True
                        break
                if specremoved:
                    if addedhc:
                        refhc=Counter(addedhc)
                        addedhc=[addedspec for addedspec in Rdata for _ in range(min([Rdata[addedspec]['count'],refhc[addedspec]])) if addedspec in addedhc]
                    refspec=Counter(addedspecies)
                    addedspecies=[addedspec for addedspec in Rdata for _ in range(min([Rdata[addedspec]['count'],refspec[addedspec]])) if addedspec in addedspecies if addedspec not in addedhc] 
                    if 'Mixture' in msg:
                        msglist=msg.split(',')
                        msglist2=copy.copy(msglist)
                        for i,msg_ in enumerate(msglist):
                            if 'Mixture' in msg_:
                                mixmsg=msg_.split(':')
                                refmixtures=mixmsg[1].split(',')
                                rmixtures={rmixture for rmixture in refmixtures if rmixture in Rdata}
                                if rmixtures:
                                    msglist2[i]=mixmsg[0]+','.join(rmixtures)
                                else:
                                    msglist2.remove(msglist2[i])
                                msg=', '.join(msglist2)
                    return balancerxn(Rdata,Pdata,Rgtdata=Rgtdata,rxnsmiles0=rxnsmiles0,first=False,usemapper=False,
                                  addedspecies=addedspecies,addedhc=addedhc,hc_prod=hc_prod,hc_react=hc_react,
                                  coefflim=coefflim,msg=msg)
                        
            # Then try help compounds
            if Rcharge==Pcharge:
                hc_prod={hcid:hc_prod[hcid] for hcid in hc_prod if hc_prod[hcid]['charge']==0} #Atom type for help compounds

            hc_prod2={hcid:hc_prod[hcid] for hcid in hc_prod if hc_prod[hcid]['atomdict']==negtype} #Exact match for 1 compound
            if not hc_prod2:
                hc_prod2={hcid:hc_prod[hcid] for hcid in hc_prod if hc_prod[hcid]['atomdict'].keys()==negtype.keys()} #Key match for 1 compound
            if hc_prod2:
                balbefore=False
            else:
                balbefore=True
            reac={}
            prod={}
            hcid=[]
            try:
                reac,prod,hcid,msg0=balance(Rdata,Pdata,hc_prod=hc_prod2,balbefore=balbefore,coefflim=coefflim,
                                            addedspecies=addedspecies,addedhc=addedhc,hc_react=hc_react)
                if msg:
                    msg=msg+', '+msg0
                else:
                    msg=msg0
            except Exception:
                if not balbefore:
                    return balancerxn(Rdata,Pdata,Rgtdata=Rgtdata,rxnsmiles0=rxnsmiles0,first=False,usemapper=False,
                                  addedspecies=addedspecies,addedhc=addedhc,hc_prod={},hc_react=hc_react,
                                  coefflim=coefflim,msg=msg)
                if msg: 
                    msg=msg+', '+'RHS species insufficient'+addedstr
                else:
                    msg='RHS species insufficient'+addedstr
                if Rcharge!=Pcharge:
                    msg=msg+', '+'Charge imbalance'
                return update_rxn(Rdata,Pdata,hcrct=addedhc,rxnsmiles0=rxnsmiles0,msg=msg)
            else:
                if Rcharge!=Pcharge:
                    if msg:
                        msg=msg+', '+'Charge imbalance'
                    else:
                        msg='Charge imbalance'
                return update_rxn(Rdata,Pdata,reac=reac,prod=prod,hcprod=hcid,hc_prod=hc_prod2,hcrct=addedhc,
                                  rxnsmiles0=rxnsmiles0,msg=msg)
    else:
        if Rcharge!=Pcharge:
            if msg:
                msg=msg+', '+'Charge imbalance'
            else:
                msg='Charge imbalance'
        return update_rxn(Rdata,Pdata,hc_prod=hc_prod,hcrct=addedhc,
                        rxnsmiles0=rxnsmiles0,msg=msg)
        


def buildrxn(Rdata,Pdata):
    '''
    Takes reaction and product data, and concatenates smiles strings, forming
    reaction smiles/smarts

    '''
    LHS=[Rdata[ID]['smiles'] for ID in Rdata for _ in range(Rdata[ID]['count'])]
    RHS=[Pdata[ID]['smiles'] for ID in Pdata for _ in range(Pdata[ID]['count'])]
    return '>>'.join([getfragments(LHS,smiles=True), getfragments(RHS,smiles=True)])


def update_stoich(stoich,compdict,hcID=None,hc_Dict=None):
    '''
    Based on balanced stoichiometry output of balance_stoichiometry function from chempy, and given a 
    dictionary of help compounds and relevant help IDs, updates species dictionary
    '''
    usedid=[]
    formdict={}
    msg=''   
    for ID in compdict:
        form=compdict[ID]['formula']
        if form not in formdict:
            formdict.update({form:[ID]})
        else:
            formdict[form].extend([ID])
#     breakpoint()
    if hcID:
        if hc_Dict is None:
            raise CustomError("Please supply help compound reference dictionary/dataframe")
        for hcid in hcID:
            form=hc_Dict[hcid]['formula']
            if form not in formdict:
                formdict.update({form:[hcid]})
            else:
                formdict[form].extend([hcid])
    for formula,coeff in stoich.items():
        if formula in formdict:
            for ID in formdict[formula]:
                if ID not in compdict:
                    compdict.update({ID:hc_Dict[ID]})
                compdict[ID]['count']=coeff
                usedid+=[ID]
        else:
            msg='Invalid balancing. Formula indicated in stoich outside compound dictionary'
            break
    if msg:
        return 'Error',msg,formdict
    else:
        valid=True
        unusedid=set(compdict.keys())-set(usedid)
        if unusedid:
            for ID in unusedid:
                if 'rxs_ids' not in compdict[ID]: #Identifying help compounds
                    valid=False
                    break
            if valid:
                for ID in unusedid:
                    del compdict[ID] # 'Reactants' that aren't reactants..these are removed
        if valid and compdict:
            return compdict,msg,formdict
        else:
#             breakpoint()
            msg='Invalid balancing. Species missing:'+','.join([str(unused) for unused in unusedid])
            return compdict,msg,formdict
        

def update_rxn(Rdata,Pdata,reac=None,prod=None,hc_prod=None,hcprod=[],hcrct=[],rxnsmiles0=None,msg=None):
    '''
    Wrapper for calling update_stoich function 
    
    '''
    stoichupdated=False
    addmsgr=''
    addmsgp=''
    if reac is not None:#Switched order
        stoichupdated=True
        Rdata1,addmsgr,formdictr=update_stoich(reac,Rdata)
    if prod is not None:
        stoichupdated=True
        Pdata1,addmsgp,formdictp=update_stoich(prod,Pdata,hcID=hcprod,hc_Dict=hc_prod)     
    if addmsgr and addmsgr!='Valid':
        if msg is not None:
            msg=addmsgr+' from LHS'+', '+msg
        else:
            msg=addmsgr+' from LHS'
    if addmsgp and addmsgp!='Valid':
        if msg is not None:
            msg=addmsgp+' from RHS'+', '+msg
        else:
            msg=addmsgp+' from RHS'
    if stoichupdated:
        if Rdata1!='Error' and Pdata1!='Error':
            Rdata=Rdata1
            Pdata=Pdata1
            try:
                balrxnsmiles=buildrxn(Rdata,Pdata)
            except Exception:
                balrxnsmiles='Error'
                msg='Invalid balancing. Species missing in database'+', '+msg #Just to make sure, although this error should never happen
        else:
            balrxnsmiles='Error'
    else:
        balrxnsmiles=buildrxn(Rdata,Pdata)
                
    if hcrct:
        LHSids=[ID for ID in Rdata if ID not in hcrct for _ in range(Rdata[ID]['count'])]
    else:
        LHSids=[ID for ID in Rdata for _ in range(Rdata[ID]['count'])]
    if hcprod:
        RHSids=[ID for ID in Pdata if ID not in hcprod for _ in range(Pdata[ID]['count'])]
    else:
        RHSids=[ID for ID in Pdata for _ in range(int(Pdata[ID]['count']))]
    if rxnsmiles0 is not None:
        return rxnsmiles0,balrxnsmiles,msg,LHSids,RHSids,hcrct,hcprod,Rdata,Pdata  #Same number and type of atoms reactant and product side and same charge ie. perfectly balanced reaction. Pretty much impossible.
    else:
        return balrxnsmiles,msg,LHSids,RHSids,hcrct,hcprod,Rdata,Pdata
    
                
def tryhelp(hc_atomtype,chempyr,chempyp,coefflim=6):
    '''
    Attempts to balance reaction with addition of help compounds in helpCompounds.py.
    
    hc_atomtype is a dictionary or list of help compounds
    chempyr is a set of LHS species formulae
    chemyp is a set of RHS species formulae
    coefflim is the maximum acceptable stoichiometric coefficient value after balancing
    
    Returns reac,prod (outputs of balance_stoichiometry function in chempy) and a list of help compounds added

    '''
#     breakpoint()
    reac={}
    prod={}
    hcid=None
    lim=len(hc_atomtype)
    keylist=[hcid for hcid in hc_atomtype]
    counter=0
    invalid=True
    while any([idx>coefflim for tup in zip(reac.values(),prod.values()) for idx in tup]) or not reac or invalid:
        if counter>lim-1:
            print('Help compounds did not help. Extra reactant atoms')
            raise CustomError("Help compounds unable to balance") 
        if hcid is not None:
            chempyp.remove(hc_atomtype[hcid]['formula'])
        hcid=keylist[counter]
        chempyp.add(hc_atomtype[hcid]['formula'])
        try:
            reac, prod = balance_stoichiometry(chempyr, chempyp,underdetermined=None,allow_duplicates=True)
            if any(idx<0 for idx in reac.values()) or any(idx<0 for idx in prod.values()): #Don't want negative stoich coefficients
                invalid=True
                raise Exception
            else:
                invalid=False
                counter+=1
                continue
        except Exception:
            counter+=1
            continue
    print('Reaction successfully balanced')
    return reac,prod,[hcid]    

def balance(Rdata,Pdata,hc_prod={},balbefore=True,coefflim=6,addedspecies=[],addedhc=[],hc_react={}):
    '''
    Balances reaction given LHS and RHS species by invoking the balance_stoichiometry function from ChemPy
    
    '''
#     breakpoint()
    chempyr=[Rdata[ID]['formula'] for ID in Rdata]
    chempyp=[Pdata[ID]['formula'] for ID in Pdata]
#     chempyr=[Rdata[ID]['formula'] for ID in Rdata for _ in range(Rdata[ID]['count'])]
#     chempyp=[Pdata[ID]['formula'] for ID in Pdata for _ in range(Pdata[ID]['count'])]
    if len(set(chempyr+chempyp))<len(chempyr+chempyp): #Isomers/same formulae present
        raise CustomError('Isomers detected. Balancing will not work')
    chempyr=set(chempyr)
    chempyp=set(chempyp)
#     highcoeff=False
    reac0={}
    prod0={}
    msg='Balanced'
    addedstr=''
    if addedspecies:
        addedstr=' with species: '+(','.join([str(species) for species in addedspecies]))
    if addedhc:
        addedstr2=' with help reactant(s): '+(','.join([hc_react[species]['formula'] for species in addedhc]))
        if addedstr:
            addedstr+', '+addedstr2
        else:
            addedstr=addedstr2
    msg+=addedstr
    
    if balbefore or not hc_prod: #Either user has indicated that first pass at balancing needs to be done or fails to specify help compounds
        try:
            reac, prod = balance_stoichiometry(chempyr, chempyp,underdetermined=None,allow_duplicates=True) #Try balancing once without adding compounds
            if any(idx<0 for idx in reac.values()) or any(idx<0 for idx in prod.values()): #Don't want negative stoich coefficients
                raise Exception
            elif any([idx>coefflim for tup in zip(reac.values(),prod.values()) for idx in tup]):
                raise Exception
 
        except Exception as e:
            if not hc_prod: #User has not supplied help compounds
                raise CustomError('Reaction does not balance. Help compounds not provided')
            pass
        else: #Reaction successfully balanced
            return reac,prod,None,msg
            
    if hc_prod: #Can try help compounds
        try:
            reac,prod,hcid=tryhelp(hc_prod,chempyr,chempyp,coefflim=coefflim)
        except Exception as e:
            raise CustomError('Reaction does not balance even with help compounds')
        else:
            hclist=','.join([hc_prod[hc]['formula'] for hc in hcid])
            return reac,prod,hcid,msg+' with help product(s): '+hclist
    

def findmatch(atomdeficit,atomdict):
    '''
    Calculates proportion of atoms mapped, based on given atom deficit dictionary and
    atom dictionary of a candidate
    '''
#     breakpoint()
    rem2=Counter()
    rem2.update(atomdict)
    rem2.subtract(Counter(atomdeficit))
    if any([val<0 for val in rem2.values()]):
        if rem2.keys()==atomdict.keys() and len(set(Counter({k:rem2[k]/atomdict[k] for k in rem2}).values()))==1:
            return 1.0 #Exact multiple
        else:
            return False
    mapprop=1-(sum(rem2.values())/sum(atomdict.values()))
    return round(mapprop,1)

def resolvecandidates(postype,Rdata,specdict,candidates,update=True,validate=True):
#     breakpoint()
    if validate:
        combinedkeys={key for candi in candidates for key in specdict[candi]['atomdict'].keys()}
        if not set(postype.keys()).issubset(combinedkeys): #Reactants/reagents/solvents cannot account for atom imbalance
            msg='LHS species insufficient'
            return Rdata,candidates,msg
            
    if len(candidates)>1: #Multiple candidates
        matches=[findmatch(postype,specdict[candi]['atomdict']) for candi in candidates]
        if 1 in matches:
            candidates=[candidates[matches.index(1)]]
        elif all(match is not False for match in matches):
            if len(postype)==1 and 'H' in postype: #Deficit is only hydrogen, so mapper will not help
                counter=Counter(matches)
                maxmatch=max(matches)
                if counter[maxmatch]==1:
                    candidates=[candidates[matches.index(maxmatch)]]
            elif len(set(matches))==len(matches):
                index_max = max(range(len(matches)), key=matches.__getitem__) #Species with maximum atoms mappable
                candidates=[candidates[index_max]]
            
    if len(postype)==1 and 'H' in postype and len(candidates)>1: #Still multiple options
        msg='Reducing agent/hydrogen carriers required and available: '+','.join([str(candi) for candi in candidates])
        return Rdata,candidates,msg      
    msg='Valid'
    if update:
        for candi in candidates:
            if candi in Rdata.keys():
                Rdata[candi]['count']+=1
            else:
                Rdata.update({candi:specdict[candi]})
    return Rdata,candidates,msg


def balance_rxn(rxnsmiles0,hc_prod={},hc_react={},coefflim=6,usemapper=True):
#     breakpoint()
    try:
        splitrxn=rxnsmiles0.split('>>')
        if len(splitrxn)==1: #Only reactants specified
            raise Exception
        rcts=splitrxn[0].split('.')
        prods=splitrxn[1].split('.')
        rcts=[Chem.MolToSmiles(molfromsmiles(rct)) for rct in rcts]
        prods=[Chem.MolToSmiles(molfromsmiles(prod)) for prod in prods]
    except Exception:
        print('Please supply valid reaction smiles. Reactant.Reactant >> Product.Product')
    rcts=Counter(rcts)
    prods=Counter(prods)
    Rdata={}
    Pdata={}
    for i,rct in enumerate(rcts):
        Rdata.update(getcompdict(ID=i,smiles=rct))
        Rdata[i]['count']=rcts[rct]
    for j,prod in enumerate(prods):
        Pdata.update(getcompdict(ID=j,smiles=prod))
        Pdata[j]['count']=prods[prod]
    addedspecies=[i for i in Rdata]
    return balancerxn(Rdata,Pdata,rxnsmiles0=rxnsmiles0,addedspecies=addedspecies,usemapper=usemapper,hc_prod=hc_prod,
                      first=False,hc_react=hc_react,coefflim=coefflim)

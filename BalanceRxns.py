# %load ./BalanceRxns.py
import rdkit
from MainFunctions import CustomError,getfragments, getcompdict, molfromsmiles
import multiprocessing
import time
# from func_timeout import func_timeout, FunctionTimedOut
from chempy import balance_stoichiometry
import copy
from collections import Counter
from decimal import Decimal, ROUND_HALF_UP
from rdkit import Chem #Importing RDKit
from rdkit.Chem import rdChemReactions #Reaction processing
# from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from MapRxns import maprxn
from math import ceil


def balance_analogue_(analoguerxns,refbalrxns=None,coefflim=6,reaxys_update=True,includesolv=True,
                     usemapper=True,helpprod=True,helpreact=False,addrctonly=False,ignoreH=False,ncpus=16,restart=True):
    '''
    Applies balance_analogue across each row of a given dataframe
    
    '''
#     breakpoint()
    if not analoguerxns.index.name and not analoguerxns.index.names:
        idxreset=True
    else:
        idxreset=False
    idxcol=[]
    if reaxys_update:
        idxcol=['ReactionID','Instance']
    else:
        idxcol=['ReactionID']
    if refbalrxns is not None:
        analoguerxns,commondf=userefrxns(analoguerxns,idxcol=idxcol,refanaloguerxns=refbalrxns)
        idxreset=False
    if not analoguerxns.empty:
        if ncpus>1:
            if restart:
                initray(num_cpus=ncpus)
            if not idxreset:
                analoguerxns.reset_index(inplace=True)
                idxreset=True
            analoguerxnsdis=mpd.DataFrame(analoguerxns)
        else:
            analoguerxnsdis=analoguerxns
        balrxns=analoguerxnsdis.apply(balance_analogue,coefflim=6,includesolv=includesolv,usemapper=usemapper,
                                   helpprod=helpprod,helpreact=helpreact,addrctonly=addrctonly,ignoreH=ignoreH,axis=1,result_type='reduce')
        balrxns=pd.Series(data=balrxns.values,index=balrxns.index) #Optional convert modin back to pandas
        analoguerxnsbal=pd.DataFrame(balrxns,columns=['rxnsmiles'])
        analoguerxnsbal[['rxnsmiles0', 'balrxnsmiles','msg','LHS','RHS','hcrct','hcprod','LHSdata','RHSdata']] = pd.DataFrame(analoguerxnsbal['rxnsmiles'].tolist(), index=analoguerxnsbal.index)
        analoguerxnsbal[['NumRefs','NumSteps','NumStages']]=analoguerxns[['NumRefs','NumSteps','NumStages']]
        cols=['NumRefs','NumSteps','NumStages','rxnsmiles0','balrxnsmiles','msg','LHS','RHS','hcrct','hcprod','LHSdata','RHSdata']
        if idxreset:
            analoguerxnsbal[idxcol]=analoguerxns[idxcol]
            cols=idxcol+cols
        analoguerxnsbal=analoguerxnsbal[cols]
        if idxreset:
            analoguerxnsbal.set_index(idxcol,inplace=True)
        if refbalrxns and not commondf.empty: #Indices need to match!
            analoguerxnsbal=pd.concat([analoguerxnsbal,commondf])
    else:
        analoguerxnsbal=commondf
        balrxns=[]
    return balrxns,analoguerxnsbal

def balance_analogue(row,basic=True,balance=True,coefflim=6,includesolv=True,
                     usemapper=True,helpprod=True,helpreact=False,addrctonly=False,ignoreH=False): #More reliable
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
        return balancerxn(Rdata,Pdata,Rgtdata=Rgtdata,Solvdata=Solvdata,rxnsmiles0=rxnsmiles0,usemapper=usemapper,
                          coefflim=coefflim,hc_prod=hc_prod,hc_react=hc_react,addrctonly=addrctonly,ignoreH=ignoreH)

def balancerxn(Rdata,Pdata,Rgtdata={},Solvdata={},rxnsmiles0=None,first=True,usemapper=True,
               addedspecies=[],addedhc=[],hc_prod={},hc_react={},coefflim=6,msg='',mandrcts={},addrctonly=False,
               ignoreH=False):
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
    mandrcts is a list of mandatory reactants (To avoid balancer from removing them)
    ignoreH refers to a boolean switch (True if all hydrogen species except H2 are ignored/not added)
    
    '''
    
    if Solvdata: #Add sovents
        Rgtdata={**Rgtdata,**Solvdata} #Newly added
        
    #%% Initialize added species/addedhc (help compound) and add small species
    if first:
        addedspecies=[]
        addedhc=[]
        if not mandrcts:
            mandrcts=copy.deepcopy(Rdata)
        if Rgtdata:
#             smallspecies=[spec for spec in Rgtdata if Rgtdata[spec]['formula'] in ['H2','O2','H2O','OH2','HOH'] if spec not in Solvdata]
#             smallspecies=[spec for spec in Rgtdata if Rgtdata[spec]['formula'] in ['H2','O2','H2O','OH2','HOH']]
            smallspecies=[spec for spec in Rgtdata if Rgtdata[spec]['formula'] in ['H2','O2']]
            if smallspecies:
                Rdata.update({spec:Rgtdata[spec] for spec in smallspecies})
                addedspecies+=smallspecies


    #%% String output handling
    addedstr=''
    if addedspecies:
        addedstr=','.join([str(species) for species in set(addedspecies) if species not in mandrcts])
        if addedstr:
            addedstr=' with species: '+addedstr
    if addedhc:
        addedstr2=','.join([hc_react[species]['formula'] for species in set(addedhc) if species not in mandrcts])
        if addedstr2:
            addedstr2=' with help reactant(s): '+addedstr2
        if addedstr:
            addedstr=addedstr+', '+addedstr2
        else:
            addedstr=addedstr2
#     if 'Mandatory' in msg or 'Smiles discrepancy' in msg:
    if 'Smiles discrepancy' in msg:
        msg=msg+addedstr
        return update_rxn(mandrcts,Pdata,hc_prod=hc_prod,hcrct=addedhc,rxnsmiles0=rxnsmiles0,msg=msg)
    if 'Hydrogen carriers' in msg:
        msg=msg+addedstr
        return update_rxn(Rdata,Pdata,hc_prod=hc_prod,hcrct=addedhc,rxnsmiles0=rxnsmiles0,msg=msg)            

    Rcount=sum([Counter(Rdata[ID]['atomdict']) for ID in Rdata for _ in range(Rdata[ID]['count'])],start=Counter()) #Sum of atom counts/types on LHS
    Rcharge=sum([Rdata[ID]['charge'] for ID in Rdata for _ in range(Rdata[ID]['count'])])
    Pcount=sum([Counter(Pdata[ID]['atomdict']) for ID in Pdata for _ in range(Pdata[ID]['count'])],start=Counter()) #Sum of atom counts/types on RHS
    Pcharge=sum([Pdata[ID]['charge'] for ID in Pdata for _ in range(Pdata[ID]['count'])])
            
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

    if postype: #Reactants, Reagents may be needed 
        status=[mandrcts[ID0]['count']>Rdata[ID0]['count'] for ID0 in mandrcts if ID0 in Rdata] #if ID0 in Rdata
        if any(status):
            if msg:
                msg=msg+', '+'Mapping error'+addedstr
            else:
                msg='Mapping error'+addedstr
            return update_rxn(Rdata,Pdata,hc_prod=hc_prod,hcrct=addedhc,rxnsmiles0=rxnsmiles0,msg=msg)
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
            Rdata,addedspecies_,msg_=resolvecandidates(postype,Rdata,Rdata,candirxt,Pdata,validate=False,rctonly=addrctonly,ignoreH=ignoreH)
        elif candirgt and not candirxt:
            Rdata,addedspecies_,msg_=resolvecandidates(postype,Rdata,Rgtdata,candirgt,Pdata,validate=False,rctonly=False,ignoreH=ignoreH)
        elif candirxt and candirgt:
            Rdata,addedspecies_,msg_=resolvecandidates(postype,Rdata,{**Rdata,**Rgtdata},candirxt+candirgt,Pdata,validate=False,rctonly=addrctonly,ignoreH=ignoreH)
        elif not candirxt and not candirgt:
            combineddict={**Rdata,**Rgtdata}
            candispec=[specid for specid in combineddict if set(combineddict[specid]['atomdict'].keys()).intersection(set(postype.keys()))]
            if not candispec:
                msg_='LHS species insufficient'
            else:
                Rdata,addedspecies_,msg_=resolvecandidates(postype,Rdata,combineddict,candispec,Pdata,rctonly=addrctonly,ignoreH=ignoreH)
        elif candihc and not candirxt and not candirgt:
            Rdata,addedhc_,msg_=resolvecandidates(postype,Rdata,hc_react,candihc,Pdata,validate=False,rctonly=addrctonly,ignoreH=ignoreH)
        if 'Hydrogen carriers' in msg_: #Multiple candidates for hydrogen carriers (mapper won't help)
            if msg:
                msg=msg+', '+msg_
            else:
                msg=msg_
            
        elif msg_!='Valid': #Atom surplus on RHS cannot be met by any LHS species
            if msg:
                msg=msg+', '+msg_+addedstr
            else:
                msg=msg_+addedstr
            return update_rxn(Rdata,Pdata,hc_prod=hc_prod,hcrct=addedhc,rxnsmiles0=rxnsmiles0,msg=msg)

        else:
            addedspecies+=addedspecies_
            if addedhc_:
                adddedhc+=addedhc_
            #New#
            if len(postype)==1 and 'H' in postype and 'With hydrogen carriers' not in msg:
                if msg:
                    msg=msg+', '+'With hydrogen carriers: '+','.join([str(addedspec) for addedspec in addedspecies_])
                else:
                    msg='With hydrogen carriers: '+','.join([str(addedspec) for addedspec in addedspecies_])
            #New#
#         breakpoint()
        return balancerxn(Rdata,Pdata,Rgtdata=Rgtdata,rxnsmiles0=rxnsmiles0,first=False,usemapper=usemapper,
                         addedspecies=addedspecies,addedhc=addedhc,hc_prod=hc_prod,
                         hc_react=hc_react,coefflim=coefflim,msg=msg,mandrcts=mandrcts,addrctonly=addrctonly,ignoreH=ignoreH) 
    
    elif negtype:
#         breakpoint()
        if usemapper and len(set(addedspecies))>1: #More than one choice or added small species, let the mapper decide
            rxnsmiles=buildrxn(Rdata,Pdata)
            mapped_rxn=maprxn([rxnsmiles])[0]
            if mapped_rxn=='Error':
                if not addrctonly:
                    addrctonly=True
                    candidates2=[candi for candi in Rdata if candi in mandrcts]
                    if candidates2 and candidates2!=list(Rdata.keys()):
                        Rdata={ID0:Rdata[ID0] for ID0 in candidates2}
                        addedspecies=[addedspec for addedspec in addedspecies if addedspec in Rdata]
                        if addedhc:
                            addedhc=[addedh for addedh in addedhc if addedh in Rdata]
                        return balancerxn(Rdata,Pdata,rxnsmiles0=rxnsmiles0,first=False,usemapper=usemapper,
                                         addedspecies=addedspecies,addedhc=addedhc,hc_prod=hc_prod,hc_react=hc_react,
                                         coefflim=coefflim,msg=msg,mandrcts=mandrcts,addrctonly=addrctonly,ignoreH=ignoreH)
                
                if all([Rdata[ID0]['count']==1 for ID0 in Rdata]) or any([Rdata[ID0]['count']>=10 for ID0 in Rdata]):
                    if msg:
                        msg=msg+', '+'Mapping error'+addedstr
                    else:
                        msg='Mapping error'+addedstr
                    return update_rxn(Rdata,Pdata,hc_prod=hc_prod,hcrct=addedhc,rxnsmiles0=rxnsmiles0,msg=msg)
                else:  
                    status=[mandrcts[ID0]['count']>=Rdata[ID0]['count'] for ID0 in mandrcts if ID0 in Rdata] #if ID0 in Rdata
                    if any(status):
                        for ID0 in Rdata:
                            Rdata[ID0]['count']=Rdata[ID0]['count']-1
                    else:
                        mandrcts=copy.deepcopy(Rdata)
                        mincount=min([Rdata[ID0]['count'] for ID0 in Rdata])
                        for ID0 in Rdata:
                            Rdata[ID0]['count']=mincount
                    return balancerxn(Rdata,Pdata,rxnsmiles0=rxnsmiles0,first=False,usemapper=usemapper,
                                    addedspecies=addedspecies,addedhc=addedhc,hc_prod=hc_prod,hc_react=hc_react,
                                    coefflim=coefflim,msg=msg,mandrcts=mandrcts,addrctonly=addrctonly,ignoreH=ignoreH)
            else:
                mappedrxn=mapped_rxn.get('mapped_rxn')
                if 'With hydrogen carriers' in msg:
                    hcarriers=[int(hcarrier) for hcarrier in msg.split('With hydrogen carriers: ')[1].split(', ')[0].split(',')]
                else:
                    hcarriers=[]
#                 breakpoint()
                LHSdata,_,msg_=checkrxn(mappedrxn,Rdata=Rdata,updateall=False,mandrcts=list(mandrcts.keys()),hcarriers=hcarriers)
                if 'Mandatory' in msg_ and not addrctonly:#Reactants are unmapped (Try again with a smaller candidate selection)
                    return balancerxn(mandrcts,Pdata,Rgtdata=Rgtdata,rxnsmiles0=rxnsmiles0,first=True,usemapper=usemapper,
                             addedspecies=addedspecies,addedhc=addedhc,hc_prod=hc_prod,
                             hc_react=hc_react,coefflim=coefflim,msg=msg,addrctonly=True,ignoreH=ignoreH)
                if msg:
                    if msg_!='Valid':
                        msg=msg_+', '+msg
                    msg='Mapper used'+', '+msg
                else:
                    if msg_!='Valid':
                        msg='Mapper used'+', '+msg_
                    else:
                        msg='Mapper used'
            if addedhc:
                addedhc=[addedh for addedh in addedhc if addedh in LHSdata]  
#             breakpoint()
            addedspecies=[addedspec for addedspec in addedspecies if addedspec in LHSdata]
            return balancerxn(LHSdata,Pdata,Rgtdata=Rgtdata,rxnsmiles0=rxnsmiles0,first=False,usemapper=False,
                             addedspecies=addedspecies,addedhc=addedhc,hc_prod=hc_prod,hc_react=hc_react,coefflim=coefflim,
                             msg=msg,mandrcts=mandrcts,addrctonly=addrctonly,ignoreH=ignoreH)
        else:
#             breakpoint()
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
                        addedhc=[addedh for addedh in addedhc if addedh in Rdata]
#                         refhc=Counter(addedhc)
#                         addedhc=[addedspec for addedspec in Rdata for _ in range(min([Rdata[addedspec]['count'],refhc[addedspec]])) if addedspec in addedhc]
                    addedspecies=[addedspec for addedspec in addedspecies if addedspec in Rdata]
#                     refspec=Counter(addedspecies)
#                     addedspecies=[addedspec for addedspec in Rdata for _ in range(min([Rdata[addedspec]['count'],refspec[addedspec]])) if addedspec in addedspecies if addedspec not in addedhc] 
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
                                  coefflim=coefflim,msg=msg,addrctonly=addrctonly,ignoreH=ignoreH)
                        
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
                                            addedspecies=[addedspec for addedspec in addedspecies if addedspec not in mandrcts],
                                            addedhc=[addedh for addedh in addedhc if addedh not in mandrcts],hc_react=hc_react)
                if msg:
                    msg=msg+', '+msg0
                else:
                    msg=msg0
            except Exception:
                if not balbefore:
                    return balancerxn(Rdata,Pdata,Rgtdata=Rgtdata,rxnsmiles0=rxnsmiles0,first=False,usemapper=False,
                                  addedspecies=addedspecies,addedhc=addedhc,hc_prod={},hc_react=hc_react,
                                  coefflim=coefflim,msg=msg,addrctonly=addrctonly,ignoreH=ignoreH)
                if Rcharge!=Pcharge:
                    if msg:
                        msg=msg+', '+'Charge imbalance'
                    else:
                        msg='Charge imbalance'
                if msg: 
                    msg=msg+', '+'RHS species insufficient'+addedstr
                else:
                    msg='RHS species insufficient'+addedstr
                
                return update_rxn(Rdata,Pdata,hcrct=addedhc,rxnsmiles0=rxnsmiles0,msg=msg)
            else:
                if Rcharge!=Pcharge:
                    if msg:
                        msg=msg+', '+'Charge imbalance'+addedstr
                    else:
                        msg='Charge imbalance'+addedstr
                return update_rxn(Rdata,Pdata,reac=reac,prod=prod,hcprod=hcid,hc_prod=hc_prod2,hcrct=addedhc,
                                  rxnsmiles0=rxnsmiles0,msg=msg)
    else:
        if Rcharge!=Pcharge:
            if msg:
                msg=msg+', '+'Charge imbalance'+addedstr
            else:
                msg='Charge imbalance'+addedstr
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
            msg='Invalid balancing. Species missing: '+','.join([str(unused) for unused in unusedid])
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
        if ('LHS species insufficient' in msg) | ('Invalid' in msg) | ('Mapping error' in msg) | ('discrepancy' in msg):
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
    

    
    
# def bal_stoichiometry(chempyr,chempyp):
#     return balance_stoichiometry(chempyr,chempyp,underdetermined=None,allow_duplicates=True)
# def bal_stoich(chempyr,chempyp):
#     try:
#         reac,prod=bal_stoichiometry(chempyr,chempyp)
#         return reac,prod
# #         reac,prod=func_timeout(5, bal_stoichiometry, args=(chempyr,chempyp))
        
#     except Exception:
#         raise Exception


    
    
                
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
            reac, prod = balance_stoichiometry(chempyr,chempyp,underdetermined=None,allow_duplicates=True)
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



def balance(Rdata,Pdata,hc_prod={},balbefore=True,coefflim=6,addedspecies=[],addedhc=[],hc_react={},mandrcts={}):
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
        addedstr=','.join([str(species) for species in set(addedspecies) if species not in mandrcts])
        if addedstr:
            addedstr=' with species: '+addedstr  
    if addedhc:
        addedstr2=' with help reactant(s): '+(','.join([hc_react[species]['formula'] for species in addedhc]))
        if addedstr:
            addedstr+', '+addedstr2
        else:
            addedstr=addedstr2
    msg+=addedstr
    
    if balbefore or not hc_prod: #Either user has indicated that first pass at balancing needs to be done or fails to specify help compounds
        try:
            reac, prod = balance_stoichiometry(chempyr,chempyp,underdetermined=None,allow_duplicates=True) #Try balancing once without adding compounds
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
#     breakpoint()   
    if hc_prod: #Can try help compounds
        try:
            reac,prod,hcid=tryhelp(hc_prod,chempyr,chempyp,coefflim=coefflim)
        except Exception as e:
            raise CustomError('Reaction does not balance even with help compounds')
        else:
            hclist=','.join([hc_prod[hc]['formula'] for hc in hcid])
            return reac,prod,hcid,msg+' with help product(s): '+hclist
    

def findmatch(atomdeficit,atomdict,strict=True,returnmultdict=True):
    '''
    Calculates proportion of atoms mapped, based on given atom deficit dictionary and
    atom dictionary of a candidate
    '''
#     breakpoint()
    if not set(atomdict.keys()).intersection(set(atomdeficit.keys())): #No match at all
        return False,False
    rem2=Counter()
    rem2.update(atomdict)
    rem2.subtract(Counter(atomdeficit))
    if any([val<0 for val in rem2.values()]):
        multdict={k:abs(atomdeficit[k]/atomdict[k]) for k in atomdeficit if k in atomdict}
        if (atomdeficit.keys()==atomdict.keys()) & (len(set(Counter(multdict).values()))==1):
            return 1.0,1 #Exact multiple
        elif returnmultdict:
            return False,multdict
        else:
            if strict:
                mult=int(ceil(max(multdict.values())))
            else:
                mult=int(ceil(min(multdict.values())))
            return False,mult
    mapprop=(sum(atomdeficit.values()))/(sum(atomdict.values()))
#     mapprop=1-(sum(rem2.values())/sum(atomdict.values()))
    return round(mapprop,1),1



def resolvecandidates(postype,Rdata,specdict,candidates,Pdata,update=True,validate=True,rctonly=False,
                      coefflim=6,ignoreH=False):
    '''
    Resolves candidates based on atom deficit (postype) and supplied candidate matches
    
    
    '''
    
#     breakpoint()
    msg=''
    if validate:
        combinedkeys={key for candi in candidates for key in specdict[candi]['atomdict'].keys()}
        if not set(postype.keys()).issubset(combinedkeys): #Reactants/reagents/solvents cannot account for atom imbalance
            msg='LHS species insufficient'
            return Rdata,candidates,msg
    if len(candidates)>1:
        matches=[]
        mult=[]
        for candi in candidates:
            match,mult_=findmatch(postype,specdict[candi]['atomdict'])
            matches+=[match]
            mult+=[mult_] 
        if 1 in matches:
            index_max=[i for i,match in enumerate(matches) if match==1]
            candidates=[candidates[idx] for idx in index_max]
            mult=[mult[idx] for idx in index_max]
        else:
            if len(candidates)>1:
                if len(postype)==1 and 'H' in postype and ignoreH:
                    msg='Hydrogen carriers: '+','.join([str(candi) for candi in candidates])
                    return Rdata,candidates,msg
                if all(match is not False for match in matches):
                    if len(postype)==1 and 'H' in postype: #Deficit is only hydrogen, so mapper will not help
                        counter=Counter(matches)
                        maxmatch=max(matches)
                        if counter[maxmatch]==1:
                            index_max=matches.index(maxmatch)
                            candidates=[candidates[index_max]]
                            mult=[mult[index_max]]
                    elif len(set(matches))==len(matches): #Possibility of mapping wrong species still there
                        index_max = max(range(len(matches)), key=matches.__getitem__) #Species with maximum atoms mappable
                        candidates=[candidates[index_max]]
                        mult=[mult[index_max]]
                else: #Higher stoichiometric coefficients needed or more than one species needed
#                     breakpoint()
                    atompop={k:[] for k in postype}
                    for candi,mult_ in zip(candidates,mult):
                        if type(mult_)==dict:
                            for k in mult_:
                                if k in atompop:
                                    atompop[k].extend([candi])
                        else:
                            for k in specdict[candi]['atomdict']:
                                if k in atompop:
                                    atompop[k].extend([candi])
#                     breakpoint()
                    for k in sorted(atompop, key=lambda k: (len(atompop[k]),postype[k])):
                        if rctonly:
                            extrarct=[rctid for rctid in Rdata if k in Rdata[rctid]['atomdict'] if rctid not in candidates]
                        else:
                            extrarct=[rctid for rctid in specdict if k in specdict[rctid]['atomdict'] if rctid not in candidates]
                        if extrarct:
                            for rctid in extrarct:
                                mult_=int(ceil(postype[k]/specdict[rctid]['atomdict'][k]))
                                if mult_<=1:
                                    mult_=1
                                    matches+=[True]
                                else:
                                    matches+=[False]
                                mult+=[mult_]
                            candidates+=extrarct
                            atompop[k].extend(extrarct)
                        if len(atompop[k])==1: #Only one candidate for atom
                            candi=atompop[k][0]
                            mult_=mult[candidates.index(candi)]
                            if type(mult_)==dict:
                                mult=[int(ceil(mult_[k]))]
                            else:
                                mult=[mult_]
                            candidates=[candi]
                            break
                        else:
                            if len(postype)==1 and 'H' in postype:
                                msg='Hydrogen carriers: '+','.join([str(candi) for candi in candidates])
                                return Rdata,candidates,msg
                            matches2=[]
                            candidates2=[]
                            for candi in atompop[k]:
                                idx=candidates.index(candi)
                                mult_=mult[idx]
                                if type(mult_)==dict:
                                    mult_=int(ceil(mult_[k]))
                                mult[idx]=mult_
                                match=matches[idx]
                                if match==False:
                                    atomdict=copy.deepcopy(specdict[candi]['atomdict'])
                                    totatomdict={j:atomdict[j]*mult_ for j in atomdict}
                                    totmapped={j:min(totatomdict[j],postype[j]) for j in postype if j in atomdict}
                                    mapprop=round((sum(totmapped.values()))/(sum(totatomdict.values())),1)
                                    matches[idx]=mapprop
                                    matches2+=[mapprop]
                                else:
                                    matches2+=[match]
                            counter=Counter(matches2)
                            maxmatch=max(matches2)
                            if rctonly:
                                candidates2=[candi for candi in atompop[k] if candi in Rdata]
                            if not candidates2:
                                candidates2=[candi for candi in atompop[k] if matches[candidates.index(candi)]==maxmatch or candi in Rdata]
                            mult2=[mult[candidates.index(candi)] for candi in candidates2]
                            candidates=candidates2
                            mult=mult2
                            break
    else:            
        mult=[1]
                                
    if len(candidates)>1  and 'Hydrogen carriers' not in msg: #or ignoreH
        if len(postype)==1 and 'H' in postype: #Still multiple options
            msg='Hydrogen carriers: '+','.join([str(candi) for candi in candidates])
            return Rdata,candidates,msg
#     breakpoint()
    msg='Valid'
    if update:
        for candi,mult_ in zip(candidates,mult):
            if candi in Rdata.keys():
                Rdata[candi]['count']+=mult_
            else:
                Rdata.update({candi:specdict[candi]})
                Rdata[candi]['count']=mult_
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

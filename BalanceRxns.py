# %load ./BalanceRxns.py

from MainFunctions import CustomError,getfragments
from chempy import balance_stoichiometry
import copy
from collections import Counter

def buildrxn(Rdata,Pdata):
    '''
    Takes reaction and product data, and concatenates smiles strings, forming
    reaction smiles/smarts

    Parameters
    ----------
    Rdata : dict
        Reaction data, formatted as output of getcompdict()
    Pdata : dict
        Product data, formatted as output of getcompdict()

    Returns
    -------
    Output: str
        Reaction smiles/smarts string

    '''
    LHS=[Rdata[ID]['smiles'] for ID in Rdata for _ in range(Rdata[ID]['count'])]
    RHS=[Pdata[ID]['smiles'] for ID in Pdata for _ in range(Pdata[ID]['count'])]
    return '>>'.join([getfragments(LHS,smiles=True), getfragments(RHS,smiles=True)])


def update_stoich(stoich,compdict,hcID=None,hc_Dict=None):
    '''
    

    Parameters
    ----------
    stoich : TYPE
        DESCRIPTION.
    compdict : TYPE
        DESCRIPTION.
    hcID : TYPE, optional
        DESCRIPTION. The default is None.
    hc_Dict : TYPE, optional
        DESCRIPTION. The default is None.

    Raises
    ------
    CustomError
        DESCRIPTION.

    Returns
    -------
    compdict : TYPE
        DESCRIPTION.
    msg : TYPE
        DESCRIPTION.

    '''
    usedid=[]
    for formula,coeff in stoich.items():
        assigned=False
        for ID in compdict:
            if compdict[ID]['formula']==formula:
                compdict[ID]['count']=coeff
                assigned=True
                usedid+=[ID] #added
                break
        if assigned:
            continue
        elif hcID:
            if hc_Dict is None:
                raise CustomError("Please supply help compound reference dictionary/dataframe")
            for hcid in hcID:
                if formula==hc_Dict[hcid]['formula']: 
                    if hcid not in compdict.keys():
                        compdict.update({hcid:hc_Dict[hcid]})
                        compdict[hcid]['count']=coeff
                        assigned=True
                        usedid+=[hcid] #added
                        break
                    else: #added
                        compdict[hcid]['count']=coeff
                        assigned=True
                        usedid+=[hcid]
                        break
        if not assigned:
            raise CustomError("Formula indicated in stoich outside compound dictionary")
    valid=True
    unusedid=set(compdict.keys())-set(usedid) #added
    if unusedid:
        for ID in unusedid:
            if 'rxs_ids' not in compdict[ID]:
                valid=False
                break
        if valid:
            for ID in unusedid:
                del compdict[ID]
    if valid:
        msg='Valid'
    else:
        msg='Invalid balancing. Species missing'
    return compdict,msg

def update_rxn(Rdata,Pdata,reac=None,prod=None,hc_prod=None,hcprod=[],hcrct=[],rxnsmiles0=None,msg=None):
    if reac is not None:#Switched order
        Rdata,addmsgr=update_stoich(reac,Rdata)
        if addmsgr!='Valid':
            if msg is not None:
                msg=addmsgr+' from LHS'+', '+msg
            else:
                msg=addmsgr+' from LHS'
    if prod is not None:
        Pdata,addmsgp=update_stoich(prod,Pdata,hcID=hcprod,hc_Dict=hc_prod)
        if addmsgp!='Valid':
            if msg is not None:
                msg=addmsgp+' from RHS'+', '+msg
            else:
                msg=addmsgp+' from RHS'
    if hcrct:
        LHSids=[ID for ID in Rdata if ID not in hcrct for _ in range(Rdata[ID]['count'])]
    else:
        LHSids=[ID for ID in Rdata for _ in range(Rdata[ID]['count'])]
    if hcprod:
        RHSids=[ID for ID in Pdata if ID not in hcprod for _ in range(Pdata[ID]['count'])]
    else:
        RHSids=[ID for ID in Pdata for _ in range(int(Pdata[ID]['count']))]
    balrxnsmiles=buildrxn(Rdata,Pdata)  
    if rxnsmiles0 is not None:
        return rxnsmiles0,balrxnsmiles,msg,LHSids,RHSids,hcrct,hcprod,Rdata,Pdata  #Same number and type of atoms reactant and product side and same charge ie. perfectly balanced reaction. Pretty much impossible.
    else:
        return balrxnsmiles,msg,LHSids,RHSids,hcrct,hcprod,Rdata,Pdata
    
                
def tryhelp(hc_atomtype,chempyr,chempyp,coefflim=5):
    '''
    Attempts to balance reaction with addition of help compounds in helpCompounds.py

    Parameters
    ----------
    hc_atomtype : Dict/List
        Dictionary/list of help compounds
    chempyr : Set
        Set of LHS species formulae
    chempyp : Set
        Set of RHS species formulae

    Returns
    -------
    TRUE, reac, prod if reaction can be balanced and stoichiometric coefficient not 1. 
    Reac and prod are dictionaries of reactants and products with respective stoichiometric coefficients

    False, errormsg if reaction cannot be balanced. Errormsg is the error message or reason for failure

    '''
    reac={}
    prod={}
    hcid=None
    lim=len(hc_atomtype)
    keylist=[hcid for hcid in hc_atomtype]
    counter=0
    while any([idx>=coefflim for tup in zip(reac.values(),prod.values()) for idx in tup]) or not reac:
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
                raise Exception
            else:
                counter+=1
                continue
        except Exception:
            counter+=1
            continue
    print('Reaction successfully balanced')
    return reac,prod,[hcid]    

def balance(Rdata,Pdata,hc_prod=None,balbefore=True,coefflim=5,addedspecies=None,addedhc=None,hc_react=None):
    chempyr={Rdata[ID]['formula'] for ID in Rdata for _ in range(Rdata[ID]['count'])}
    chempyp={Pdata[ID]['formula'] for ID in Pdata for _ in range(Pdata[ID]['count'])} 
    highcoeff=False
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
    
    if balbefore or hc_prod is None: #Either user has indicated that first pass at balancing needs to be done or fails to specify help compounds
        try:
            reac, prod = balance_stoichiometry(chempyr, chempyp,underdetermined=None,allow_duplicates=True) #Try balancing once without adding compounds
            if any([idx>=coefflim for tup in zip(reac.values(),prod.values()) for idx in tup]):
                reac0=copy.copy(reac)
                prod0=copy.copy(prod)
                highcoeff=True
                raise Exception
            elif any(idx<0 for idx in reac.values()) or any(idx<0 for idx in prod.values()): #Don't want negative stoich coefficients
                raise Exception
        except Exception as e:
            if hc_prod is None: #User has not supplied help compounds
                if highcoeff:
                    return reac0,prod0,None,'Warning. Coeffs high'
                else:
                    raise CustomError('Reaction does not balance. Help compounds not provided')
            pass
        else: #Reaction successfully balanced
            return reac,prod,None,msg
            
    if hc_prod is not None: #Can try help compounds
        try:
            reac,prod,hcid=tryhelp(hc_prod,chempyr,chempyp,coefflim=coefflim)
        except Exception as e:
            if highcoeff:
                return reac0,prod0,None,'Warning. Coeffs high'
            else:
                raise CustomError('Reaction does not balance even with help compounds')
        else:
            hclist=','.join([hc_prod[hc]['formula'] for hc in hcid])
            return reac,prod,hcid,msg+' with help product(s): '+hclist
    

def findmatch(postype,atomdict):
    '''
    Calculates proportion of atoms mapped
    
    '''
    rem2=Counter()
    rem2.update(atomdict)
    rem2.subtract(Counter(postype))
    if any([val<0 for val in rem2.values()]):
        return False
    mapprop=1-(sum(rem2.values())/sum(atomdict.values()))
    return round(mapprop,1)

def balancerxn(Rdata,Pdata,Rgtdata={},rxnsmiles0=None,numrefs=None,first=True,
               addedspecies=[],addedhc=[],hc_prod=None,hc_react=None,coefflim=5):
    
    #%% Initialize added species/addedhc
    if first:
        addedspecies=[]
        addedhc=[]
    #%% Settle diatomic species (H2,O2 and H2O if present add them as they are probably involved)
    if first and Rgtdata:
        ds=[did for did in Rgtdata if Rgtdata[did]['smiles'] in [hc_react[hcid]['smiles'] for hcid in [1,2]]]
#         ds=[did for did in Rgtdata if Rgtdata[did]['smiles'] in [hc_react[hcid]['smiles'] for hcid in [0,1,2]]]
        if ds:
            Rdata.update({did:Rgtdata[did] for did in ds})
            addedspecies+=ds
        
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
            addedstr+', '+addedstr2
        else:
            addedstr=addedstr2
            
    #%% If reaction is balanced already
    if Rcount==Pcount and Rcharge==Pcharge:
        print('Reaction is fully balanced')
        if first:
            if rxnsmiles0 is not None:
                balrxnsmiles=rxnsmiles0
            else:
                balrxnsmiles=buildrxn(Rdata,Pdata)
            msg='Already balanced'
        else:
            balrxnsmiles=buildrxn(Rdata,Pdata)
            msg='Balanced'
        if addedstr:
            msg+=addedstr
        return update_rxn(Rdata,Pdata,hcrct=addedhc,rxnsmiles0=rxnsmiles0,msg=msg)
        
    
#%% Testing for stoichiometric multipliers (If compounds have same molecular formulae but different structural formulae)
    multrhs=True
    multlhs=False
    balcharge=True
    if Rcount.keys()==Pcount.keys():
        try:
            multset={Rcount[key]/Pcount[key] for key in Rcount} #Is there a common stoichiometric multiplier
            if len(multset)>1:
                raise Exception
        except Exception as e:
            pass
        else:
            mult=list(multset)[0]
            if mult<1:
                mult=1/mult
                multlhs=True
                multrhs=False
                if Pcharge!=mult*Rcharge:
                    balcharge=False
            elif Rcharge!=mult*Pcharge:
                balcharge=False
            if mult.is_integer() and balcharge:
                if multrhs:
                    for ID in Pdata:
                        Pdata[ID]['count']=int(mult)
                else:
                    for ID in Rdata:
                        Rdata[ID]['count']=int(mult)
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
    
    if postype: #Reactants, Reagents may be needed (postype or (postype and negtype))
        if hc_prod is not None:
            if Rcharge==Pcharge:
                hc_prod={hcid:hc_prod[hcid] for hcid in hc_prod if hc_prod[hcid]['charge']==0} #Atom type for help compounds
        #%% Initializing variables
        candirxt=[]
        candirgt=[]
        candihc=[]
        matches=[]
        #%% Get reactant match first
        candirxt=[rctid for rctid in Rdata if set(postype.keys()).issubset(set(Rdata[rctid]['atomdict'].keys()))]
        #%% Get reagent match
        candirgt=[rgtid for rgtid in Rgtdata if set(postype.keys()).issubset(set(Rgtdata[rgtid]['atomdict'].keys()))]
        if candirgt and numrefs is not None and numrefs>1: #Too many references can lead to many reagents being added at the same time (limitation of our method)
            candihc=[ID for ID in candirgt if Rgtdata[ID]['smiles'] in [hc_react[hcid]['smiles'] for hcid in [0]]] #H2O
#             candihc=[ID for ID in candirgt if Rgtdata[ID]['smiles'] in [hc_react[hcid]['smiles'] for hcid in [0,1,2]]] #H2, O2, H2O
            if candihc:
                candirgt=list(set(candihc+[candirgt[0]]))
            else:
                candirgt=[candirgt[0]]
        #%% Get help compound match
        if hc_react is not None:
            candihc=[hcid for hcid in hc_react if set(postype.keys()).issubset(set(hc_react[hcid]['atomdict'].keys()))]
            
#         breakpoint() 
        
        if candirxt and not candirgt: #Only reactant matches
            if len(candirxt)>1:
                matches=[findmatch(postype,Rdata[candi]['atomdict']) for candi in candirxt]
                if all([match for match in matches]) and len(set(matches))==len(matches): #Means deficit can be made up by all matches and no duplicates  
                    index_max = max(range(len(matches)), key=matches.__getitem__) #Species with maximum atoms mappable
                    candirxt=[candirxt[index_max]]
            try:
                if set(candirxt).issubset(set(addedspecies)): #If reactant has been added before, don't bother balancing
                    raise Exception
                reac,prod,hcid,msg=balance(Rdata,Pdata,hc_prod=hc_prod,coefflim=coefflim,
                                           addedspecies=addedspecies,addedhc=addedhc,hc_react=hc_react)
            except Exception:
                first=False
                for candi in candirxt:
                    Rdata[candi]['count']+=1
                return balancerxn(Rdata,Pdata,Rgtdata=Rgtdata,rxnsmiles0=rxnsmiles0,numrefs=numrefs,
                        first=first,addedspecies=addedspecies+candirxt,addedhc=addedhc,hc_prod=hc_prod,hc_react=hc_react) #Recursion
            else:
                return update_rxn(Rdata,Pdata,reac=reac,prod=prod,hcprod=hcid,hc_prod=hc_prod,rxnsmiles0=rxnsmiles0,msg=msg)
            
        elif candirgt and not candirxt:
            if len(candirgt)>1:
                matches=[findmatch(postype,Rgtdata[candi]['atomdict']) for candi in candirgt]
                if all([match for match in matches]) and len(set(matches))==len(matches): #Means deficit can be made up by all matches and no duplicates  
                    index_max = max(range(len(matches)), key=matches.__getitem__) #Species with maximum atoms mappable
                    candirgt=[candirgt[index_max]]
                
            for candi in candirgt:
                if candi in Rdata.keys():
                    Rdata[candi]['count']+=1
                else:
                    Rdata.update({candi:Rgtdata[candi]})
                AddedSpecies=candirgt
            try:
                if set(AddedSpecies).issubset(set(addedspecies)):
                    raise Exception
                reac,prod,hcid,msg=balance(Rdata,Pdata,hc_prod=hc_prod,coefflim=coefflim,
                                          addedspecies=addedspecies+AddedSpecies,addedhc=addedhc,hc_react=hc_react)
            except Exception:
                first=False
                return balancerxn(Rdata,Pdata,Rgtdata=Rgtdata,rxnsmiles0=rxnsmiles0,numrefs=numrefs,
                        first=first,addedspecies=addedspecies+AddedSpecies,addedhc=addedhc,hc_prod=hc_prod,
                        hc_react=hc_react) #Recursion
            else:
                return update_rxn(Rdata,Pdata,reac=reac,prod=prod,hcprod=hcid,hc_prod=hc_prod,rxnsmiles0=rxnsmiles0,msg=msg)
        
        elif candirxt and candirgt: #Choice presents issues..mapper may be incorrect in selecting proper reactant
            finalrt={}
            finmatch={}
            finalrt.update({candir:Rdata[candir] for candir in set(candirxt)})
            finalrt.update({candirg:Rgtdata[candirg] for candirg in set(candirgt)})
            matches=[findmatch(postype,finalrt[candi]['atomdict']) for candi in finalrt]
            if all([match for match in matches]) and len(set(matches))==len(matches): #Means deficit can be made up by all matches and no duplicates  
                index_max = max(range(len(matches)), key=matches.__getitem__) #Species with maximum atoms mappable
                finmatch=[list(finalrt.keys())[index_max]]
            else:
                finmatch=list(finalrt.keys())
            for candi in finmatch:
                if candi in Rdata.keys():
                    Rdata[candi]['count']+=1
                else:
                    Rdata.update({candi:Rgtdata[candi]})
            try:
                if set(finmatch).issubset(set(addedspecies)):
                    raise Exception
                reac,prod,hcid,msg=balance(Rdata,Pdata,hc_prod=hc_prod,coefflim=coefflim,
                                          addedspecies=addedspecies+list(finmatch),addedhc=addedhc,hc_react=hc_react)
            except Exception:
                first=False
                return balancerxn(Rdata,Pdata,Rgtdata=Rgtdata,rxnsmiles0=rxnsmiles0,numrefs=numrefs,
                                first=first,addedspecies=addedspecies+list(finmatch),addedhc=addedhc,hc_prod=hc_prod,
                                hc_react=hc_react) #Recursion
            else:
                return update_rxn(Rdata,Pdata,reac=reac,prod=prod,hcprod=hcid,hc_prod=hc_prod,rxnsmiles0=rxnsmiles0,msg=msg)
                
        elif candihc and not candirxt and not candirgt:
            for candi in candihc:
                if candi in Rdata.keys():
                    Rdata[candi]['count']+=1
                else:
                    Rdata.update({candi:hc_react[candi]})
            try:
                if set(candihc).issubset(set(addedhc)):
                    raise Exception
                reac,prod,hcid,msg=balance(Rdata,Pdata,hc_prod=hc_prod,coefflim=coefflim,
                                           addedspecies=addedspecies,addedhc=addedhc+candihc,hc_react=hc_react)
            except Exception:
                first=False
                return balancerxn(Rdata,Pdata,Rgtdata=Rgtdata,rxnsmiles0=rxnsmiles0,numrefs=numrefs,
                        first=first,addedspecies=addedspecies,addedhc=addedhc+candihc,hc_prod=hc_prod,
                        hc_react=hc_react) #Recursion
            else:
                return update_rxn(Rdata,Pdata,reac=reac,prod=prod,hcprod=hcid,hc_prod=hc_prod,hcrct=addedhc+candihc,
                                  rxnsmiles0=rxnsmiles0,msg=msg)
        
        else: #No matches
            if rxnsmiles0 is not None:
                return rxnsmiles0,'Error','Reactants/reagents/help compounds insufficient',list(Rdata.keys()),list(Pdata.keys()),[],[],Rdata,Pdata
            else:
                return 'Error','Reactants/reagents/help reactants insufficient',list(Rdata.keys()),list(Pdata.keys()),[],[],Rdata,Pdata
    elif negtype and not postype:
        hc_list=[]
        hc_prod2={}
#         breakpoint()
        if first:
            try:
                hc_list=[hc for hc in hc_prod if hc_prod[hc]['atomdict']==negtype]
                if not hc_list:
                    hc_list=[hc for hc in hc_prod if hc_prod[hc]['atomdict'].keys()==negtype.keys()]
                if not hc_list: #No match
                    raise Exception
                hc_prod2={hc:hc_prod[hc] for hc in hc_list} #Narrow down list of help compounds
                reac,prod,hcid,msg=balance(Rdata,Pdata,hc_prod=hc_prod2,coefflim=coefflim,addedspecies=addedspecies,
                                          addedhc=addedhc,hc_react=hc_react)
            except Exception:
                return update_rxn(Rdata,Pdata,hcrct=addedhc,rxnsmiles0=rxnsmiles0,msg='Imbalanced'+addedstr) 
            
            else:
                return update_rxn(Rdata,Pdata,reac=reac,prod=prod,hcprod=hcid,hc_prod=hc_prod2,hcrct=addedhc,rxnsmiles0=rxnsmiles0,
                                 msg=msg)
        else:
            return update_rxn(Rdata,Pdata,hcrct=addedhc,rxnsmiles0=rxnsmiles0,msg='Imbalanced'+addedstr) 

def balance_analogue(row,basic=True,balance=True,coefflim=5): #More reliable
    Rdata={}
    Pdata={}
    Rgtdata={}
    hc_prod={}
    hc_react={}
    Rdata=copy.deepcopy(row['Rdata'])
    Pdata=copy.deepcopy(row['Pdata'])
    Rgtdata=copy.deepcopy(row['Rgtdata'])
    hc_prod=copy.deepcopy(row['hc_prod'])
    hc_react=copy.deepcopy(row['hc_react'])
    numrefs=row['NumRefs']
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
        return balancerxn(Rdata,Pdata,Rgtdata=Rgtdata,rxnsmiles0=rxnsmiles0,numrefs=numrefs,coefflim=coefflim,
                      hc_prod=hc_prod,hc_react=hc_react)
                     

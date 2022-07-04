# %load ./GenTempl.py
import copy,itertools
from rdkit import Chem
from MainFunctions import molfromsmiles
from rdkit.Chem import rdChemReactions
from rdkit.Chem import Descriptors

# def replace_deuterated(smi):
#     return re.sub('\[2H\]', r'[H]', smi)

def gentemplate(analoguerxnsfinal,ncpus=16,restart=True,specificity='loose'):
    if ncpus>1:
        if restart:
            initray(num_cpus=ncpus)
        analoguerxnsfinaldis=mpd.DataFrame(analoguerxnsfinal)
    else:
        analoguerxnsfinaldis=analoguerxnsfinal
    templts=analoguerxnsfinaldis.apply(gen_template_row,specificity=specificity,axis=1,result_type='reduce')
    templtser=pd.Series(data=templts.values,index=templts.index) #Optional convert modin back to pandas
    templtdf=pd.DataFrame(data=templtser.tolist(),index=templtser.index,columns=['template','LHSdata','RHSdata','msg4','farfg','unusedprod'])
    analoguerxnstempl=copy.deepcopy(analoguerxnsfinal)
    analoguerxnstempl[['template','LHSdata','RHSdata','msg4','farfg','unusedprod']]=templtdf
    return analoguerxnstempl

def gen_template_row(row,specificity='loose'):
    LHSdata=copy.deepcopy(row.LHSdata)
    RHSdata=copy.deepcopy(row.RHSdata)
    specmap=copy.deepcopy(row.specmap)
    rnbmap=copy.deepcopy(row.rnbmap)
    outfrag=copy.deepcopy(row.outfrag)
    unusedanalogue=row.unusedanalogue
    return gen_template(LHSdata,RHSdata,specmap,rnbmap=rnbmap,outfrag=outfrag,unusedanalogue=unusedanalogue,specificity=specificity)

# def testrow(row):
#     LHSdata=copy.deepcopy(row.LHSdata)
#     hatomssuspect=set()
#     for analoguecompd in LHSdata:
#         LHSdata_=LHSdata[analoguecompd]
#         if 'templatefragidx' not in LHSdata_ or 'templatefrag' not in LHSdata_:
#             if 'reacfrag' not in LHSdata_ or not LHSdata_['reacfrag']:
#                 continue
#             for inst,fraginf in LHSdata_['reacfrag'].items():
#                 reacfragidx={atomidx for frag,matchidx in fraginf.items() for match in matchidx for atomidx in list(LHSdata_['fragloc'][inst][frag]['corrmatches'])[match]}
#                 if type(inst)==tuple:
#                     reacmol=molfromsmiles(LHSdata_['mappedsmiles'][inst[0]][inst[1]])
#                 else:
#                     reacmol=molfromsmiles(LHSdata_['mappedsmiles'][inst])
#                 hatomssuspect.update(gen_template_fragment(reacfragidx,reacmol))
#     return hatomssuspect
        

def gen_template_fragment(fragidx,mappedmol,specificity='loose',atomsymbols=OrderedDict({}),funcgroupmapnum=set()):
#     breakpoint()
    if not atomsymbols:
#         breakpoint()
        mappedmol=Chem.RemoveHs(mappedmol) #Hydrogens become implicit
        clear_isotope(mappedmol) #Redefining isotopes
        hatomsidx={atom.GetIdx() for atom in mappedmol.GetAtoms() if atom.GetSymbol()=='H'}
        if fragidx.issubset(hatomsidx):
            hatomsremaining={}
        else:
            hatomsremaining=fragidx.intersection(hatomsidx)
        if hatomsremaining: #Likely isotopes
            affectedatomidx={nb.GetIdx() for hatom in hatomsremaining for nb in mappedmol.GetAtomWithIdx(hatom).GetNeighbors()}
        else:
            affectedatomidx={}
        fragidx={idx for idx in fragidx if idx<len(mappedmol.GetAtoms()) if idx not in hatomsremaining}
        for atom in Chem.AddHs(mappedmol).GetAtoms(): #Functional group atoms described in detail, everything else has hydrogens stripped
#             breakpoint()
            numHs_=None
            degree_=None
            if atom.GetSymbol()=='H' and atom.GetIdx() not in hatomsidx:
                continue
            refatom=mappedmol.GetAtomWithIdx(atom.GetIdx())
            if affectedatomidx and atom.GetIdx() in affectedatomidx:
                Hslist=[nb.GetSymbol() for nb in atom.GetNeighbors()]
                numHs_=Hslist.count('H')
                degree_=atom.GetDegree()-numHs_
            if specificity=='loose':
                if funcgroupmapnum and atom.HasProp('molAtomMapNumber') and atom.GetAtomMapNum() in funcgroupmapnum: 
#                 if refatom.GetDegree()==1: #Terminal atom
                    symbol=get_strict_smarts_for_atom(refatom,numHs_=numHs_,degree_=degree_)
                else:
                    symbol=atom.GetSmarts(isomericSmiles=False) #General SMARTS to yield a generalizable template
            else:
                symbol=get_strict_smarts_for_atom(refatom,numHs_=numHs_,degree_=degree_)
            atomsymbols.update({atom.GetIdx():symbol})
    if specificity=='loose':
        frag=Chem.MolFragmentToSmiles(mappedmol,fragidx,atomSymbols=list(atomsymbols.values()))
    else: #Closer to rdchiral template
        frag=Chem.MolFragmentToSmiles(mappedmol,fragidx,atomSymbols=list(atomsymbols.values()),allHsExplicit=True,allBondsExplicit=True)   
    return frag,atomsymbols,fragidx

def resolvelargetemplate(fragidx,mol):
#     breakpoint()
    atomidxbound=set()
    for atom in [mol.GetAtomWithIdx(idx) for idx in fragidx]:
        if any([nb.GetIdx() not in fragidx for nb in atom.GetNeighbors()]) or atom.GetDegree()==1:
            atomidxbound.add(atom.GetIdx())
    if len(atomidxbound)>2:
        fragcurrlist=[]
        fragidxlist=[]
        fragidxold=copy.deepcopy(fragidx)
        comblist=list(itertools.combinations(atomidxbound,2))
        for comb in comblist:
            nblist=[nb.GetIdx() for nb in mol.GetAtomWithIdx(comb[0]).GetNeighbors()]
            if comb[1] in nblist:
                continue
            fragidx=set(Chem.GetShortestPath(mol,comb[0],comb[1]))
            fragidx=fragidx.union(fragidxold)
            fragcurr=Chem.MolFragmentToSmiles(mol,fragidx)
            fragidxlist+=[fragidx]
            fragcurrlist+=[fragcurr]
            if not len(fragcurr.split('.'))>1:
                return fragidx
        numsplitfrags=[len(fragcurr.split('.')) for fragcurr in fragcurrlist]
        minfragnum=min(numsplitfrags)
        minidx=[idx for idx,val in enumerate(numsplitfrags) if val==minfragnum]
        if len(minidx)>1:
            fragidxlist=[fragidxlist[idx] for idx in minidx]
            fragcurrlist=[fragcurrlist[idx] for idx in minidx]
            lenfrags=[len(fragcurr) for fragcurr in fragcurrlist]
            minfraglen=min(lenfrags)
            minidx=[idx for idx,val in enumerate(lenfrags) if val==minfraglen]
        fragidx=fragidxlist[minidx[0]]
#         breakpoint()
        return resolvelargetemplate(fragidx,mol)
    return fragidx
    
    
            
def gen_template(LHSdata,RHSdata,specmap,rnbmap={},outfrag={},unusedanalogue=[],specificity='loose'): #Invalid reactions outside fragments may yield general templates
    farfg=[]
    unusedprod=[]
    addstr=[]
    alloutfrag=[]
    funcgroupmapnum=set()
    
    for analoguecompd in LHSdata:
        LHSdata_=LHSdata[analoguecompd]
        if 'templatefragidx' not in LHSdata_ or 'templatefrag' not in LHSdata_:
            LHSdata_['templatefragidx']={}
            LHSdata_['templatefrag']={}
#         breakpoint()
        if 'reacfrag' not in LHSdata_ or not LHSdata_['reacfrag']: #Hydrogen will not be mapped and won't show up in reaction center
            if (analoguecompd not in unusedanalogue) and (analoguecompd not in outfrag): #Hydrogen
                for inst in LHSdata_['fragloc']:
                    if type(inst)==tuple:
                        mappedmol=molfromsmiles(LHSdata_['mappedsmiles'][inst[0]][inst[1]])
                    else:
                        mappedmol=molfromsmiles(LHSdata_['mappedsmiles'][inst])
                    fragidx={atom.GetIdx() for atom in mappedmol.GetAtoms()}
                    fragcurr,ratomsymbols,fragidx=gen_template_fragment(fragidx,mappedmol,specificity=specificity,atomsymbols=OrderedDict({}))
                    LHSdata_['templatefragidx'].update({inst:fragidx})
                    LHSdata_['templatefrag'].update({inst:fragcurr})
            elif analoguecompd in outfrag: # Reaction center completely outside fragments
                alloutfrag+=[analoguecompd]
            continue
        for inst,fraginf in LHSdata_['reacfrag'].items():
            reacfragidx={atomidx for frag,matchidx in fraginf.items() for match in matchidx for atomidx in list(LHSdata_['fragloc'][inst][frag]['corrmatches'])[match]}
            try:
                funcgroupids={atomidx for frag,matchidx in fraginf.items() for match in matchidx for atomidx in list(LHSdata_['fragloc'][inst][frag]['funcgroupids'])[match]}
            except IndexError:
                funcgroupids={}
            if outfrag and analoguecompd in outfrag: #Reaction center partially outside fragment
                combinedmap={**specmap,**rnbmap}
                reacfragidx=reacfragidx.union({combinedmap[mapnum][2] for mapnum in outfrag[analoguecompd]})
                addstr.append('Template expanded to include atoms in RC outside relevant fragments')
            if type(inst)==tuple:
                reacmol=molfromsmiles(LHSdata_['mappedsmiles'][inst[0]][inst[1]])
            else:
                reacmol=molfromsmiles(LHSdata_['mappedsmiles'][inst])
            if specificity=='loose':
                funcgroupmapnum.update({reacmol.GetAtomWithIdx(funcgroupid).GetAtomMapNum() for funcgroupid in funcgroupids if reacmol.GetAtomWithIdx(funcgroupid).HasProp('molAtomMapNumber')})
            reacfragcurr,ratomsymbols,reacfragidx=gen_template_fragment(reacfragidx,reacmol,specificity=specificity,atomsymbols=OrderedDict({}),funcgroupmapnum=funcgroupmapnum)
            if len(reacfragcurr.split('.'))>1:
                reacfragidx=resolvelargetemplate(reacfragidx,reacmol)
#                 reacfragidx=set(Chem.GetShortestPath(reacmol,min(reacfragidx),max(reacfragidx)))
                reacfragcurr,ratomsymbols,reacfragidx=gen_template_fragment(reacfragidx,reacmol,specificity=specificity,atomsymbols=ratomsymbols,funcgroupmapnum=funcgroupmapnum)
#             breakpoint()
            if len(reacfragcurr.split('.'))>1 and analoguecompd not in farfg: # Fragments involved in the reaction are too far away, and are decomposed into more than 1 species. Invalid template.
                farfg+=[analoguecompd]
            LHSdata_['templatefragidx'].update({inst:reacfragidx})
            LHSdata_['templatefrag'].update({inst:reacfragcurr})
            for atomidx in reacfragidx:
#                 breakpoint()
                RHSdata=updateRHSdata(RHSdata,analoguecompd,inst,atomidx,specmap)
    for prodid in RHSdata:
#         breakpoint()
        RHSdata_=RHSdata[prodid]
        if 'templatefragidx' not in RHSdata_ or 'templatefrag' not in RHSdata_:
            RHSdata_['templatefragidx']={}
            RHSdata_['templatefrag']={}
        if 'rxnatomidx' not in RHSdata_:
            if RHSdata_['formula']=='H2': # Template may NOT be balanced!
                for i,inst in enumerate(RHSdata_['mappedsmiles']):
                    if type(inst)==tuple:
                        for j,inst1 in enumerate(inst):
                            mappedmol=molfromsmiles(inst1)
                            prodfragidx={atom.GetIdx() for atom in mappedmol.GetAtoms()}
                            prodfragcurr,patomsymbols,prodfragidx=gen_template_fragment(prodfragidx,mappedmol,specificity=specificity,atomsymbols=OrderedDict({}),funcgroupmapnum=funcgroupmapnum)
                            RHSdata_['templatefragidx'].update({(i,j):prodfragidx})
                            RHSdata_['templatefrag'].update({(i,j):prodfragcurr})
                    else:
                        mappedmol=molfromsmiles(inst)
                        prodfragidx={atom.GetIdx() for atom in mappedmol.GetAtoms()}
                        prodfragcurr,patomsymbols,prodfragidx=gen_template_fragment(prodfragidx,mappedmol,specificity=specificity,atomsymbols=OrderedDict({}),funcgroupmapnum=funcgroupmapnum)
                        RHSdata_['templatefragidx'].update({i:prodfragidx})
                        RHSdata_['templatefrag'].update({i:prodfragcurr})
            elif prodid not in unusedprod:
                unusedprod+=[prodid]
            continue
        for inst,prodfragidx in RHSdata_['rxnatomidx'].items(): #Ideally need to keep looping
#             breakpoint()
            if type(inst)==tuple:
                prodmol=molfromsmiles(RHSdata_['mappedsmiles'][inst[0]][inst[1]])
            else:
                prodmol=molfromsmiles(RHSdata_['mappedsmiles'][inst])
            prodfragcurr,patomsymbols,prodfragidx=gen_template_fragment(prodfragidx,prodmol,specificity=specificity,atomsymbols=OrderedDict({}),funcgroupmapnum=funcgroupmapnum)
#             breakpoint()
            if len(prodfragcurr.split('.'))>1:
                prodfragidx=resolvelargetemplate(prodfragidx,prodmol)
#                 prodfragidx=set(Chem.GetShortestPath(prodmol,min(prodfragidx),max(prodfragidx)))
                prodfragcurr,patomsymbols,prodfragidx=gen_template_fragment(prodfragidx,prodmol,specificity=specificity,atomsymbols=patomsymbols,funcgroupmapnum=funcgroupmapnum)
                if len(prodfragcurr.split('.'))>1 and prodid not in farfg: # Product fragment is still split
                    farfg+=[prodid]
                prodfragmapidx={atom.GetAtomMapNum() for frag in prodfragcurr.split('.') for atom in Chem.MolFromSmarts(frag).GetAtoms()}
                reacfragmapidx={atom.GetAtomMapNum() for analoguecompd in LHSdata 
                                for reacfrag in LHSdata[analoguecompd]['templatefrag'].values()
                                for atom in Chem.MolFromSmarts(reacfrag).GetAtoms() if atom.HasProp('molAtomMapNumber')}
                if prodfragmapidx!=reacfragmapidx: #Extra atoms in product
#                     breakpoint()
                    missingmapidx=prodfragmapidx-reacfragmapidx
                    for missingmap in missingmapidx:
                        if missingmap not in specmap:
                            addstr.append('Mapped atom '+str(missingmap)+' in product '+str(prodid)+' cannot be traced back to reactant')
                            break
                        analoguecompd=specmap[missingmap][0]
                        inst=specmap[missingmap][1]
                        idx=specmap[missingmap][2]
                        LHSdata[analoguecompd]['templatefragidx'][inst].add(idx)
                        reacmol=molfromsmiles(LHSdata[analoguecompd]['mappedsmiles'][inst])
                        reacfragcurr=gen_template_fragment(LHSdata[analoguecompd]['templatefragidx'][inst],reacmol,specificity=specificity,atomsymbols=OrderedDict({}),funcgroupmapnum=funcgroupmapnum)[0]
                        LHSdata[analoguecompd]['templatefrag'][inst]=reacfragcurr
            RHSdata_['templatefragidx'].update({inst:prodfragidx})
            RHSdata_['templatefrag'].update({inst:prodfragcurr})
    if alloutfrag:
        msg4='Species '+', '.join([str(ID) for ID in alloutfrag])+' react completely outside relevant fragments'
        return 'Error',LHSdata,RHSdata,msg4,farfg,unusedprod
    if farfg:
        addstr.append('Reacting functional groups are too far away in '+'species '+', '.join([str(rctid) for rctid in farfg]))
    if unusedprod:
        addstr.append('Product '+', '.join([str(ID) for ID in unusedprod])+' is not produced from reacting fragments')
    if not addstr:
        msg4='Valid'
    else:
        msg4=', '.join(addstr)
    reacfrag=[reacfragcurr for analoguecompd in LHSdata for reacfragcurr in LHSdata[analoguecompd]['templatefrag'].values() if reacfragcurr]
    prodfrag=[prodfragcurr for prodid in RHSdata for prodfragcurr in RHSdata[prodid]['templatefrag'].values() if prodfragcurr]
    return '>>'.join(['.'.join(reacfrag),'.'.join(prodfrag)]),LHSdata,RHSdata,msg4,farfg,unusedprod

def updateRHSdata(RHSdata,analoguecompd,inst,atomidx,specmap):
    for val in specmap.values(): #Don't consider rnbmap as unmapped reactant atoms will not appear on RHS
        if val[0]==analoguecompd and val[1]==inst and val[2]==atomidx:
    #                         breakpoint()
            prodid=val[3]
            prodinst=val[4]
            prodatomidx=val[5]
            if 'rxnatomidx' not in RHSdata[prodid].keys():
                RHSdata[prodid]['rxnatomidx']={prodinst:{prodatomidx}}
                break
            elif prodinst in RHSdata[prodid]['rxnatomidx']:
                RHSdata[prodid]['rxnatomidx'][prodinst].add(prodatomidx)
                break
            else:
                RHSdata[prodid]['rxnatomidx'].update({prodinst:{prodatomidx}})
                break
    return RHSdata
    

def matchtemplaterow(row,samplerxn):
    LHSdata=copy.deepcopy(row.LHSdata)
    template=copy.deepcopy(row.template)
    return matchtemplate(samplerxn,LHSdata,template)
    

def matchtemplate(samplerxn,LHSdata,template,exact=False): #Pass in similarly sized fragments
#     breakpoint()
    samplerxn_=rdChemReactions.ReactionFromSmarts(samplerxn,useSmiles=True)
    samplercts=[Chem.MolToSmiles(samplerct) for samplerct in samplerxn_.GetReactants()]
    samplerctscounter=Counter(samplercts)
    sampleprods=[Chem.MolToSmiles(sampleprod) for sampleprod in samplerxn_.GetProducts()]
    sampleprodscounter=Counter(sampleprods)
    templaterxn=rdChemReactions.ReactionFromSmarts(template,useSmiles=True)
    rdChemReactions.RemoveMappingNumbersFromReactions(templaterxn)
    templprods=[Chem.MolToSmiles(prod) for prod in templaterxn.GetProducts()]
    templprodscounter=Counter(templprods)
    valid=True
    for ID0 in LHSdata:
        if LHSdata[ID0]['templatefrag']:
            if 'reacfrag' not in LHSdata[ID0] or not LHSdata[ID0]['reacfrag']:
                fraglist=[Chem.MolToSmiles(Chem.MolFromSmarts(frag)) for inst in LHSdata[ID0]['fragloc'] for frag in LHSdata[ID0]['fragloc'][inst]]
            else:
                fraglist=[Chem.MolToSmiles(Chem.MolFromSmarts(frag)) for inst in LHSdata[ID0]['reacfrag'] for frag in LHSdata[ID0]['reacfrag'][inst]]
            fraglistcounter=Counter(fraglist)
            if set(fraglistcounter.keys()).issubset(set(samplerctscounter.keys())):
                rem=Counter()
                rem.update(samplerctscounter)
                rem.subtract(fraglistcounter)
                samplerctscounter={i:rem[i] for i in rem if rem[i]>0}
            else:
                valid=False
                break
    if samplerctscounter:
        valid=False
    if not valid:
        return valid
    rem=Counter()
    rem.update(sampleprodscounter)
    rem.subtract(templprodscounter)
    sampleprodscounter2={i:rem[i] for i in rem if rem[i]>0}
    templprodscounter2={i:rem[i] for i in rem if rem[i]<0}
    if sampleprodscounter2:
        if not exact:
            for sampleprod in sampleprodscounter2:
                for templprod in templprodscounter2:
                    match=findfrag(templprod,sampleprod,fragment=True,returnindices=False)
                    if match:
                        sampleprodscounter2[sampleprod]-=templprodscounter2[templprod]
                if sampleprodscounter2[sampleprod]>0:
                    valid=False
                    return valid
    return valid


def get_strict_smarts_for_atom(atom,numHs_=None,degree_=None,charge_=None):
    '''
    For an RDkit atom object, generate a SMARTS pattern that
    matches the atom as strictly as possible, taken from rdChiral
    '''
    USE_STEREOCHEMISTRY=True
    symbol = atom.GetSmarts()
#     if atom.GetSymbol() == 'H':
#         symbol = '[#1]'
    if numHs_ is not None:
        numHs=numHs_
    else:
        numHs=atom.GetTotalNumHs()
    if degree_ is not None:
        degree=degree_
    else:
        degree=atom.GetDegree()
    if charge_ is not None:
        charge=charge_
    else:
        charge=atom.GetFormalCharge()

    if '[' not in symbol:
        symbol = '[' + symbol + ']'

    # Explicit stereochemistry - *before* H
    if USE_STEREOCHEMISTRY:
        if atom.GetChiralTag() != Chem.rdchem.ChiralType.CHI_UNSPECIFIED:
            if '@' not in symbol:
                # Be explicit when there is a tetrahedral chiral tag
                if atom.GetChiralTag() == Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW:
                    tag = '@'
                elif atom.GetChiralTag() == Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW:
                    tag = '@@'
                if ':' in symbol:
                    symbol = symbol.replace(':', ';{}:'.format(tag))
                else:
                    symbol = symbol.replace(']', ';{}]'.format(tag))

    if 'H' not in symbol:
        H_symbol = 'H{}'.format(numHs)
        # Explicit number of hydrogens: include "H0" when no hydrogens present
        if ':' in symbol: # stick H0 before label
            symbol = symbol.replace(':', ';{}:'.format(H_symbol))
        else:
            symbol = symbol.replace(']', ';{}]'.format(H_symbol))
      
    # Explicit degree
    if ':' in symbol:
        symbol = symbol.replace(':', ';D{}:'.format(degree))
    else:
        symbol = symbol.replace(']', ';D{}]'.format(degree))

    # Explicit formal charge
    if '+' not in symbol and '-' not in symbol:
        charge_symbol = '+' if (charge >= 0) else '-'
        charge_symbol += '{}'.format(abs(charge))
        if ':' in symbol: 
            symbol = symbol.replace(':', ';{}:'.format(charge_symbol))
        else:
            symbol = symbol.replace(']', ';{}]'.format(charge_symbol))

    return symbol

def clear_isotope(mol):
    [a.SetIsotope(0) for a in mol.GetAtoms()]


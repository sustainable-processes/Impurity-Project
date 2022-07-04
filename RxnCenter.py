#%load ./RxnCenter.py

from MainFunctions import molfromsmiles
from FindFunctionalGroups import identify_functional_groups as IFG
import copy
from rdkit import Chem #Importing RDKit

def reactioncenter(analoguerxnsassignedfilt,ncpus=16,restart=True):
    '''
    Finds reaction center of a reaction
    '''
    if ncpus>1:
        if restart:
            initray(num_cpus=ncpus)
        analoguerxnsassignedfiltdis=mpd.DataFrame(analoguerxnsassignedfilt)
    else:
        analoguerxnsassignedfiltdis=analoguerxnsassignedfilt
    rxncenterser=analoguerxnsassignedfiltdis.apply(getrxncenterrow,axis=1,result_type='reduce')
    rxncenterser=pd.Series(data=rxncenterser.values,index=rxncenterser.index) #Optional convert modin back to pandas
    rxncenterdf=pd.DataFrame(data=rxncenterser.tolist(),index=rxncenterser.index,columns=['specmap','rnbmap','rxncentermapnum','rxncenter'])
    analoguerxnscent=copy.deepcopy(analoguerxnsassignedfilt)
    analoguerxnscent[['specmap','rnbmap','rxncentermapnum','rxncenter']]=rxncenterdf
    return analoguerxnscent

def validreactioncenter(analoguerxnscentfilt,ncpus=16,restart=True):
    if ncpus>1:
        if restart:
            initray(num_cpus=ncpus)
        analoguerxnscentfiltdis=mpd.DataFrame(analoguerxnscentfilt)
    else:
        analoguerxnscentfiltdis=analoguerxnscentfilt
    validrxncenterser=analoguerxnscentfiltdis.apply(validrxncenterrow,axis=1,result_type='reduce')
    validrxncenterser=pd.Series(data=validrxncenterser.values,index=validrxncenterser.index)
    validrxncenterdf=pd.DataFrame(data=validrxncenterser.tolist(),index=validrxncenterser.index,columns=['LHSdata','msg3','outfrag','outfg','outneighbor','unusedanalogue'])
    analoguerxnsvalid=copy.deepcopy(analoguerxnscentfilt)
    analoguerxnsvalid[['LHSdata','msg3','outfrag','outfg','outneighbor','unusedanalogue']]=validrxncenterdf
    return analoguerxnsvalid
    
    

#%% Reaction center identification

def bond_to_label(bond):
    '''
    This function takes an RDKit bond and creates a label describing
    the beginning atom (atomic number and atom map number) and ending atom
    (atomic number and atom map number)

    Parameters
    ----------
    bond : RDKit bond

    Returns
    -------
    output: str
        String label containing information about beginning and end atom linked
        by bond
        

    '''
    a1_label = str(bond.GetBeginAtom().GetAtomicNum())
    a2_label = str(bond.GetEndAtom().GetAtomicNum())
    if bond.GetBeginAtom().HasProp('molAtomMapNumber'):
        a1_label += bond.GetBeginAtom().GetProp('molAtomMapNumber')
    if bond.GetEndAtom().HasProp('molAtomMapNumber'):
        a2_label += bond.GetEndAtom().GetProp('molAtomMapNumber')
    atoms = sorted([a1_label, a2_label])
    # return '{}{}{}'.format(atoms[0], bond.GetSmarts(), atoms[1])
    return '{}{}'.format(atoms[0], atoms[1])



def atoms_are_different(atom1, atom2, level = 1,usesmarts=True):
    '''
    Compares two RDKit atoms based on common properties (Smarts, number of
    bonded hydrogens, charge, degree, radical electrons, neighbors etc.).
    From rdchiral.

    Parameters
    ----------
    atom1 : RDKit atom
    atom2 : RDKit atom
    level : Optional
            The default is 1.
    usesmarts : Optional (boolean)
                Option to use or ignore differences in atom smarts

    Returns
    -------
    True if the two atoms are different and False if they are similar.

    '''
    # import pdb; pdb.set_trace()
    if usesmarts:
        if atom1.GetSmarts() != atom2.GetSmarts(): return True # should be very general
    if atom1.GetAtomicNum() != atom2.GetAtomicNum(): return True # must be true for atom mapping
    if atom1.GetTotalNumHs() != atom2.GetTotalNumHs(): return True
    if atom1.GetFormalCharge() != atom2.GetFormalCharge(): return True
    if atom1.GetDegree() != atom2.GetDegree(): return True
    if atom1.GetIsAromatic() != atom2.GetIsAromatic(): return True
    #if atom1.IsInRing() != atom2.IsInRing(): return True
    if atom1.GetNumRadicalElectrons() != atom2.GetNumRadicalElectrons(): return True
    # TODO: add # pi electrons like ICSynth? Account for chirality
    # Check bonds and nearest neighbor identity
    if level >= 1:
        bonds1 = sorted([bond_to_label(bond) for bond in atom1.GetBonds()])
        bonds2 = sorted([bond_to_label(bond) for bond in atom2.GetBonds()])
        if bonds1 != bonds2: return True
    return False

def getrxncenterrow(row):
    '''
    Applies getrxncenter to each row of a dataframe

    '''
    LHSdata=row.LHSdata
    RHSdata=row.RHSdata
    specmap,rnbmap=parsemap(LHSdata,RHSdata)
    rxncenter,msg=getrxncenter(specmap,LHSdata,RHSdata)
    return specmap,rnbmap,rxncenter,msg

def storeatommap(mappedsmiles,specid=0,idx=0,atommap={},neighbormap={}):
    mappedmol=molfromsmiles(mappedsmiles) #Hydrogens will never be mapped
    for atom in mappedmol.GetAtoms():
        if atom.HasProp('molAtomMapNumber'):
            mnum = atom.GetAtomMapNum()
            atommap[mnum] = (specid,idx,atom.GetIdx())
#             breakpoint()
            neighbors=set([nb.GetIdx() for nb in atom.GetNeighbors() if not nb.HasProp('molAtomMapNumber')])
            if neighbors:
                if neighbormap:
                    startidx=max([int(key.split('n')[1]) for key in neighbormap])+1
                else:
                    startidx=0
                for i,nb in enumerate(neighbors,start=startidx):
                    neighbormap['n'+str(i)]=(specid,idx,nb)
    return atommap,neighbormap

    

def parsemap(LHSdata,RHSdata):
    '''
    Takes reaction and product data, analyzes mapped smiles
    and generates a dictionary with keys as atom map numbers, and values as a tuple 
    containing the reactant ID, reactant instance, reactant atom index, product ID,
    product instance,and product atom index

    Parameters
    ----------
    LHSdata : dict
        Reactant data formatted as output of getcompdict with 'mappedsmiles' key under
        each reactant ID with value as mapped smiles string
    RHSdata : dict
        Product data formatted as output of getcompdict with 'mappedsmiles' key under
        each product ID with value as mapped smiles string

    Returns
    -------
    specmap : dict
        Dictionary with keys as atom map numbers, and values as tuples 
        containing the owning reactant ID, reactant instance, reactant atom index, product ID,
        product instance,and product atom index

    '''
#     breakpoint()
    reactantMap = {}
    rnbmap={}
    for rctid in LHSdata:
        for ridx,mappedsmiles in enumerate(LHSdata[rctid]['mappedsmiles']):
            if type(mappedsmiles)==tuple: #Mixture detected
                for ridx2,mappedsmile in enumerate(mappedsmiles):
                    reactantMap,rnbmap=storeatommap(mappedsmile,specid=rctid,idx=(ridx,ridx2),atommap=reactantMap,
                                                    neighbormap=rnbmap)
            else:
                reactantMap,rnbmap=storeatommap(mappedsmiles,specid=rctid,idx=ridx,atommap=reactantMap,
                                                neighbormap=rnbmap)          
    productMap={}
    pnbmap={}
    for prodid in RHSdata:
        for pidx,mappedsmiles in enumerate(RHSdata[prodid]['mappedsmiles']):
            if type(mappedsmiles)==tuple: #Mixture detected
                for pidx2,mappedsmile in enumerate(mappedsmiles):
                    productMap,pnbmap=storeatommap(mappedsmile,specid=prodid,idx=(pidx,pidx2),atommap=productMap,
                                                   neighbormap=pnbmap)
            else:
                productMap,pnbmap=storeatommap(mappedsmiles,specid=prodid,idx=pidx,atommap=productMap,
                                               neighbormap=pnbmap)
#         if pnbmap: Error will never happen
            
    specmap={mapnum:tuple([elem for map_ in [reactantMap,productMap] for elem in map_[mapnum]]) for mapnum in reactantMap}              
    return specmap,rnbmap


def getrxncenter(specmap,LHSdata,RHSdata):
    '''
    Takes specmap (mapdict), reactant data, product data, optionally neighbor information and returns a set of
    changed atoms, as well as a message indicating if reaction center exists
    Note reaction center will ignore hydrogen reactants as these aren't mapped in products

    Parameters
    ----------
    specmap : dict
        Output of parsemap() function. Dictionary with keys as atom map numbers, and values as tuples 
        containing the owning reactant ID, reactant instance, reactant atom index, product ID,
        product instance,and product atom index
    LHSdata : dict
        DESCRIPTION.
    RHSdata : dict
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.
    msg : TYPE
        DESCRIPTION.

    '''
    
    changed_atoms=[]
    changed_mapnum=[]
    msg=True

    for mapnum,val in specmap.items():
        rctid=val[0]
        ridx=val[1]
        ratomidx=val[2]
        prodid=val[3]
        pidx=val[4]
        patomidx=val[5]
        
        if type(ridx)==tuple: #Mixture detected
            ridx1=ridx[0]
            ridx2=ridx[1]
            rtempl=molfromsmiles(LHSdata[rctid]['mappedsmiles'][ridx1][ridx2])
        else:
            rtempl=molfromsmiles(LHSdata[rctid]['mappedsmiles'][ridx])
        if type(pidx)==tuple: #Mixture detected
            pidx1=pidx[0]
            pidx2=pidx[1]
            ptempl=molfromsmiles(RHSdata[prodid]['mappedsmiles'][pidx1][pidx2])
        else:
            ptempl=molfromsmiles(RHSdata[prodid]['mappedsmiles'][pidx])
        ratom=rtempl.GetAtomWithIdx(ratomidx)
        patom=ptempl.GetAtomWithIdx(patomidx)
        if mapnum not in changed_mapnum:
            if atoms_are_different(patom,ratom):
                changed_atoms.append(ratom)
                changed_mapnum.append(mapnum)
    if not changed_mapnum:
        msg=False
    return set(changed_mapnum),msg 

#%% Valid reaction center

def validrxncenterrow(row):
    specmap=row.specmap
    rxncenter=row.rxncentermapnum
    LHSdata=copy.deepcopy(row.LHSdata)
    rnbmap=row.rnbmap
    return validrxncenter(specmap,rxncenter,LHSdata,rnbmap=rnbmap) 

def validrxncenter(specmap,rxncenter,LHSdata,rnbmap={}):
#     breakpoint()
    RCs=[copy.copy(rxncenter),set(rnbmap.keys())]
    mapdicts=[specmap,rnbmap]
    outfrag={}
#     outfrag=[] #Species not directly involved in the reaction at carrier fragments
    outfg={} #Species not directly involved in the reaction at functional group
    outneighbor={} #Species not involved at a neighboring atom to a functional group
#     unusedanalogue={analoguecompd for analoguecompd in LHSdata}
    for RC,mapdict in zip(RCs,mapdicts):
        for changemapnum in RC:
#             breakpoint()
            assigned=False
            fg=False
            nb=False
            rctid=mapdict[changemapnum][0]
            ridx=mapdict[changemapnum][1]
            idxr=mapdict[changemapnum][2]
            if type(ridx)==tuple:
                rctmol=molfromsmiles(LHSdata[rctid]['mappedsmiles'][ridx[0]][ridx[1]])
            else:
                rctmol=molfromsmiles(LHSdata[rctid]['mappedsmiles'][ridx])
            neidx=[atom.GetIdx() for atom in rctmol.GetAtomWithIdx(idxr).GetNeighbors()]
            if 'reacfrag' not in LHSdata[rctid]:
                LHSdata[rctid].update({'reacfrag':{}})
            reacfrag=LHSdata[rctid]['reacfrag']
#             unusedanalogue-={rctid}
            for carrierfrag,loc in LHSdata[rctid]['fragloc'][ridx].items():
                for matchidx,corrmatch in enumerate(loc['corrmatches']):
    #                 breakpoint()
                    if {idxr}.issubset(corrmatch):
                        if loc['funcgroupids']: # Some carrier fragments don't have functional groups
                            if {idxr}.issubset(loc['funcgroupids'][matchidx]):
                                fg=True #changed atom lies within functional group
                                nb=True #To bypass neighbor check
                            # New code for neighbors
                            else:
                                if any([atomidx in loc['funcgroupids'][matchidx] for atomidx in neidx]):
                                    nb=True
                        else:
                            fg=True
                            nb=True
                        if ridx not in reacfrag:
                            reacfrag.update({ridx:{carrierfrag:[matchidx]}})
                        else: #already initialized
                            if ridx in reacfrag:
                                if carrierfrag in reacfrag[ridx]:
                                    if matchidx not in reacfrag[ridx][carrierfrag]:
                                        reacfrag[ridx][carrierfrag].extend([matchidx])
                                else:
                                    reacfrag[ridx].update({carrierfrag:[matchidx]})
                            else:
                                reacfrag.update({ridx:{carrierfrag:[matchidx]}})
                        assigned=True
#             breakpoint()
            if not assigned: #changed atom not in any fragments
                if rctid not in outfrag:
                    outfrag.update({rctid:[changemapnum]})
                else:
                    outfrag[rctid].extend([changemapnum])
#                     outfrag+=[rctid]
            if not fg: #changed atom not in any functional groups
                if rctid not in outfg:
                    outfg.update({rctid:[changemapnum]})
                else:
                    outfg[rctid].extend([changemapnum])
            if not nb: #changed atom neighbors functional groups
                if rctid not in outneighbor:
                    outneighbor.update({rctid:[changemapnum]})
                else:
                    outneighbor[rctid].extend([changemapnum])
    #%% New (SUSPECT)
    unusedanalogue=[ID for ID in LHSdata if LHSdata[ID]['formula']!='H2' if 
                    ('reacfrag' not in LHSdata[ID] or set(LHSdata[ID]['reacfrag'].keys())!=set(LHSdata[ID]['fragloc'].keys()))
                   if ID not in outfrag if ID not in outfg if ID not in outneighbor]
    #%% New
    addstr=[]
    if outfrag:
        addstr.append('Species '+', '.join([str(ID) for ID in outfrag])+' not reacting at carrier fragment')
    if outfg:
        addstr.append('Species '+', '.join([str(ID) for ID in outfg])+' not reacting at functional group')
    if unusedanalogue:
        addstr.append('Species '+', '.join([str(ID) for ID in unusedanalogue])+' does not participate in reaction')
    if not addstr:
        msg='Valid'
    else:
        msg=', '.join(addstr)
    return LHSdata,msg,outfrag,outfg,outneighbor,unusedanalogue
        

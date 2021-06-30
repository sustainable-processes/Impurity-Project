#%load ./RxnCenter.py

from MainFunctions import molfromsmiles
from FindFunctionalGroups import identify_functional_groups as IFG
import copy
from rdkit import Chem #Importing RDKit

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
    res=parsemap(LHSdata,RHSdata)
    rxncenter,msg=getrxncenter(res,LHSdata,RHSdata)
    return res,rxncenter,msg
                
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
    res : dict
        Dictionary with keys as atom map numbers, and values as tuples 
        containing the owning reactant ID, reactant instance, reactant atom index, product ID,
        product instance,and product atom index

    '''
    res = {}
    reactantMap = {}
    for rctid in LHSdata:
        for ridx,mappedsmiles in enumerate(LHSdata[rctid]['mappedsmiles']):
            rtempl=molfromsmiles(mappedsmiles) #Add hydrogens to keep index constant for later substructure search (no need)
            for atom in rtempl.GetAtoms():
                if atom.HasProp('molAtomMapNumber'):
                    mnum = atom.GetAtomMapNum()
                    reactantMap[mnum] = (rctid,ridx,atom.GetIdx())
                else:
                    continue      
    for prodid in RHSdata:
        for pidx,mappedsmiles in enumerate(RHSdata[prodid]['mappedsmiles']):
            ptempl=molfromsmiles(mappedsmiles)
            for atom in ptempl.GetAtoms():
                if atom.HasProp('molAtomMapNumber'):
                    mnum = atom.GetAtomMapNum()
                    if mnum in reactantMap:
                        rctid,ridx,idxr = reactantMap[mnum]
                        res[mnum] = (rctid,ridx,idxr,prodid,pidx,atom.GetIdx())
                else:
                    continue         
    return res

def getrxncenter(res,LHSdata,RHSdata):
    '''
    Takes res (mapdict), reactant data, product data, and returns a set of
    changed atoms, as well as a message indicating if reaction center exists
    Note reaction center will ignore hydrogen reactants as these aren't mapped in products

    Parameters
    ----------
    res : dict
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
    changed_mapidx=[]
    msg=True

    for mapidx,val in res.items():
        rctid=val[0]
        ridx=val[1]
        idxr=val[2]
        prodid=val[3]
        pidx=val[4]
        idxp=val[5]

        rtempl=molfromsmiles(LHSdata[rctid]['mappedsmiles'][ridx])
        ptempl=molfromsmiles(RHSdata[prodid]['mappedsmiles'][pidx])

        ratom=rtempl.GetAtomWithIdx(idxr)
        patom=ptempl.GetAtomWithIdx(idxp)
        if mapidx not in changed_mapidx:
            if atoms_are_different(patom,ratom):
                changed_atoms.append(ratom)
                changed_mapidx.append(mapidx)
    if not changed_mapidx:
        msg=False
    return set(changed_mapidx),msg 


#%% Substructure matching
        
def get_matches(mol,patt,checkresults=True):
    '''
    Returns atom indices for substructure matches of a pattern in a molecule
    
        Parameters
        ----------
        mol : RDKit mol
            Molecule to check pattern
        patt : RDKit mol
            Pattern fragment
        checkresults: bool
            Optional, True if strict match with the pattern is needed including hydrogens

        Returns
        -------
       corr_matches: Set
           Set of tuples containing atom indices for every pattern match, after verification
           that correct functional group is present. Returned only if checkresults is True. 
       
       funcgroupids: List
           List of atom ids corresponding to functional groups in each correct pattern match
           
        matches: Set
            Set of tuples containing atom indices for every pattern match. Returned only if checkresults is False.
    '''
    # import pdb; pdb.set_trace()
    matches=mol.GetSubstructMatches(patt)
    if not matches:
        return False,False
    elif checkresults:
#         breakpoint()
        funcgroupmol=IFG(mol) #Functional groups of RDKit reactant
        funcgrouppatt=IFG(patt) #Functional groups of carrier fragment
        funcids=set() #Store functional groups that are of the same type as the carrier fragment
        for funcgroup in funcgrouppatt:
            matchtype=[molgroup for molgroup in funcgroupmol if molgroup.type==funcgroup.type] #change to .atoms if not working
            for molgroup in matchtype:
                if not any([atoms_are_different(mol.GetAtomWithIdx(atomid),patt.GetAtomWithIdx(pattid),usesmarts=False) 
                            for atomid,pattid in zip(molgroup.atomIds,funcgroup.atomIds)]): #BUGGY
                    funcids.update({atomid for atomid in molgroup.atomIds})

        corr_matches=[match for match in matches if set(match).intersection(funcids)]
        # funcgroupids={atomid for match in corr_matches for atomid in set(match).intersection(funcids)}
        funcgroupids=[set(match).intersection(funcids) for match in corr_matches]
        return corr_matches,funcgroupids
    else:
        return matches

#%% Valid reaction center

def validrxncenterrow(row,fragdict):
    res=row.mapdict
    rxncenter=row.rxncentermapnum
    LHSdata=copy.deepcopy(row.LHSdata)
    print(row.index)
    return validrxncenter(res,rxncenter,LHSdata,fragdict) 

def validrxncenter(res,rxncenter,LHSdata,fragdict):
#     breakpoint()
    RC=copy.copy(rxncenter)
    nofg=[] #Species with no functional groups
    noanalogue=[] #Species that don't have correct/desired functional group
    outfrag=[] #Species not directly involved in the reaction at carrier fragments
    outfg={} #Species not directly involved in the reaction at functional group
    outneighbor={} #Species not involved at a neighboring atom to a functional group
#     unused=[] #Species not directly involved in the reaction
#     extratom=[] #Species involved at functional group but reacts at other atoms outside as well
    for analoguecompd in LHSdata:
        fragloc={}
        for idx,cleanmol0 in enumerate(LHSdata[analoguecompd]['cleanmol']):
#         cleanmol=Chem.AddHs(LHSdata[analoguecompd]['cleanmol'])#molfromsmiles(LHSdata[analoguecompd]['smiles']))
            cleanmol=Chem.AddHs(cleanmol0)
            if not IFG(cleanmol): #Compound has no active fragment/functional group (including dihydrogen)
                nofg+=[analoguecompd]
#             continue
            matched=False
            for carrierfrag in LHSdata[analoguecompd]['querycompds']:
                patt=Chem.MolFromSmiles(fragdict[carrierfrag][0].smiles,sanitize=False)
                patt.UpdatePropertyCache(strict=False)
#                 patt=fragdict[carrierfrag][0].mol
                if analoguecompd in nofg: #Compound has no active fragment/functional group (including dihydrogen)
                    corr_matches=get_matches(cleanmol,patt,checkresults=False)
                    corr_matches=[match for match in corr_matches]
                    funcgroupids=[{elem for elem in match} for match in corr_matches]
                else:
                    corr_matches,funcgroupids=get_matches(cleanmol,patt) #funcgroupids refers to active fragment
                if not corr_matches:
                    fragloc.update({idx:{}})
                    continue
                corr_matches=[tuple((idx for idx in co if idx<len(molfromsmiles(LHSdata[analoguecompd]['smiles']).GetAtoms()))) for co in corr_matches]
                matched=True
                if idx not in fragloc.keys():
                    fragloc.update({idx:{carrierfrag:{'corrmatches':corr_matches,'funcgroupids':funcgroupids}}})
                else:
                    fragloc[idx].update({carrierfrag:{'corrmatches':corr_matches,'funcgroupids':funcgroupids}})
        if not matched:
            noanalogue+=[analoguecompd] #Species not analogue
        LHSdata[analoguecompd]['fragloc']=fragloc
        LHSdata[analoguecompd]['reacfrag']={}
#     return LHSdata,fragloc,matched
        
    unusedanalogue={analoguecompd for analoguecompd in LHSdata}
    for changemapnum in RC:
#         breakpoint()
        assigned=False
#         reacfrag={}
        fg=False
        nb=False
        rctid=res[changemapnum][0]
        ridx=res[changemapnum][1]
        idxr=res[changemapnum][2]
        rctmol=molfromsmiles(LHSdata[rctid]['mappedsmiles'][ridx])
        neidx=[atom.GetIdx() for atom in rctmol.GetAtomWithIdx(idxr).GetNeighbors()]
        reacfrag=LHSdata[rctid]['reacfrag']
        unusedanalogue-={rctid}
        
        for carrierfrag,loc in LHSdata[rctid]['fragloc'][ridx].items():
            for matchidx,corrmatch in enumerate(loc['corrmatches']):
#                 breakpoint()
                if {idxr}.issubset(corrmatch):
                    if {idxr}.issubset(loc['funcgroupids'][matchidx]):
                        fg=True #changed atom lies within functional group
                        nb=True #To bypass neighbor check
                    # New code for neighbors
                    else:
                        if any([atomidx in loc['funcgroupids'][matchidx] for atomidx in neidx]):
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
        if not assigned and (rctid not in noanalogue) and (rctid not in nofg): #changed atom not in any fragments
            if rctid not in outfrag:
                outfrag+=[rctid]
            
        if not fg and (rctid not in noanalogue) and (rctid not in nofg): #changed atom not in any functional groups
            if rctid not in outfg:
                outfg.update({rctid:[changemapnum]})
            else:
                outfg[rctid].extend([changemapnum])
        if not nb and (rctid not in noanalogue) and (rctid not in nofg): #changed atom neighbors functional groups
            if rctid not in outneighbor:
                outneighbor.update({rctid:[changemapnum]})
            else:
                outneighbor[rctid].extend([changemapnum])
    #%% New (SUSPECT)
    if unusedanalogue: #Reactants not involved in the reaction
        unusedanalogue={ID for ID in unusedanalogue if LHSdata[ID]['formula']!='H2'}
    #%% New
    addstr=[]
    if nofg:
        addstr.append("No functional groups in "+'species '+', '.join([str(ID) for ID in nofg]))
    if noanalogue:
        addstr.append('Species '+', '.join([str(ID) for ID in noanalogue])+' not analogue')
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
    return LHSdata,msg,nofg,noanalogue,outfrag,outfg,outneighbor,unusedanalogue
        

          
# =============================================================================
#FROM COLEY (OPTIONAL)
# def get_special_groups(mol):
#     '''
#     This retrieves special groups from input molecules that should be part of
#     the final template fragments. Only for specific templates. If general templates are
#     preferred, then don't call this function. NOT IN USE AT THE MOMENT.
# 
#     Parameters
#     ----------
#     mol : TYPE
#         DESCRIPTION.
# 
#     Returns
#     -------
#     None.
# 
#     '''
# 
# # Define templates
#     group_templates = [
#         'C(=O)Cl', # acid chloride
#         'C(=O)[O;H,-]', # carboxylic acid
#         '[$(S-!@[#6])](=O)(=O)(Cl)', # sulfonyl chloride
#         '[$(B-!@[#6])](O)(O)', # boronic acid
#         '[$(N-!@[#6])](=!@C=!@O)', # isocyanate
#         '[N;H0;$(N-[#6]);D2]=[N;D2]=[N;D1]', # azide
#         'O=C1N(Br)C(=O)CC1', # NBS brominating agent
#         'C=O', # carbonyl
#         'ClS(Cl)=O', # thionyl chloride
#         '[Mg][Br,Cl]', # grinard (non-disassociated)
#         '[#6]S(=O)(=O)[O]', # RSO3 leaving group
#         '[O]S(=O)(=O)[O]', # SO4 group
#         '[N-]=[N+]=[C]', # diazo-alkyl
#         ]
#     # Build list
#     groups = []
#     for template in group_templates:
#         matches = mol.GetSubstructMatches(Chem.MolFromSmarts(template))
#         groups.extend(list(matches))
#     return groups
# 
#     # group_templates = [
#     #     (range(3), '[OH0,SH0]=C[O,Cl,I,Br,F]',), # carboxylic acid / halogen
#     #     (range(3), '[OH0,SH0]=CN',), # amide/sulfamide
#     #     (range(4), 'S(O)(O)[Cl]',), # sulfonyl chloride
#     #     (range(3), 'B(O)O',), # boronic acid/ester
#     #     ((0,), '[Si](C)(C)C'), # trialkyl silane
#     #     ((0,), '[Si](OC)(OC)(OC)'), # trialkoxy silane, default to methyl
#     #     (range(3), '[N;H0;$(N-[#6]);D2]-,=[N;D2]-,=[N;D1]',), # azide
#     #     (range(8), 'O=C1N([Br,I,F,Cl])C(=O)CC1',), # NBS brominating agent
#     #     (range(11), 'Cc1ccc(S(=O)(=O)O)cc1'), # Tosyl
#     #     ((7,), 'CC(C)(C)OC(=O)[N]'), # N(boc)
#     #     ((4,), '[CH3][CH0]([CH3])([CH3])O'), #
#     #     (range(2), '[C,N]=[C,N]',), # alkene/imine
#     #     (range(2), '[C,N]#[C,N]',), # alkyne/nitrile
#     #     ((2,), 'C=C-[*]',), # adj to alkene
#     #     ((2,), 'C#C-[*]',), # adj to alkyne
#     #     ((2,), 'O=C-[*]',), # adj to carbonyl
#     #     ((3,), 'O=C([CH3])-[*]'), # adj to methyl ketone
#     #     ((3,), 'O=C([O,N])-[*]',), # adj to carboxylic acid/amide/ester
#     #     (range(4), 'ClS(Cl)=O',), # thionyl chloride
#     #     (range(2), '[Mg,Li,Zn,Sn][Br,Cl,I,F]',), # grinard/metal (non-disassociated)
#     #     (range(3), 'S(O)(O)',), # SO2 group
#     #     (range(2), 'N~N',), # diazo
#     #     ((1,), '[!#6;R]@[#6;R]',), # adjacency to heteroatom in ring
#     #     ((2,), '[a!c]:a:a',), # two-steps away from heteroatom in aromatic ring
#     #     #((1,), 'c(-,=[*]):c([Cl,I,Br,F])',), # ortho to halogen on ring - too specific?
#     #     #((1,), 'c(-,=[*]):c:c([Cl,I,Br,F])',), # meta to halogen on ring - too specific?
#     #     ((0,), '[B,C](F)(F)F'), # CF3, BF3 should have the F3 included
#     # ]
#     # # Build list
#     # groups = []
#     # for (add_if_match, template) in group_templates:
#     #     matches = mol.GetSubstructMatches(Chem.MolFromSmarts(template), useChirality=False)
#     #     for match in matches:
#     #         add_if = []
#     #         for pattern_idx, atom_idx in enumerate(match):
#     #             if pattern_idx in add_if_match:
#     #                 add_if.append(atom_idx)
#     #         groups.append((add_if, match))
#     # return groups
# 
# 
# 
# def get_fragments_for_changed_atoms(mapdict,changed_mapidx, category, radius = 0):
#     """
#     This builds fragments around changed atoms. NOT IN USE AT THE MOMENT.  
#     
#     Parameters
#     ----------
#     mapdict : TYPE
#         DESCRIPTION.
#     changed_mapidx : TYPE
#         DESCRIPTION.
#     category : Specify either reactant or product
#         Category to generate fragments for template smarts
#     radius : TYPE, optional
#         DESCRIPTION. The default is 0.
# 
#     Returns
#     -------
#     None.
# 
#     """
#     fragments=''
#     mol_done=[]
#     atoms_to_use=[] #List for fragment construction
#     symbols=[]
# 
#     for mapidx in changed_mapidx:
#         val=mapdict[mapidx]
#         prod=val[0]
#         react=val[1]
#         idxp=val[2]
#         idxr=val[3]
#         if category=='reactants':
#             mol=react
#             idx=idxr
#         else:
#             mol=prod
#             idx=idxp
# 
#         atoms_to_use.append(idx) #Adding changed atom to atoms list for fragment construction
# 
#         if mol not in mol_done:
#             for atom in mol.GetAtoms():
#                 symbol=atom.GetSmarts()
#                 if atom.GetTotalNumHs() == 0:
#                 # Be explicit when there are no hydrogens
#                     if ':' in symbol: # stick H0 before label
#                         symbol = symbol.replace(':', ';H0:')
#                     else: # stick before end
#                         symbol = symbol.replace(']', ';H0]')
#                     # print('Being explicit about H0!!!!')
#         if atom.GetFormalCharge() == 0:
#                 # Also be explicit when there is no charge
#             if ':' in symbol:
#                 symbol = symbol.replace(':', ';+0:')
#             else:
#                 symbol = symbol.replace(']', ';+0]')
#         if symbol !=atom.GetSmarts():
#                 symbol_replacements.append()
#                 symbols=[atom.GetSmarts() for atom in mol.GetAtoms()]
# 
#         mol_done.append(mol)
# 
#             # CUSTOM SYMBOL CHANGES
#         if atom.GetTotalNumHs() == 0:
#                 # Be explicit when there are no hydrogens
#             if ':' in symbol: # stick H0 before label
#                 symbol = symbol.replace(':', ';H0:')
#             else: # stick before end
#                 symbol = symbol.replace(']', ';H0]')
#                     # print('Being explicit about H0!!!!')
#         if atom.GetFormalCharge() == 0:
#                 # Also be explicit when there is no charge
#             if ':' in symbol:
#                 symbol = symbol.replace(':', ';+0:')
#             else:
#                 symbol = symbol.replace(']', ';+0]')
#         if symbol !=atom.GetSmarts():
#                 symbol_replacements.append()
# 
# =============================================================================

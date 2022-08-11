# %load ./MainFunctions.py
import ray
from rdkit import Chem #Importing RDKit
from rdkit.Chem import Draw #For drawing molecules/reactions
from rdkit.Chem import rdChemReactions #Reaction processing
from rdkit.Chem.Draw import rdMolDraw2D #Drawing 2D molecules/reactions
from IPython.display import SVG  #For SVG support
from PIL import Image #Working with images
import io
import os #Working with the OS
try:
    import pickle5 as pickle #Only if pickle doesn't work
except Exception:
    import pickle
import cairosvg
import shutil
from rdkit.Chem import PeriodicTable, GetPeriodicTable
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

#%% Generating mol files

def molfromsmiles(SMILES):
    """
    Converts a smiles string into a mol object for RDKit use.

    Args:
        SMILES (Str): Smiles string of molecule

    Returns:
        mol: RDKit mol object
    """
    mol=Chem.MolFromSmiles(SMILES)
    Chem.SanitizeMol(mol)
    mol.UpdatePropertyCache(strict=False)
    return mol

#%% Drawing reactions/mols


def mol_with_atom_index(mol):
    """
    Draw a molecule with atom index numbers (set by RDKit)

    Args:
        mol (RDKit mol): RDKit mol plain

    Returns:
        mol: RDKit mol with atom indices
    """
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(atom.GetIdx())
    return mol


def moveAtomMapsToNotes(m):
    '''
    Make note of the atom mapping

    Parameters
    ----------
    m : RDKit mol
      Mapped RDKit mol

    Returns
    -------
    RDKit mol with map numbers set as notes

    '''
    for at in m.GetAtoms():
        if at.GetAtomMapNum():
            at.SetProp("atomNote",str(at.GetAtomMapNum()))
            #at.SetAtomMapNum(0)

def drawReaction(rxn):
    '''
    Produces an SVG of a mapped reaction

    Parameters
    ----------
    rxn : String
        Reaction SMARTS string

    Returns
    -------
    SVG image

    '''
    trxn = rdChemReactions.ChemicalReaction(rxn)
    for m in trxn.GetReactants():
        moveAtomMapsToNotes(m)
    for m in trxn.GetProducts():
        moveAtomMapsToNotes(m)
    # if filetype=='svg':
#     d2d = rdMolDraw2D.MolDraw2DSVG(4000,500)
    d2d=rdMolDraw2D.MolDraw2DSVG(1000,300)
#     d2d.drawOptions().minFontSize = 20
    # else:
    #     d2d=rdMolDraw2D.MolDraw2DCairo(4000,500)
    d2d.DrawReaction(trxn,highlightByReactant=True)
    d2d.FinishDrawing()
    img=d2d.GetDrawingText()
    # if filetype=='svg':
    return SVG(img)
    # else:
    #     im=Image.open(io.BytesIO(img))
    #     return im
        # d2d.WriteDrawingText('text.png') #Only works with cairo, plus atom notes not included

def drawMol(m,filetype):
    '''
    Produces an SVG of RDKit mol
    
    Parameters
    ----------
    m : RDKit mol

    Returns
    -------
    SVG image

    '''
    if filetype=='svg':
        d=rdMolDraw2D.MolDraw2DSVG(400,300)
    else:
        d=rdMolDraw2D.MolDraw2DCairo(400,300)

    Draw.PrepareAndDrawMolecule(d,m)
    d.FinishDrawing()
    img=d.GetDrawingText()

    if filetype=='svg':
        return SVG(img)
    else:
        im=Image.open(io.BytesIO(img))
        return im
def highlightsubstruct(smi,pattlist=[]):
    m=molfromsmiles(smi)
    m=Chem.AddHs(m)
    hits_atsfin=[]
    hits_bondfin=[]
    for patt in pattlist:
        patt=Chem.MolFromSmarts(patt)
        hit_ats=[list(match) for match in m.GetSubstructMatches(patt)]
#                 hit_ats=[at for match in m.GetSubstructMatches(patt) for at in match]
        if not hit_ats:
            continue
        hit_bonds = []
        for bond in patt.GetBonds():
            for match in hit_ats:
                aid1 = match[bond.GetBeginAtomIdx()]
                aid2 = match[bond.GetEndAtomIdx()]
                hit_bonds.append(m.GetBondBetweenAtoms(aid1,aid2).GetIdx())
        hits_atsfin+=list(set([at for match in hit_ats for at in match]))
        hits_bondfin+=list(set(hit_bonds))
    colors=[(0.7,1,0.7)]
    atom_cols={}
    for i, at in enumerate(hits_atsfin):
        atom_cols[at] = colors[0]
    bond_cols = {}
    for i, bd in enumerate(hits_bondfin):
        bond_cols[bd] = colors[0]
    d = rdMolDraw2D.MolDraw2DSVG(500, 500) # or MolDraw2DCairo to get PNGs
    rdMolDraw2D.PrepareAndDrawMolecule(d, m, highlightAtoms=hits_atsfin,
                                        highlightBonds=hits_bondfin,highlightAtomColors=atom_cols,highlightBondColors=bond_cols)
    d.FinishDrawing()
    return SVG(d.GetDrawingText())

#%% Utility Functions
        
class CustomError(Exception):
    '''
    For error handling
    '''
    pass

def openpickle(filename):
    '''
    Takes a pickled file and loads it

    Parameters
    ----------
    filename : str
        File path of .pickle file

    Returns
    -------
    loadfile : Python object
        

    '''
    infile=open(filename,'rb')
    infile.seek(0)
    loadfile=pickle.load(infile)
    infile.close()
    return loadfile
    # with open('filename', 'rb') as handle: #Faster way
    # return pickle.load(handle)

def writepickle(pkl,directory,filename):
    '''
    Parameters
    ----------
    pkl : Python object
    
    filename : String
        Directory to save in (exclude extension)

    Returns
    -------
    None.

    '''
    if not os.path.exists(directory):
        os.makedirs(directory)
    with open(directory+filename+'.pickle', 'wb') as handle:
        pickle.dump(pkl, handle, protocol=pickle.HIGHEST_PROTOCOL)


def getlist(refdict,key):
    '''
    Gets list from a reference dictionary, given a key.
    Assumes refdict is a dictionary within a dictionary.
    
    Parameters
    ----------
    refdict : Dictionary
        Reference dictionary with keys as Reaxys ID
    key : String
        The reaction property that need to be retrieved

    Returns
    -------
    list: List of elements

    '''
    return [rxn[key] for rxn in refdict.values() if key in rxn.keys()]


def writetofile(rxnimg,directory):
    '''
    Takes a reaction image and saves it to a specified directory (usually after calling drawReaction)

    Parameters
    ----------
    rxnimg : Reaction sketch
    
    directory : String
        File path to write the image to

    Returns
    -------
    None.

    '''
    open(directory,'w').write(rxnimg.data)

def delcontents(directory):
    '''
    Deletes contents of a directory
    
    Parameters
    ----------
    directory : String
        Directory to delete

    Returns
    -------
    None.
    
    '''
    for filename in os.listdir(directory):
        file_path = os.path.join(directory, filename)
        try:
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.unlink(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)
        except Exception as e:
            print('Failed to delete %s. Reason: %s' % (file_path, e))



def convSVGtoPNG(filename, filenameout):
    '''
    Converts SVG images to PNG images for easier copy pasting (eg./ to Word)

    Parameters
    ----------
    filename : String
        File name of SVG image 
        
    filenameout : String
        File name of PNG image

    Returns
    -------
    None.

    '''
    cairosvg.svg2png(url=filename+".svg", write_to=filenameout+".png")

#%% Analyzing mol files

def atomtypes(mol):
    """
    Generates an atom type dictionary with counts for a given mol file. Also
    returns overall change of mol

    Parameters
    ----------
    mol : RDKit mol

    Returns
    -------
    typedict : dict
        Dictionary with keys as atom element and value as count
    charge : int
        Overall charge of supplied molecule

    """

    typedict={}
    mol2=Chem.AddHs(mol) #Assumes hydrogens are added
    charge=0
    for atom in mol2.GetAtoms():
        elem=PeriodicTable.GetElementSymbol(GetPeriodicTable(),atom.GetAtomicNum())
        charge+=atom.GetFormalCharge()
        if elem not in typedict.keys():
            typedict[elem]=1
        else:
            typedict[elem]+=1
    return typedict, charge

def getcompdict(ID=1,mol=None,smiles=None,formula=None,FragDB=None):
    """
    Wrapper for atomtypes. ID must be input whether Reaxys ID or random number. Default is 1. 
    Smiles and mol can also be specified. Can be converted to an object instance if intended.
    
    Args:
        ID (int, optional): ID of compound/mol. Defaults to 1.
        mol (RDKit mol, optional): RDKit mol. Defaults to None.
        smiles (str, optional): mol smiles. Defaults to None.
        formula (str, optional): mol formula. Defaults to None.
        FragDB (Pandas dataframe, optional): Pandas dataframe to find smiles of ID. Defaults to None.

    Raises:
        CustomError: Raises error if smiles/mol not specified, only ID specified
        without reference fragment dataframe. Also if anything specified
        is not valid/cannot be processed.

    Returns:
        dict: Compound dictionary with ID as key and another dictionary with the
        following keys: 'atomdict' (Output of atomtypes, dictionary with atom
        elements and count), 'charge' (Overall charge of mol),'smiles'
        (smiles of mol),'formula' (chemical formula),'count' (instance number, always 1)
    """

    if smiles is None and mol is None:
        if FragDB is None:
            raise CustomError("Please provide a reference fragment database")
        try:
            smiles=locrecord(ID,FragDB,smiles=True)
        except Exception as e:
            raise CustomError("Compound "+str(ID)+" isn't in dataframe")
        try:
            mol=molfromsmiles(smiles)
        except Exception as e:
            raise CustomError("Compound "+str(ID)+" smiles is invalid")
    elif smiles:
        try:
            mol=molfromsmiles(smiles)
        except Exception as e:
            raise CustomError("Compound "+str(ID)+" smiles is invalid")
    if formula is None:
        formula=CalcMolFormula(mol)
    atomdata=atomtypes(mol)
    compddict={ID:{'atomdict':atomdata[0],'charge':atomdata[1],'smiles':smiles,'formula':formula,'count':1}}
    return compddict

#%% Database functions

def locrecord(ID,DB,smiles=False,fragsmarts=False,fragsmiles=False,mixture=False):

    '''
    Returns desired information by searching fragment database by ID

    Parameters
    ----------
    ID : int
        Substance ID to search for
    DB : pandas dataframe
        Fragment database with primary index as fragment smiles and secondary index as substance ID
    smiles : bool, optional
        Specify True if only smiles needs to be returned. The default is False.
    fragsmarts : bool, optional
        Specify True if only fragment smarts need to be returned. The default is False.
    fragsmiles : bool, optional
        Specify True if only fragment smiles needs to be returned . The default is False.
    mixture : bool, optional
        Specify True if only mixture presence needs to be returned. The default is False.

    Returns
    -------
    Output : pandas row or str depending on input preference. If everything is False, row will be returned

    '''

    record=DB.xs(ID,level=1)
    if smiles:
        return record.Smiles.unique()[0]
    if fragsmarts:
        return record.FragmentSmarts.values
    if fragsmiles:
        return record.index.values
    if mixture:
        return record['>1 Compound'].unique()[0]
    return record

#%% Joining smiles

def getfragments(chemlist,smiles=False,ref=None):
    '''
    Concatenates smiles strings in a list separated by '.' If ID's are given instead,
    looks up IDs in reference dictionary or database.

    Parameters
    ----------
    chemlist : list
        List of smiles strings OR substance IDs
    ref : dict/database (Optional if list of smiles, mandatory if list of IDs)
        If database provided should be multiindexed with smiles string attached to the ID as key. 
        At the moment only fragment database allowed.
    smiles: bool
        Indicate if chemlist is smiles or not

    Returns
    -------
    frag : str
        String joining all smiles given in/calculated from chemlist with '.' operator

    '''
    if chemlist==[]:
        raise CustomError('ERROR: NO CHEMICALS IN LIST PROVIDED')
    if smiles: #List of smiles provided
        return '.'.join(chemlist)
    else: #List of IDs given
        if ref is None:
            raise CustomError("Please supply reference (fragment) database multiindexed with ID as second level and containing a Smiles column")
        try:
            frag='.'.join([locrecord(ID,ref,smiles=True) for ID in chemlist])
        except Exception:
            raise CustomError('One or more compounds missing from database')
        else:
            return frag
#%% Parallel processing functions

def initray(restart=True,num_cpus=16,log_to_driver=False):
    """
    Initiates cluster of CPUs for parallel processing. If restart is True, will restart cluster if already running.

    Args:
        restart (bool, optional): Controls if ray cluster already exists, will restart. Defaults to True.
        num_cpus (int, optional): Number of CPUs to run in parallel. Defaults to 16.
        log_to_driver (bool, optional): Controls if more fine-grained status information on parallel execution is outputted to console. Defaults to False.
    """
    if restart:
        ray.shutdown()
    ray.init(num_cpus=num_cpus,log_to_driver=log_to_driver)


def chunks(lst, s=None,k=None):
    """Yield successive s-sized chunks or k chunks from lst"""
    if s is not None:
        return [lst[i:i+s] for i in range(0,len(lst),s)]
    elif k is not None:
        n=len(lst)
        return [lst[i * (n // k) + min(i, n % k):(i+1) * (n // k) + min(i+1, n % k)] for i in range(k)]

from rdkit import Chem

# help compounds
# website to draw structure and compute SMILES: 
# http://www.cheminfo.org/flavor/malaria/Utilities/SMILES_generator___checker/index.html
# or ChemDraw
hc_smilesDict = {
                 'H2': '[H][H]', 
                 'O2': 'O=O', 
                 'N2': 'N#N',
    
                 'Na+': '[Na+]',
                 'K+': '[K+]',
                 'NO2+': 'O=[N+]=O',
                 'SO3-': 'O=S(O)[O-]',

                 'CH4': 'C',
                 'H2O': 'O', 
                 'HCl': 'Cl',
                 'HBr': 'Br',
                 'HI': 'I',
                 'HF': 'F',
                 'PH3': 'P',
                 'NH3': 'N',
                 'H2S': 'S',
                 
                 'H3PO3': 'O=P(O)O',
                 'H2SO4': 'O=S(=O)(O)O',
                 'HNO3': 'O=N(=O)O',
                 'HClO3': 'O=Cl(=O)O',
                 
                 'CH3OH': 'CO',
                 'CH3CH2OH': 'CCO', # ethanol
                 'CH3CH2CH2OH': 'CCCO',
                 'CH3CH2CH2CH2OH': 'CCCCO',
                 
                 'CH3COOH': 'CC(=O)O',
                 'CH3CH2COOH': 'CCC(=O)O',
                 'CH3CH2CH2COOH': 'CCCC(=O)O',
    
                 'CH3OCH2COOH': 'COCC(=O)O',
                 
                 'C6H6': 'c1ccccc1', # benzene
                 'C6H5CH3': 'Cc1ccccc1', # toluene
                 'C6H5OH':'Oc1ccccc1', # phenol
                 'C6H5Cl': 'Clc1ccccc1'            
                }

hc_molDict = {hc: Chem.AddHs(Chem.MolFromSmiles(hc_smilesDict[hc])) for hc in hc_smilesDict}
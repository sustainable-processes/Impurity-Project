# -*- coding: utf-8 -*-
"""
Created on Tue Sep 22 12:10:08 2020

@author: E0014
"""

import cirpy # CIRpy 

# Check SMILES string
# cirpy.Molecule('p-nitrochlorobenzene','smiles')
    
s1=cirpy.resolve('(CH3)(CH2)C(O)OC(O)(CH2)(CH3)','smiles')
s2=cirpy.resolve('CCC(=O)OC(=O)CC','formula')


def getiupacnames(smiles):
    iupacnames=[]
    for smles in smiles:
        iupacnames.append(cirpy.resolve(smles,'iupac_name'))
    return iupacnames



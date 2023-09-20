#!/usr/bin/python
"""
find a compound by name, smiles, or formula 

"""
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Draw

def mol_infos(compound):
    mols_list,imgs_list,infos_list = [],[],[]
    for i in range(len(compound)):
        smiles  = compound[i].canonical_smiles
        formula = compound[i].molecular_formula
        weight  = compound[i].molecular_weight
        name    = compound[i].iupac_name
        synonyms= compound[i].synonyms
        infos = {
            "smiles":smiles,
            "formula":formula,
            "weight":weight,
            "name":name,
            "synonyms":synonyms
        }
        mol = Chem.MolFromSmiles(smiles)
        mols_list.append(mol)
        img = Draw.MolToImage(mol)
        imgs_list.append(img)
        infos_list.append(infos)
        # img.save("pvp.png")
        # img.show()

    return mols_list,imgs_list,infos_list

def is_smiles(smiles_str):
    try:
        molecule = Chem.MolFromSmiles(smiles_str)
        if molecule is not None:
            return True
        else:
            return False
    except:
        return False

def get_type_compounds(name):
    if is_smiles(name):
        compound = pcp.get_compounds(name,"smiles")
    else:
        compound = pcp.get_compounds(name,"name")

    return compound


if __name__ == '__main__':
    
    compound = get_type_compounds(name="pvp")
    mols_list,imgs_list,infos_list = mol_infos(compound)
    print(mols_list)
    print(infos_list[0]["smiles"])


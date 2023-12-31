from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.DataStructs.cDataStructs import TanimotoSimilarity

def calculate_tanimoto(smiles1, smiles2):
    # Convert the SMILES strings to RDKit molecule objects
    mol1 = Chem.MolFromSmiles(smiles1)
    mol2 = Chem.MolFromSmiles(smiles2)

    # Generate the fingerprints for each molecule
    fp1 = AllChem.GetMorganFingerprint(mol1, 2)
    fp2 = AllChem.GetMorganFingerprint(mol2, 2)

    # Calculate the Tanimoto coefficient
    tanimoto_coeff = TanimotoSimilarity(fp1, fp2)
    return tanimoto_coeff

with open('SMILES.txt') as file:
    lines = [line.rstrip() for line in file]

index = 0
with open('Enamine_Diverse_REAL_drug-like_48.2M_cxsmiles.cxsmiles') as file:
    for line in file:
        if index > 0:
            index2 = 0
            while index2 < len(lines):
                print(calculate_tanimoto(lines[index2], line.split()[0].strip()))
                index2 += 1
        index += 1

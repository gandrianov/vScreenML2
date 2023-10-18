from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdMolDescriptors import CalcFractionCSP3, CalcTPSA
from rdkit.Chem.Lipinski import NumHAcceptors, NumHDonors
from rdkit.Chem.rdFreeSASA import CalcSASA, classifyAtoms
from rdkit.Chem.Crippen import MolLogP
from rdkit.Chem.rdchem import AtomMonomerType

class RDKitCalculator:
    def __init__(self):
        pass

    def CalcFCsp3(self, mol):
        return CalcFractionCSP3(mol)

    def CalcNumHAcceptors(self, mol):
        return NumHAcceptors(mol)

    def CalcNumHDonors(self, mol):
        return NumHDonors(mol)

    def CalcMolLogP(self, mol):
        return MolLogP(mol)

    def CalcTPSA(self, mol):
        return CalcTPSA(mol)

def calculate_features(ligand_mol):
   
    from rdkit.Chem import rdMolDescriptors 

    features = {}
    descriptor_names = list(rdMolDescriptors.Properties.GetAvailableProperties())
    del descriptor_names[descriptor_names.index("NumAtomStereoCenters")]
    del descriptor_names[descriptor_names.index("NumUnspecifiedAtomStereoCenters")]

    get_descriptors = rdMolDescriptors.Properties(descriptor_names)
    descriptors = list(get_descriptors.ComputeProperties(ligand_mol))

    features = {k:v for k, v in zip(descriptor_names, descriptors)}

    rdkit_calculator = RDKitCalculator()
    #features["FCsp3"] = rdkit_calculator.CalcFCsp3(ligand_mol)
    features["NumHAcceptors"] = rdkit_calculator.CalcNumHAcceptors(ligand_mol)
    features["NumHDonors"] = rdkit_calculator.CalcNumHDonors(ligand_mol)
    features["MolLogP"] = rdkit_calculator.CalcMolLogP(ligand_mol)
    #features["TPSA"] = rdkit_calculator.CalcTPSA(ligand_mol)

    return features

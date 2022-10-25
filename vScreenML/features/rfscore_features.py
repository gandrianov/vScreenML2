"""Summary
"""
from oddt.toolkits.rdk import Molecule
from oddt.scoring.descriptors import close_contacts_descriptor


class RFScoreCalculator:

    """Summary

    Attributes:
        cutoff (TYPE): Description
        ligand_atoms (TYPE): Description
        protein_atoms (TYPE): Description
    """

    def __init__(self,
                 ligand_atoms=[6, 7, 8, 9, 15, 16, 17, 35, 53],
                 protein_atoms=[6, 7, 8, 16],
                 cutoff=12):
        """Summary

        Args:
            ligand_atoms (list, optional): Atomic numbers of ligand atoms allowed for calculations
            protein_atoms (list, optional): Atomic numbers of protein atoms allowed for calculations
            cutoff (int, optional): Distance cutoff for determining ligand-protein contacts
        """
        self.cutoff = cutoff
        self.ligand_atoms = ligand_atoms
        self.protein_atoms = protein_atoms

    def GetScore(self, ligand_mol, protein_mol):
        """Summary

        Args:
            ligand_mol (TYPE): Ligand structure
            protein_mol (TYPE): Protein structure

        Returns:
            dict: Return dictionary of number of contacts between `ligand_atoms` and `protein_atoms`
        """
        ligand_mol = Molecule(ligand_mol)
        protein_mol = Molecule(protein_mol)

        descriptors = close_contacts_descriptor(protein_mol,
                                                cutoff=self.cutoff,
                                                protein_types=self.protein_atoms,
                                                ligand_types=self.ligand_atoms,
                                                aligned_pairs=False)

        features = descriptors.build(ligand_mol).tolist()[0]
        return dict(zip(descriptors.titles, features))



def calculate_features(ligand_mol, protein_mol):
    
    rfscore_calculator = RFScoreCalculator()
    return rfscore_calculator.GetScore(ligand_mol, protein_mol)

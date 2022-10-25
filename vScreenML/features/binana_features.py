from oddt.toolkits.extras.rdkit import MolToPDBQTBlock
from vScreenML.external.binana.load_ligand_receptor import from_texts
from vScreenML.external.binana.interactions import get_all_interactions


class BinanaCalculator:
    def __init__(self):
        pass

    def CalcInteractions(self, ligand_pdbqt, protein_pdbqt):
        self.ligand, self.protein = from_texts(ligand_pdbqt, protein_pdbqt)
        self.interactions = get_all_interactions(self.ligand, self.protein)

    def GetSideFlexAlpha(self):
        return self.interactions["active_site_flexibility"]["counts"].get("SIDECHAIN_ALPHA", 0)

    def GetSideFlexBeta(self):
        return self.interactions["active_site_flexibility"]["counts"].get("SIDECHAIN_BETA", 0)

    def GetSideFlexOther(self):
        return self.interactions["active_site_flexibility"]["counts"].get("SIDECHAIN_OTHER", 0)

    def GetBackFlexAlpha(self):
        return self.interactions["active_site_flexibility"]["counts"].get("BACKBONE_ALPHA", 0)

    def GetBackFlexBeta(self):
        return self.interactions["active_site_flexibility"]["counts"].get("BACKBONE_BETA", 0)

    def GetBackFlexOther(self):
        return self.interactions["active_site_flexibility"]["counts"].get("BACKBONE_OTHER", 0)

    def GetPiPi(self):
        return self.interactions["pi_pi"]["counts"].get("pi_stacking", 0)

    def GetTStacking(self):
        for key in self.interactions["cat_pi"]["counts"].keys():
            self.interactions["pi_pi"]["counts"][key] = self.interactions["cat_pi"]["counts"][key]

        return self.interactions["pi_pi"]["counts"].get("T_stacking", 0)

    def GetCationPi(self):

        total_cat_pi = self.interactions["cat_pi"]["counts"]
        total_cat_pi = [v for v in total_cat_pi.values()]
        total_cat_pi = sum(total_cat_pi)

        return total_cat_pi

    def GetSaltBridge(self):

        total_salt_bridges = self.interactions["salt_bridges"]["counts"]
        total_salt_bridges = [v for v in total_salt_bridges.values()]
        total_salt_bridges = sum(total_salt_bridges)

        return total_salt_bridges

    def GetTotalElec(self):
        
        total_elec = self.interactions["electrostatic_energies"]["counts"].values()
        total_elec = sum(total_elec)

        return total_elec

    def GetTotalHBond(self):
        
        total_hbond = self.interactions["hydrogen_bonds"]["counts"].values()
        total_hbond = sum(total_hbond)

        return total_hbond

    def GetTotalHphobics(self):

        total_hphobics = self.interactions["hydrophobics"]["counts"].values()
        total_hphobics = sum(total_hphobics)

        return total_hphobics


def calculate_features(ligand_mol, protein_mol):

    features = {}

    # create PDBQt objects from RDKit mol
    ligand_pdbqt = MolToPDBQTBlock(ligand_mol, flexible=True, addHs=False, computeCharges=True)
    protein_pdbqt = MolToPDBQTBlock(protein_mol, flexible=False, addHs=False, computeCharges=True)

    binana_calculator = BinanaCalculator()
    
    binana_calculator.CalcInteractions(ligand_pdbqt, protein_pdbqt)

    features["SideFlexAlpha"] = binana_calculator.GetSideFlexAlpha()
    features["SideFlexBeta"] = binana_calculator.GetSideFlexBeta()
    features["SideFlexOther"] = binana_calculator.GetSideFlexOther()
    features["BackFlexAlpha"] = binana_calculator.GetBackFlexAlpha()
    features["BackFlexBeta"] = binana_calculator.GetBackFlexBeta()
    features["BackFlexOther"] = binana_calculator.GetBackFlexOther()
    features["PiPi"] = binana_calculator.GetPiPi()
    features["TStacking"] = binana_calculator.GetTStacking()
    features["CationPi"] = binana_calculator.GetCationPi()
    features["SaltBridge"] = binana_calculator.GetSaltBridge()
    features["TotalElec"] = binana_calculator.GetTotalElec()
    features["TotalHBond"] = binana_calculator.GetTotalHBond()
    features["TotalHphobics"] = binana_calculator.GetTotalHphobics()

    return features



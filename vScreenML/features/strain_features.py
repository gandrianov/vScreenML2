import pandas as pd
from vScreenML.utils.strain_utils import *

class StrainCalculator():

    def __init__(self, library_xml=None):

        if library_xml is None:
            dirname = "/".join(__file__.split("/")[:-2])
            dirname = dirname + "/external/STRAIN/"
            xml = dirname + "TL_2.1_VERSION_6.xml"

        self.library = get_torsion_library(xml)

    def calculate_energy(self, mol):

        data = []

        positions = mol.GetConformer().GetPositions()

        for class_name, patterns in self.library.items():
            for p in patterns:

                matches = mol.GetSubstructMatches(p["smarts_mol"])
                for m in matches:

                    if len(m) > 4:
                        continue

                    a1_idx, a2_idx, a3_idx, a4_idx = m

                    a1_symbol = mol.GetAtomWithIdx(a1_idx).GetSymbol()
                    a4_symbol = mol.GetAtomWithIdx(a4_idx).GetSymbol()

                    if a1_symbol == 'H' or a4_symbol == 'H':
                        continue

                    a1 = positions[a1_idx]
                    a2 = positions[a2_idx]
                    a3 = positions[a3_idx]
                    a4 = positions[a4_idx]

                    theta = dihedral(a1, a2, a3, a4)

                    d = process_match(theta, p)
                    d["Atom IDs"] = ",".join(map(str, m))
                    d["ID.1"] = m if m[1] > m[2] else m[::-1]
                    d["ID.2"] = d["ID.1"][1:3]
                    d["ID.1"] = ",".join(map(str, d["ID.1"]))
                    d["ID.2"] = ",".join(map(str, d["ID.2"]))

                    d["Theta"] = theta
                    d["Smarts"] = p["smarts"]
                    d["Class"] = "general" if class_name == "GG" else "specific"
                    d["Method"] = p["method"]
                    d["Rank"] = p["rank"]

                    data.append(d)

        data = pd.DataFrame(data)
        data = data.rename(columns={"energy":"strain_energy",
                                    "ci_lower":"strain_energy_ci_lower",
                                    "ci_upper":"strain_energy_ci_upper"})

        if data.shape[0] == 0:
            data = {"strain_energy":None, "strain_energy_ci_lower": None, "strain_energy_ci_upper":None}
            return data

        data = data.sort_values(["ID.1", "Rank"]).drop_duplicates(["ID.1"])
        data = data.sort_values(["ID.2", "Class", "strain_energy"], ascending=False)
        data = data.drop_duplicates(["ID.2"])
        
        data = data.sort_values("strain_energy", ascending=False)
        data = data.drop(["ID.1", "ID.2", "Rank"], axis=1)

        data_sums = data[["strain_energy", "strain_energy_ci_lower", "strain_energy_ci_upper"]].sum()

        return data_sums.to_dict()


def calculate_features(mol):

    calc = StrainCalculator()

    features = calc.calculate_energy(mol)
    
    return features
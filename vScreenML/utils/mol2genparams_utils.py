import sys
import argparse
from rdkit import Chem

def run_mol2params():

    root = __file__
    root = root[:root.rfind("/")]
    root = root[:root.rfind("/")]

    genpot_path = f"{root}/external/generic_potential"
    sys.path.append(genpot_path)

    from vScreenML.external.generic_potential.mol2genparams import main
    from vScreenML.external.generic_potential.BasicClasses import OptionClass

    option = OptionClass(sys.argv)
    main(option)

    # temporary solution

    try:
        mol = Chem.MolFromMol2File(option.opt.inputs[0], removeHs=False, sanitize=True)
        smiles = Chem.MolToSmiles(mol, kekuleSmiles=True)
    except:
        smiles = "None"

    with open(option.opt.prefix + ".params", "a") as fwr:
        fwr.write(f"\n#SMILES: {smiles}\n")
        fwr.close()
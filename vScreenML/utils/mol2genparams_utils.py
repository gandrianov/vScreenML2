import io
import sys
import argparse

from vScreenML.utils.rdkit_utils import read_mol2

def run_mol2params():

    root = __file__
    root = root[:root.rfind("/")]
    root = root[:root.rfind("/")]

    genpot_path = f"{root}/external/generic_potential"
    sys.path.append(genpot_path)

    from vScreenML.external.generic_potential.Molecule import MoleculeClass
    from vScreenML.external.generic_potential.BasicClasses import OptionClass

    option = OptionClass(sys.argv)

    input_mol2 = option.opt.inputs[0]
    output_pdb = option.opt.prefix + "_0001.pdb"
    output_params = option.opt.prefix + ".params"

    smiles = str(read_mol2(input_mol2))
    
    molecule = MoleculeClass(input_mol2, option)

    pdb = io.StringIO()
    params = io.StringIO()
    
    molecule.report_pdbfile(pdb, as_string=True)
    molecule.report_paramsfile(params, as_string=True)

    with open(output_pdb, "w") as fwr:
        fwr.write(pdb.getvalue())
        fwr.close()

    with open(output_params, "w") as fwr:
        fwr.write(f"#SMILES: {smiles}\n")
        fwr.write(params.getvalue())
        fwr.close()


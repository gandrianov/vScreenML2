import sys
import argparse

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
import argparse
from vScreenML.utils import pyrosetta_utils
from vScreenML.utils import rdkit_utils
from vScreenML.utils import xgboost_utils

from vScreenML.features import pyrosetta_features
from vScreenML.features import binana_features
from vScreenML.features import rfscore_features
from vScreenML.features import rdkit_features
from vScreenML.features import strain_features
from vScreenML.features import luna_features

from PocketDruggability.cmds import pocket_features

def minimize_pdb_complex():

    parser = argparse.ArgumentParser()

    parser.add_argument("-ligands", nargs="+", default=[], required=False,
                        help="PDB file of ligand generated by generic_potential/mol2genparams.py")
    parser.add_argument("-params", nargs="+", default=[], required=True,
                        help="Params file of ligand generated by generic_potential/mol2genparams.py")
    parser.add_argument("-protein", type=str, required=True,
                        help="PDB file of input protein structure")
    parser.add_argument("-output", type=str, required=True,
                        help="File name for minimized protein-ligand complex")
    parser.add_argument("-flexible-radii", type=float,
                        help="Radii in angstroms for defining flexible residues during minimization")

    args = parser.parse_args()

    # initialize variables
    ligand_pdbs = args.ligands
    ligand_params = args.params
    protein = args.protein
    output_pdb = args.output
    minimizeable_radii = args.flexible_radii
    
    # initialize pyrosetta
    params_filenames = " -extra_res_fa ".join(ligand_params)
    pyrosetta_utils.init_pyrosetta(f"-beta -extra_res_fa {params_filenames}")

    # compile pdb into one pdb strings
    ligand_pdb_strings = [open(p, "r").read() for p in ligand_pdbs]
    protein_pdb_string = open(protein, "r").read() 
    pdbstring = protein_pdb_string + "\n" + "\n".join(ligand_pdb_strings)

    # read input structure in rosetta
    bound_pose = pyrosetta_utils.load_pdbstring(pdbstring)
    
    # minimize bound pose
    if minimizeable_radii is not None:
        flexible_residues = pyrosetta_utils.get_around_residues(bound_pose, radii=minimizeable_radii)
    else:
        flexible_residues = []

    bound_pose = pyrosetta_utils.minimize_complex(bound_pose, flexible_residues=flexible_residues)

    # convert pose to pdb string
    bound_pdb = pyrosetta_utils.export_pdbstring(bound_pose)

    # write to file
    with open(output_pdb, "w") as fwr:
        fwr.write(bound_pdb)
        fwr.close()


def calculate_features():

    parser = argparse.ArgumentParser()
    parser.add_argument("-pdb-complex", type=str, required=True,
                        help="PDB file of protein/ligand complex minimized in Rosetta/PyRosetta")
    parser.add_argument("-ligand-params", nargs="+", default=[], required=True,
                        help="Params file of ligand generated by generic_potential/mol2genparams.py")
    parser.add_argument("-ligand-code", type=str, required=False, default="LG1",
                        help="Params file of ligand generated by generic_potential/mol2genparams.py")
    parser.add_argument("-output", type=str, required=False, default=None,
                        help="If file name specified, the calculated features will be written in CSV format")

    args = parser.parse_args()

    # initialize variables
    pdb_complex = args.pdb_complex
    ligand_params = args.ligand_params
    ligand_code = args.ligand_code
    output = args.output

    params_filenames = " -extra_res_fa ".join(ligand_params)
    params_strings = [open(p, "r").read() for p in ligand_params]

    # initialize pyrosetta
    pyrosetta_utils.init_pyrosetta(f"-beta -extra_res_fa {params_filenames}") 

    # read input structure in rosetta
    pdbstring = open(pdb_complex).read()
    bound_pose = pyrosetta_utils.load_pdbstring(pdbstring)

    # ligand residue number in bound structure
    ligand_resi = pyrosetta_utils.get_residue_index(bound_pose, name=ligand_code)

    # manipulation with bound structure
    unbound__pose = pyrosetta_utils.unbound_pose(bound_pose, residue_id=ligand_resi)
    ligand_pose, protein_pose = pyrosetta_utils.decompose_pose(bound_pose, residue_id=ligand_resi)

    # export disassebled ligand and protein to pdb strings
    ligand_pdb = pyrosetta_utils.export_pdbstring(ligand_pose)
    protein_pdb = pyrosetta_utils.export_pdbstring(protein_pose)

    # read ligand and protein structure in rdkit
    ligand_mol = rdkit_utils.load_pdbstring(ligand_pdb, params=params_strings)
    protein_mol = rdkit_utils.load_pdbstring(protein_pdb, params=params_strings)

    # features calculation
    
    pyrosetta_feats = pyrosetta_features.calculate_features(bound_pose, unbound__pose, ligand_pose, ligand_resi)
    binana_feats = binana_features.calculate_features(ligand_mol, protein_mol)
    rfscore_feats = rfscore_features.calculate_features(ligand_mol, protein_mol)
    rdkit_feats = rdkit_features.calculate_features(ligand_mol)
    strain_feats = strain_features.calculate_features(ligand_mol)
    luna_feats = luna_features.calculate_features(pdbstring, params_strings)
    pocket_feats = pocket_features(pdbstring.split("\n"), "LG1", 4.0)
    del pocket_feats["PDBid"]

    # buns
    polar_groups = rdkit_utils.find_polar_groups(ligand_mol) 
    burunsat_groups = rdkit_utils.calculate_burunsat_group(pyrosetta_feats["LigandInterfaceUnsat"], polar_groups)
    pyrosetta_feats.update(burunsat_groups)
    del pyrosetta_feats["LigandInterfaceUnsat"]

    # storage for features
    features = {"name": pdb_complex}
    features.update(pyrosetta_feats)
    features.update(binana_feats)
    features.update(rfscore_feats)
    features.update(rdkit_feats)
    features.update(strain_feats)
    features.update(luna_feats)
    features.update(pocket_feats)
    

    for k,v in features.items():
        if type(v) is not str:
            if v is not None:
                features[k] = str(round(v, 2))
            else:
                features[k] = ""

    formatted_output = ",".join(list(features.keys())) + "\n"
    formatted_output += ",".join(list(features.values()))

    if output is not None:
        with open(output, "w") as fwr:
            fwr.write(formatted_output)
            fwr.close()

    else:
        print(formatted_output)











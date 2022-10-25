import json, collections
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdchem import AtomMonomerType
from rdkit import RDLogger

RDLogger.DisableLog('rdApp.*')

root = __file__
root = root[:root.rfind("/")]
AMINOACID_CONNECTIVITY = json.load(open(f"{root}/aminoacid_connectivity.txt", "r"))


BOND_TYPES = {1:Chem.BondType.SINGLE,
              2:Chem.BondType.DOUBLE,
              3:Chem.BondType.TRIPLE,
              4:Chem.BondType.AROMATIC}


def parse_structure(mol):


    structure = {}

    for i, atom in enumerate(mol.GetAtoms()):
        pdbinfo = atom.GetPDBResidueInfo()

        atom_name = pdbinfo.GetName().strip()
        residue_num = pdbinfo.GetResidueNumber()
        residue_name = pdbinfo.GetResidueName().strip()
        chain = pdbinfo.GetChainId().strip()

        if chain not in structure:
            structure[chain] = {}

        if f"{residue_name}{residue_num}" not in structure[chain]:
            structure[chain][f"{residue_name}{residue_num}"] = {}

        if atom_name not in structure[chain][f"{residue_name}{residue_num}"]:
            structure[chain][f"{residue_name}{residue_num}"][atom_name] = i

    return structure


def assign_bonds(mol, mol_structure, params):

    mol_editable = Chem.RWMol(mol)

    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtomIdx()
        a2 = bond.GetEndAtomIdx()

        mol_editable.RemoveBond(a1, a2)

    for chain, residues in mol_structure.items():

        previous_residue_idx = None
        previous_residue_number = None

        for residue_idx in sorted(residues, key=lambda x: int(x[3:])):

            residue_name = residue_idx[:3]
            residue_number = int(residue_idx[3:])

            if residue_name in AMINOACID_CONNECTIVITY:

                if residue_name == "LYS":
                    mol_editable.GetAtomWithIdx(mol_structure[chain][residue_idx]["NZ"]).SetFormalCharge(1)

                if residue_name == "ARG":
                    mol_editable.GetAtomWithIdx(mol_structure[chain][residue_idx]["NE"]).SetFormalCharge(1)

                if residue_name == "HIS":
                    mol_editable.GetAtomWithIdx(mol_structure[chain][residue_idx]["ND1"]).SetFormalCharge(1)

                # add connections between backbone N of residue i and backbone C of residue i - 1

                if previous_residue_number is None:
                    pass

                elif previous_residue_number + 1 == residue_number:
                    mol_editable.AddBond(mol_structure[chain][residue_idx]["N"], mol_structure[chain][previous_residue_idx]["C"], BOND_TYPES[1])

                previous_residue_idx = residue_idx
                previous_residue_number = residue_number

                # add connections inside residue i

                for atom_name, atom_idx in residues[residue_idx].items():

                    # if amino acid is terminal then add terminal atoms

                    if atom_name in ["1H", "2H", "3H"]:
                        mol_editable.AddBond(atom_idx, mol_structure[chain][residue_idx]["N"], BOND_TYPES[1])
                        mol_editable.GetAtomWithIdx(mol_structure[chain][residue_idx]["N"]).SetFormalCharge(1)
                        continue

                    if atom_name == "OXT":
                        mol_editable.AddBond(atom_idx, mol_structure[chain][residue_idx]["C"], BOND_TYPES[1])
                        continue

                    elif atom_name in AMINOACID_CONNECTIVITY[residue_name]:
                        for atom2_name, bond_order in AMINOACID_CONNECTIVITY[residue_name][atom_name]:
                            if atom2_name in mol_structure[chain][residue_idx]:
                                mol_editable.AddBond(atom_idx, mol_structure[chain][residue_idx][atom2_name], BOND_TYPES[bond_order])

            elif residue_name in params:

                for atom_name, atom_idx in residues[residue_idx].items():

                    if atom_name in params[residue_name]:
                        for atom2_name, bond_order in params[residue_name][atom_name]:
                            if atom2_name in mol_structure[chain][residue_idx]:
                                mol_editable.AddBond(atom_idx, mol_structure[chain][residue_idx][atom2_name], BOND_TYPES[bond_order])

            else:
                print(f"Residue {residue_name} {residue_number} is not recognized")

    
    return mol_editable.GetMol()


def parse_params(paramsstring):

    name = None
    connectivity = []

    for s in paramsstring.split("\n"):

        if s.startswith("IO_STRING"):
            name = s.split(" ")[1]

        if s.startswith("BOND"):
            s = list(filter(lambda x: x != "", s.split(" ")))
            connectivity.append([s[1], s[2], int(s[3].replace("#ORGBND", ""))])

    nodes = set([cc for c in connectivity for cc in c[:2]])

    reshaped_connectivity = {}

    for node in nodes:
        reshaped_connectivity[node] = []
        
        to_delete = []
        for i, c in enumerate(connectivity):

            if node not in c:
                continue
            
            node2 = c[1] if c.index(node) == 0 else c[0]
            bond_order = c[-1]

            reshaped_connectivity[node].append([node2, bond_order])

            to_delete.append(i)

        for i in sorted(to_delete, reverse=True):
            del connectivity[i]

    reshaped_connectivity = {k:v for k, v in reshaped_connectivity.items() if len(v) != 0}

    return name, reshaped_connectivity


def load_pdbstring(pdbstring, params=[]):

    mol = Chem.MolFromPDBBlock(pdbstring, removeHs=False, proximityBonding=False, sanitize=False)
    mol_structure = parse_structure(mol)

    params = [parse_params(p) for p in params]
    params = dict(params)

    mol = assign_bonds(mol, mol_structure, params)
    
    mol.UpdatePropertyCache()
    Chem.SanitizeMol(mol)
   
    return mol



def prepare_amino_acids_connectivity():

    root = pkgutil.get_loader("pyrosetta").path
    root = root[:root.rfind("/")] + "/database/chemical/residue_type_sets/fa_standard/residue_types"
    
    aas = {}
    
    for f in glob.glob(f"{root}/l-caa/*params"):
        name, connectivity = parse_canonical_params(open(f,"r").read())
        aas[name] = connectivity

        if "DELOCALIZED" in set([aa[2] for aa in connectivity]):
            print(f, name)

    print(set([aa[2] for a in aas.values() for aa in a]))

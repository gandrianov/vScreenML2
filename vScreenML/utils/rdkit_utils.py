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
        insertion_code = pdbinfo.GetInsertionCode()

        if chain not in structure:
            structure[chain] = {}

        if f"{residue_name}{residue_num}{insertion_code}" not in structure[chain]:
            structure[chain][f"{residue_name}{residue_num}{insertion_code}"] = {}

        if atom_name not in structure[chain][f"{residue_name}{residue_num}{insertion_code}"]:
            structure[chain][f"{residue_name}{residue_num}{insertion_code}"][atom_name] = i

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

        for residue_idx in sorted(residues, key=lambda x: int(x[3:-1])):

            residue_name = residue_idx[:3]
            residue_number = int(residue_idx[3:-1])

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

                bond_props = params[residue_name][0]
                atom_props = params[residue_name][1]

                for atom_name, atom_idx in residues[residue_idx].items():
                    if atom_name in bond_props:
                        for atom2_name, bond_order in bond_props[atom_name]:
                            if atom2_name in mol_structure[chain][residue_idx]:
                                mol_editable.AddBond(atom_idx, mol_structure[chain][residue_idx][atom2_name], BOND_TYPES[bond_order])

                for atom_name, atom_idx in residues[residue_idx].items():
                    if atom_name in atom_props:
                        atom = mol_editable.GetAtomWithIdx(atom_idx)
                        atom.UpdatePropertyCache(strict=False)
                        atom.SetFormalCharge(atom_props[atom_name])
                        atom.SetNoImplicit(True)

            else:
                print(f"Residue {residue_name} {residue_number} is not recognized")

    return mol_editable.GetMol()


def extract_structure_properties(mol):

    BOND_TYPES = {1.0:1,
                  2.0:2,
                  3.0:3,
                  1.5:4}

    atom_properties = {}
    bond_properties = {}

    for i, a in enumerate(mol.GetAtoms()):
        formal_charge = a.GetFormalCharge()
        atom_properties[i] = formal_charge

    for i, b in enumerate(mol.GetBonds()):
        begin_idx = b.GetBeginAtomIdx()
        end_idx = b.GetEndAtomIdx()
        b_order = b.GetBondTypeAsDouble()
        bond_properties[i] = [begin_idx, end_idx, BOND_TYPES[b_order]]

    return atom_properties, bond_properties


def create_skeleton_mol(connectivity):

    chem_symbols = ["Al", 'Si', 'Se', 'Ca', 'Mg', 'Mn', 'Fe', 'Zn', 'Co', 'Cu', 'Ni', 'Cd', 'Br', 'Cl', 'C', 'H', 'O', 'N', 'P', 'S', 'I', 'F', 'B']

    atoms = list(set([cc for c in connectivity for cc in c[:2]]))
    chem_symbols = [[s for s in chem_symbols if "".join([aa for aa in a if not aa.isdigit()]).capitalize().startswith(s)][0] for a in atoms]

    if len(atoms) != len(chem_symbols):
        raise Exception("Restoring connectivity failed")

    skeleton_mol = Chem.RWMol()

    for a, symbol  in zip(atoms,chem_symbols):
        idx = skeleton_mol.AddAtom(Chem.Atom(symbol))
        skeleton_mol.GetAtomWithIdx(idx).SetProp("_Name", a)

    for begin_idx, end_idx, bond_order in connectivity:
        begin_idx = atoms.index(begin_idx)
        end_idx   = atoms.index(end_idx)
        skeleton_mol.AddBond(begin_idx, end_idx, order=BOND_TYPES[bond_order])

    return skeleton_mol


def parse_params(paramsstring):

    smiles = "None"
    name = None
    connectivity = []

    for s in paramsstring.split("\n"):

        if s.startswith("#SMILES: "):
            smiles = s.split(" ")[1]

        if s.startswith("IO_STRING"):
            name = s.split(" ")[1]

        if s.startswith("BOND"):
            s = list(filter(lambda x: x != "", s.split(" ")))
            connectivity.append([s[1], s[2], int(s[-1].replace("#ORGBND", ""))])

    smiles_params = Chem.SmilesParserParams()
    smiles_params.removeHs = False
    template = Chem.MolFromSmiles(smiles, smiles_params)
    
    if template is None:
        smiles = "None"
    
    skeleton_mol = create_skeleton_mol(connectivity)
    
    if smiles != "None":
        skeleton_mol = AllChem.AssignBondOrdersFromTemplate(template, skeleton_mol)

        if template.HasSubstructMatch(skeleton_mol):
            skeleton_atom_props  = {} 
            skeleton_bonds_props = {}

            match = template.GetSubstructMatch(skeleton_mol)

            for i, j in enumerate(match):
                ref_atom = template.GetAtomWithIdx(j)
                target_atom = skeleton_mol.GetAtomWithIdx(i)

                target_atom.SetFormalCharge(ref_atom.GetFormalCharge())

        elif skeleton_mol.HasSubstructMatch(template):
            skeleton_atom_props  = {} 
            skeleton_bonds_props = {}

            match = skeleton_mol.GetSubstructMatch(template)

            for j, i in enumerate(match):
                ref_atom = template.GetAtomWithIdx(j)
                target_atom = skeleton_mol.GetAtomWithIdx(i)

                target_atom.SetFormalCharge(ref_atom.GetFormalCharge())


    atom_names = [atom.GetProp("_Name") for atom in skeleton_mol.GetAtoms()]
    skeleton_atom_props, skeleton_bonds_props = extract_structure_properties(skeleton_mol)

    skeleton_atom_props = {atom_names[idx]:props for idx, props in skeleton_atom_props.items()}
    skeleton_bonds_props = [[atom_names[props[0]], atom_names[props[1]], props[2]] for _, props in skeleton_bonds_props.items()]

    reshaped_skeleton_bonds_props = {}

    for node in atom_names:
        reshaped_skeleton_bonds_props[node] = []
        
        to_delete = []
        for i, c in enumerate(skeleton_bonds_props):

            if node not in c:
                continue
            
            node2 = c[1] if c.index(node) == 0 else c[0]
            bond_order = c[-1]

            reshaped_skeleton_bonds_props[node].append([node2, bond_order])

            to_delete.append(i)

        for i in sorted(to_delete, reverse=True):
            del skeleton_bonds_props[i]

    reshaped_skeleton_bonds_props = {k:v for k, v in reshaped_skeleton_bonds_props.items() if len(v) != 0}

    return name, reshaped_skeleton_bonds_props, skeleton_atom_props


def load_pdbstring(pdbstring, params=[]):

    mol = Chem.MolFromPDBBlock(pdbstring, removeHs=False, proximityBonding=False, sanitize=False)
    mol_structure = parse_structure(mol)

    params = [parse_params(p) for p in params]
    params = {name:[bonds_props, atom_props] for name, bonds_props, atom_props in params}

    mol = assign_bonds(mol, mol_structure, params)
    # mol.UpdatePropertyCache()

    problems = Chem.DetectChemistryProblems(mol)

    Chem.SanitizeMol(mol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE)
   
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


def read_mol2(mol2fname):

    bonds_dict = {"ar":Chem.BondType.AROMATIC,
                  "am":Chem.BondType.SINGLE,
                  "1" :Chem.BondType.SINGLE,
                  "2" :Chem.BondType.DOUBLE,
                  "3" :Chem.BondType.TRIPLE}


    mol2 = open(mol2fname, "r").read()

    # obabel case
    if len(mol2.split("\n")[1].split(" ")) == 2:
        smiles = mol2.split("\n")[1].split(" ")[-1]
        try:
            mol = Chem.MolFromSmiles(smiles)
            if Chem.MolFromSmiles(smiles) is not None:
                mol = Chem.AddHs(mol)
                return smiles
        except:
            mol = Chem.RWMol()
    else:
        mol = Chem.RWMol()

    # openeye case

    atoms_mol2 = [r for r in mol2.split("@") if r.startswith("<TRIPOS>ATOM\n")][0]
    atoms_mol2 = atoms_mol2.split("<TRIPOS>ATOM\n")[-1]

    atoms_mol2 = [[rr for rr in r.split(" ") if rr != ""] for r in atoms_mol2.split("\n") if len(r) != 0]
    atoms_mol2 = {r[0]:r[5].split(".")[0] for r in atoms_mol2}

    atoms_rdkit = {idx:mol.AddAtom(Chem.Atom(symbol)) for idx, symbol in atoms_mol2.items()}
    atoms_rdkit_reverse = {v:k for k,v in atoms_rdkit.items()}

    bonds = [r for r in mol2.split("@") if r.startswith("<TRIPOS>BOND\n")][0]
    bonds = bonds.split("<TRIPOS>BOND\n")[-1]

    bonds = [[rr for rr in r.split(" ") if rr != ""] for r in bonds.split("\n") if len(r) != 0]
    bonds = [r[1:] for r in bonds]

    for begin_idx, end_idx, b_order in bonds:
        begin_idx = atoms_rdkit[begin_idx]
        end_idx   = atoms_rdkit[end_idx]
        
        mol.AddBond(begin_idx, end_idx, order=bonds_dict[b_order])

    explicit_valence = {"O":2, "N":3, "C":4, "H":1}

    for i, a in enumerate(mol.GetAtoms()):
        a_symbol = a.GetSymbol()
        a.UpdatePropertyCache(strict=False)
        a.SetNoImplicit(True)

        if a_symbol in explicit_valence:
            formal_charge = a.GetExplicitValence() - explicit_valence[a_symbol]
            a.SetFormalCharge(formal_charge)

    if mol.GetNumAtoms() != Chem.AddHs(mol).GetNumAtoms():
        return None

    if Chem.SanitizeMol(mol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE, catchErrors=True) != 0:
        return None
        
    return Chem.MolToSmiles(mol)


def find_polar_groups(mol):

    atom_names = [a.GetPDBResidueInfo().GetName().strip() for a in mol.GetAtoms()]

    templates = {"NH1":"[#7H][H]", "NH2":"[#7H2]([H])[H]", "NH3":"[#7H3]([H])([H])[H]", "OH":"[OH][H]", "carbonyl O":"[#6]([!O])([!O])=[O&!R]", "carboxylate O":"[#6](=[#8&!R])[#8&!R&D1]"}    
    templates = {k:Chem.MolFromSmarts(v) for k, v in templates.items()}
        
    groups = {k:list(mol.GetSubstructMatches(v)) for k, v in templates.items()}

    for i, g in enumerate(groups["NH1"]):
        g = [idx for idx in g if mol.GetAtomWithIdx(idx).GetSymbol() != "N"]
        groups["NH1"][i] = [atom_names[a] for a in g]

    for i, g in enumerate(groups["NH2"]):
        g = [idx for idx in g if mol.GetAtomWithIdx(idx).GetSymbol() != "N"]
        groups["NH2"][i] = [atom_names[a] for a in g]

    for i, g in enumerate(groups["NH3"]):
        g = [idx for idx in g if mol.GetAtomWithIdx(idx).GetSymbol() != "N"]
        groups["NH3"][i] = [atom_names[a] for a in g]
    
    for i, g in enumerate(groups["carbonyl O"]):
        g = [idx for idx in g if mol.GetAtomWithIdx(idx).GetSymbol() == "O"]
        groups["carbonyl O"][i] = [atom_names[a] for a in g]

    for i, g in enumerate(groups["carboxylate O"]):
        g = [idx for idx in g if mol.GetAtomWithIdx(idx).GetSymbol() == "O"]
        groups["carboxylate O"][i] = [atom_names[a] for a in g]

    return groups


def calculate_burunsat_group(burunsat_atoms, polar_groups):


    burunsat_groups = {}

    for grp_name, grp_atoms in polar_groups.items():

        burunsat_groups[grp_name] = []

        for g in grp_atoms:

            atoms = set(burunsat_atoms).intersection(g)

            if len(atoms) == len(g):
                burunsat_groups[grp_name].append(atoms)

        if len(grp_atoms) == 0:
            burunsat_groups[grp_name] = 0
        else:
            burunsat_groups[grp_name] = len(burunsat_groups[grp_name]) / len(grp_atoms)

    return burunsat_groups
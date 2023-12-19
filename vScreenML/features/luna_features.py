from biopandas.pdb import PandasPdb

from luna.mol.entry import MolEntry
from luna.mol.groups import AtomGroupPerceiver
from luna.mol.features import FeatureExtractor

from luna.util.default_values import ATOM_PROP_FILE

from luna.MyBio.util import get_entity_from_entry
from luna.MyBio.PDB.PDBParser import PDBParser

from luna.interaction.calc import InteractionCalculator
from luna.interaction.contact import get_contacts_with
from luna.interaction.filter import InteractionFilter

from vScreenML.utils.rdkit_utils import load_pdbstring

from rdkit.Chem import ChemicalFeatures

def process_pdb(pdb_file):

    pdb = PandasPdb().read_pdb_from_list(pdb_file)
    atoms = pdb.df['ATOM']

    idx = 1
    res_name = None
    resnumbers = []

    for chain, res_number in atoms[["insertion", "residue_number"]].values:

        if res_name is None:
            res_name = str(res_number) + chain


        if res_name == (str(res_number) + chain):
            resnumbers.append(idx)
        else:
            res_name = str(res_number) + chain
            idx += 1
            resnumbers.append(idx)

    atoms["residue_number"] = resnumbers
    pdb.df['ATOM'] = atoms

    pdb = pdb.to_pdb_stream()

    return pdb


def parse_pdb(pdb_file):
    
    pdb_data = {}
    
    current_chain = None
    current_residue = None

    connect = []

    file = process_pdb(pdb_file).read().split("\n")

    for line in file:

        if line.startswith("ATOM") or line.startswith("HETATM"):
            
            chain_id = line[21]
            residue_number = int(line[22:26].strip())
            pdb_string = line.rstrip()

            if chain_id != current_chain:
                current_chain = chain_id
                pdb_data[current_chain] = {}

            if residue_number != current_residue:
                current_residue = int(residue_number)
                pdb_data[current_chain][current_residue] = []

            pdb_data[current_chain][current_residue].append(pdb_string)

        elif line.startswith("CONECT"):
            connect.append(line)

            
    pdb_data = {c:{n:"\n".join(l) for n,l in r.items()} for c, r in pdb_data.items()}
    connect = "".join(connect)

    return pdb_data, connect

def process_structure(pdb_file, params_strings, ligand_resname, ligand_chain, ligand_num):

    pdb_strings, connect = parse_pdb(pdb_file.split("\n"))
    pdb_id = "TEMP"

    # entry of ligand
    entry = MolEntry(pdb_id, ligand_chain, ligand_resname, ligand_num)

    # pdb parser
    pdb_parser = PDBParser(QUIET=True)

    structure = pdb_parser.get_structure(pdb_id, process_pdb(pdb_file.split("\n")))

    ligand = get_entity_from_entry(structure, entry)
    ligand.set_as_target(is_target=True)

    nb_pairs = get_contacts_with(structure[0], ligand, level='R', radius=6.2)
    nb_pairs = [pp for p in nb_pairs for pp in p]
    nb_pairs = set(nb_pairs)

    nb_mol_objs = {}

    for n in nb_pairs:
        _, _, chain, (res_name, res_number, _) = n.get_full_id()
        pdb_string = pdb_strings[chain][res_number]
        
        if res_name == "H_" + ligand_resname:

            pdb_string += "\nTER\n" + connect 

            # template = Chem.MolFromSmiles(pdb_smi)
            # template = Chem.AddHs(template)

            try:
                mol = load_pdbstring(pdb_string, params=params_strings)

                # mol = Chem.MolFromPDBBlock(pdb_string, removeHs=False, proximityBonding=False)
                # rdDetermineBonds.DetermineBonds(mol, charge=0)
                # mol = AllChem.AssignBondOrdersFromTemplate(template, mol)
                nb_mol_objs[n.get_full_id()] = mol
            except:
                try:
                    mol = Chem.MolFromPDBBlock(pdb_string, removeHs=False, proximityBonding=True)
                    mol = AllChem.AssignBondOrdersFromTemplate(template, mol)
                    nb_mol_objs[n.get_full_id()] = mol
                except:
                    print("Failed for:", pdb_file, pdb_smi)
                    return

        else:
            # print(n.get_full_id())
            # print(pdb_string)
            # print()
            nb_mol_objs[n.get_full_id()] = load_pdbstring(pdb_string)

    feature_factory = ChemicalFeatures.BuildFeatureFactory(ATOM_PROP_FILE)
    feature_extractor = FeatureExtractor(feature_factory)

    perceiver = AtomGroupPerceiver(feature_extractor, expand_selection=False)
    atm_grps_mngr = perceiver.perceive_atom_groups(nb_pairs, mol_objs_dict=nb_mol_objs)

    inter_calc = InteractionCalculator(inter_filter=InteractionFilter.new_pli_filter(), add_proximal=False)
    interactions_mngr = inter_calc.calc_interactions(atm_grps_mngr.atm_grps)

    # atm_grps_mngr.merge_hydrophobic_atoms(interactions_mngr)

    interactions = {"Proximal": 0, "Hydrogen bond": 0, "Ionic": 0, "Salt bridge": 0, "Cation-pi": 0, "Hydrophobic": 0, "Halogen bond": 0, "Repulsive": 0, "Water-bridged hydrogen bond": 0, "Amide-aromatic stacking": 0, "Weak hydrogen bond": 0, "Covalent bond": 0, "Atom overlap": 0, "Van der Waals clash": 0, "Van der Waals": 0, "Chalcogen bond": 0, "Chalcogen-pi": 0, "Halogen-pi": 0, "Orthogonal multipolar": 0, "Parallel multipolar": 0, "Antiparallel multipolar": 0, "Tilted multipolar": 0, "Multipolar": 0, "Cation-nucleophile": 0, "Anion-electrophile": 0, "Unfavorable anion-nucleophile": 0, "Unfavorable cation-electrophile": 0, "Unfavorable nucleophile-nucleophile": 0, "Unfavorable electrophile-electrophile": 0, "Pi-stacking": 0, "Face-to-face pi-stacking": 0, "Face-to-edge pi-stacking": 0, "Face-to-slope pi-stacking": 0, "Edge-to-edge pi-stacking": 0, "Edge-to-face pi-stacking": 0, "Edge-to-slope pi-stacking": 0, "Displaced face-to-face pi-stacking": 0, "Displaced face-to-edge pi-stacking": 0, "Displaced face-to-slope pi-stacking": 0}

    for inter in interactions_mngr.interactions:
        grp1 = ";".join(sorted(["/".join(a.full_atom_name.split("/")) for a in inter.src_grp.atoms]))
        grp2 = ";".join(sorted(["/".join(a.full_atom_name.split("/")) for a in inter.trgt_grp.atoms]))

        grp1, grp2 = sorted([grp1, grp2])

        if ligand_resname in grp1 and ligand_resname in grp2:
            continue

        interactions[inter.type] += 1

    return interactions


def calculate_features(pdbstring, params_strings):

    for l in pdbstring.split("\n"):
        if l.startswith("HETATM") and "LG1" in l[17:21]:

            ligand_chain = l[21:22]
            ligand_num = int(l[22:27])
            
            break

    features = process_structure(pdbstring, params_strings, ligand_resname="LG1", ligand_chain="X", ligand_num=1)
    return features

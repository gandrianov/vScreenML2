import importlib

if importlib.util.find_spec("pyrosetta") is None:
    raise ModuleNotFoundError("PyRosetta is not installed")

import pyrosetta
from pyrosetta.rosetta.core.import_pose import pose_from_pdbstring
from pyrosetta.rosetta.protocols.rigid import RigidBodyTransMover
from pyrosetta.rosetta.protocols.grafting import delete_region

from pyrosetta.rosetta.core.scoring import get_score_function
from pyrosetta.rosetta.core.kinematics import MoveMap
from pyrosetta.rosetta.protocols.minimization_packing import MinMover


def init_pyrosetta(args):
    pyrosetta.init(args)

def load_pdbstring(pdbstring):

    pose = pyrosetta.Pose()    
    pose_from_pdbstring(pose, pdbstring)

    return pose


def unbound_pose(pose, residue_id=None):
    
    unbound = pose.clone()

    if residue_id is None:
        
        # if residue is not specified, then pick last jump (chain)
        # in most cases last jump is target ligand (depends how PDB of complex was prepared)
        jump_id = unbound.num_jump()
    
    else:
        # define residue_id belongs to what jumps
        n_residues = unbound.size()
        jump = unbound.fold_tree().get_residue_edge(residue_id)
        jump_id = jump.label() #e.is_jump()            
        
    trans_mover = RigidBodyTransMover(unbound, jump_id)
    trans_mover.step_size(1000000.0)
    trans_mover.apply(unbound)

    return unbound


def get_residue_index(pose, name):

    for i in range(1, pose.size()+1):
        if pose.residue(i).name() == name:
            return i


def decompose_pose(bound, residue_id):

    ligand = bound.clone()
    protein = bound.clone()

    # remove protein from bound
    
    # if ligand residue is first
    if residue_id == 1:
        delete_region(ligand, 2, bound.size())

    # if ligand residue is last
    elif residue_id == bound.size():
        delete_region(ligand, 1, bound.size() - 1)

    # if ligand residue is in the middle
    else:
        delete_region(ligand, residue_id + 1, bound.size())
        delete_region(ligand, 1, residue_id - 1)

    # remove ligand from bound
    delete_region(protein, residue_id, residue_id)


    return ligand, protein


def export_pdbstring(pose):

    """Prepares PDB string"""
    
    buffer = pyrosetta.rosetta.std.stringbuf()
    pose.dump_pdb(pyrosetta.rosetta.std.ostream(buffer))
    
    return buffer.str()


def minimize_complex(pose):

    """Minimizes merged protein-ligand complex"""
    
    # all manipulations were done on cloned pose
    minimized_pose = pose.clone()

    sfxn = get_score_function()
    
    # create move map 
    mmap = MoveMap()
    # allow flexibility on side chains
    mmap.set_chi(True)
    # allow flexibility on back bones
    mmap.set_bb(True)
    # allow flexibility on chains
    mmap.set_jump(True)
    
    # minimizer
    minimizer = MinMover(movemap_in=mmap,
                         scorefxn_in=sfxn,
                         min_type_in="dfpmin_armijo_nonmonotone",
                         tolerance_in=0.000001,
                         use_nb_list_in=True,
                         deriv_check_in=False,
                         deriv_check_verbose_in=False)

    minimizer.max_iter(10000)
    minimizer.apply(minimized_pose)
    
    return minimized_pose
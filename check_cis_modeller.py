import sys
from modeller import *
from modeller.scripts import complete_pdb
from modeller.optimizers import ConjugateGradients, MolecularDynamics
from modeller.automodel import autosched
import math
import numpy as np
import random
import time

# Mapping of single-letter to three-letter residue codes
residue_map = {
    'A': 'ALA', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU', 'F': 'PHE', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
    'K': 'LYS', 'L': 'LEU', 'M': 'MET', 'N': 'ASN', 'P': 'PRO', 'Q': 'GLN', 'R': 'ARG', 'S': 'SER',
    'T': 'THR', 'V': 'VAL', 'W': 'TRP', 'Y': 'TYR'
}

def convert_residues_and_atoms(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                resname = line[17:20].strip()
                if len(resname) == 1:
                    line = line[:17] + residue_map[resname].ljust(3) + line[20:]
                # Set atom name to "CA"
                line = line[:12] + ' CA '.ljust(4) + line[16:]
            outfile.write(line)

def reflect_point_across_line(P, A, B):
    """
    Reflect point P across the line defined by points A and B.
    """
    # Direction vector of the line
    d = B - A
    # Normalize the direction vector
    u = d / np.linalg.norm(d)
    # Vector from A to P
    v = P - A
    # Projection of v onto u
    v_parallel = np.dot(v, u) * u
    # Perpendicular component
    v_perpendicular = v - v_parallel
    # Reflected vector
    v_reflected = v_parallel - v_perpendicular
    # New position of P
    P_reflected = A + v_reflected
    return P_reflected

def correct_cis_peptides(mdl, move_N=True, move_C=True, optimization=True):
    at = mdl.atoms
    residues = mdl.residues

    reslist=[]
    for i in range(1, len(residues)):
        chain_id = residues[i].chain.name
        try:
            # Calculate omega dihedral angle correctly
            atom_CA_prev = at['CA:%d:%s' % (i-1, chain_id)]
            atom_C_prev = at['C:%d:%s' % (i-1, chain_id)]
            atom_N_curr = at['N:%d:%s' % (i, chain_id)]
            atom_CA_curr = at['CA:%d:%s' % (i, chain_id)]
            
            omega = features.Dihedral(atom_CA_prev, atom_C_prev, atom_N_curr, atom_CA_curr)
            omega_value = omega.get_value()

            # Check if omega is in cis configuration
            if math.radians(-89) < omega_value < math.radians(89) or omega_value > math.radians(271) or omega_value < math.radians(-271):
                reslist.append(i-1)
                reslist.append(i)

                # Coordinates of the atoms
                CA_prev_coords = np.array([atom_CA_prev.x, atom_CA_prev.y, atom_CA_prev.z])
                CA_curr_coords = np.array([atom_CA_curr.x, atom_CA_curr.y, atom_CA_curr.z])
               
                if move_N == True and move_C == False:
                    N_curr_coords = np.array([atom_N_curr.x, atom_N_curr.y, atom_N_curr.z])
                    # Reflect the N atom across the line defined by CA_prev and CA_curr
                    N_curr_new_coords = reflect_point_across_line(N_curr_coords, CA_prev_coords, CA_curr_coords)
                    atom_N_curr.x, atom_N_curr.y, atom_N_curr.z = N_curr_new_coords
                elif move_N == False and move_C == True:
                    C_prev_coords = np.array([atom_C_prev.x, atom_C_prev.y, atom_C_prev.z])
                    # Reflect the C atom across the line defined by CA_prev and CA_curr
                    C_prev_new_coords = reflect_point_across_line(C_prev_coords, CA_prev_coords, CA_curr_coords)
                    atom_C_prev.x, atom_C_prev.y, atom_C_prev.z = C_prev_new_coords
                elif move_N == True and move_C == True:
                    # Pick one randomly, but this does not work as good as sticking to the same atom.
                    if random.choice([True, False]):
                        N_curr_coords = np.array([atom_N_curr.x, atom_N_curr.y, atom_N_curr.z])
                        N_curr_new_coords = reflect_point_across_line(N_curr_coords, CA_prev_coords, CA_curr_coords)
                        atom_N_curr.x, atom_N_curr.y, atom_N_curr.z = N_curr_new_coords
                    else:
                        C_prev_coords = np.array([atom_C_prev.x, atom_C_prev.y, atom_C_prev.z])
                        C_prev_new_coords = reflect_point_across_line(C_prev_coords, CA_prev_coords, CA_curr_coords)
                        atom_C_prev.x, atom_C_prev.y, atom_C_prev.z = C_prev_new_coords
                else:
                    continue
                
                print(f"Cis bond at residue {i}. Omega angle is changed from "  + str(omega_value) + " to " + str(omega.get_value()))
        except KeyError:
            # Skip if atoms are not found
            continue

    # Local minimization only for the moved residues
    if optimization==True and len(reslist)>0:
        atsel = Selection([residues[r] for r in reslist])
                    
        cg = ConjugateGradients()
        cg.optimize(atsel, max_iterations=200, min_atom_shift=0.001)

    return len(reslist)
    
def check_cis_peptides(mdl):
    at = mdl.atoms
    residues = mdl.residues

    reslist=[]

    for i in range(1, len(residues)):
        chain_id = residues[i].chain.name
        
        atom_CA_prev = at['CA:%d:%s' % (i, chain_id)]
        atom_C_prev = at['C:%d:%s' % (i, chain_id)]
        atom_N_curr = at['N:%d:%s' % (i+1, chain_id)]
        atom_CA_curr = at['CA:%d:%s' % (i+1, chain_id)]
        
        omega = features.Dihedral(atom_CA_prev, atom_C_prev, atom_N_curr, atom_CA_curr)
        omega_value = omega.get_value()

        print(f"Residue {i} -- {i+1}. Omega angle is "  +  str(math.degrees(omega_value)), end=" " )

        # Check if omega is in cis configuration
        if math.radians(-89) < omega_value < math.radians(89) or omega_value > math.radians(271) or omega_value < math.radians(-271):
            reslist.append(i)
            reslist.append(i+1)         
            print(f"CIS detected.", end=" " )

        print( " " )

    return len(reslist)/2

def optimize(atmsel, sched):
    #conjugate gradient
    for step in sched:
        step.optimize(atmsel, max_iterations=200, min_atom_shift=0.001)
    #md
    refine(atmsel)
    cg = ConjugateGradients()
    cg.optimize(atmsel, max_iterations=200, min_atom_shift=0.001)


#molecular dynamics
def refine(atmsel):
    # at T=1000, max_atom_shift for 4fs is cca 0.15 A.
    md = MolecularDynamics(cap_atom_shift=0.39, md_time_step=4.0,
                           md_return='FINAL')
    init_vel = True
    for (its, equil, temps) in ((200, 20, (150.0, 250.0, 400.0, 700.0, 1000.0)),
                                (200, 600,
                                 (1000.0, 800.0, 600.0, 500.0, 400.0, 300.0))):
        for temp in temps:
            md.optimize(atmsel, init_velocities=init_vel, temperature=temp,
                         max_iterations=its, equilibrate=equil)
            init_vel = False


if len(sys.argv) != 2:
    print("Usage: python complete_structure_modeller.py input_pdb")
    sys.exit(1)

input_pdb = sys.argv[1]

log.verbose()
env = environ()


# Read topology and parameter libraries
env.libs.topology.read(file='$(LIB)/top_heav.lib')
env.libs.parameters.read(file='$(LIB)/par.lib')

# Create a model object
mdl = Model(env)

# Read in the PDB file
mdl.read(file=input_pdb)

# Function to list all atoms
def list_all_atoms(model):
    for chain in model.chains:
        for residue in chain.residues:
            for atom in residue.atoms:
                print(f'Chain: {chain.name}, Residue: {residue.num} {residue.name}, Atom: {atom.name}, Type: {atom.type}, Coordinates: ({atom.x}, {atom.y}, {atom.z})')

# List all atoms in the model
#list_all_atoms(mdl)

###
print("Checking CIS")
ncis = check_cis_peptides(mdl)
print("Number of CIS = ", str(ncis))

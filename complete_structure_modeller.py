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
        r_prev = residues[i-1].num
        r_curr = residues[i].num
        #print(r_prev, r_curr, i-1, i) 
        
        # Calculate omega dihedral angle correctly
        atom_CA_prev = at['CA:%s:%s' % (r_prev, chain_id)]
        atom_C_prev = at['C:%s:%s' % (r_prev, chain_id)]
        atom_N_curr = at['N:%s:%s' % (r_curr, chain_id)]
        atom_CA_curr = at['CA:%s:%s' % (r_curr, chain_id)]
        
        omega = features.Dihedral(atom_CA_prev, atom_C_prev, atom_N_curr, atom_CA_curr)
        omega_value = omega.get_value()

        # Check if omega is in cis configuration
        if math.radians(-89) < omega_value < math.radians(89) or omega_value > math.radians(271) or omega_value < math.radians(-271):
            # Coordinates of the atoms
            CA_prev_coords = np.array([atom_CA_prev.x, atom_CA_prev.y, atom_CA_prev.z])
            CA_curr_coords = np.array([atom_CA_curr.x, atom_CA_curr.y, atom_CA_curr.z])
           
            if move_N == True and move_C == False:
                N_curr_coords = np.array([atom_N_curr.x, atom_N_curr.y, atom_N_curr.z])
                # Reflect the N atom across the line defined by CA_prev and CA_curr
                N_curr_new_coords = reflect_point_across_line(N_curr_coords, CA_prev_coords, CA_curr_coords)
                atom_N_curr.x, atom_N_curr.y, atom_N_curr.z = N_curr_new_coords
                reslist.append(i-1)
            elif move_N == False and move_C == True:
                C_prev_coords = np.array([atom_C_prev.x, atom_C_prev.y, atom_C_prev.z])
                # Reflect the C atom across the line defined by CA_prev and CA_curr
                C_prev_new_coords = reflect_point_across_line(C_prev_coords, CA_prev_coords, CA_curr_coords)
                atom_C_prev.x, atom_C_prev.y, atom_C_prev.z = C_prev_new_coords
                reslist.append(i)
            elif move_N == True and move_C == True:
                # Pick one randomly, but this does not work as good as sticking to the same atom.
                if random.choice([True, False]):
                    N_curr_coords = np.array([atom_N_curr.x, atom_N_curr.y, atom_N_curr.z])
                    N_curr_new_coords = reflect_point_across_line(N_curr_coords, CA_prev_coords, CA_curr_coords)
                    atom_N_curr.x, atom_N_curr.y, atom_N_curr.z = N_curr_new_coords
                    reslist.append(i-1)
                else:
                    C_prev_coords = np.array([atom_C_prev.x, atom_C_prev.y, atom_C_prev.z])
                    C_prev_new_coords = reflect_point_across_line(C_prev_coords, CA_prev_coords, CA_curr_coords)
                    atom_C_prev.x, atom_C_prev.y, atom_C_prev.z = C_prev_new_coords
                    reslist.append(i)
            else:
                continue
            
            print(f"Cis bond at residue {r_prev} -- {r_curr}. Omega angle is changed from "  + str(omega_value) + " to " + str(omega.get_value()))

    # Local minimization only for the moved residues
    if optimization==True and len(reslist)>0:
        atsel = Selection([residues[r] for r in reslist])
                    
        cg = ConjugateGradients()
        cg.optimize(atsel, max_iterations=200, min_atom_shift=0.001)

    return len(reslist)
    

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


if len(sys.argv) != 3:
    print("Usage: python complete_structure_modeller.py input_pdb output_pdb")
    sys.exit(1)

input_pdb = sys.argv[1]
intermediate_pdb = 'intermediate.pdb'
output_pdb = sys.argv[2]

# Convert single-letter residues to three-letter residues and set atom names to "CA"
convert_residues_and_atoms(input_pdb, intermediate_pdb)

log.verbose()
env = environ()
# Set the random seed using the current wall clock time
seed = int(time.time() % 10000)  # Modeller's random seed should be an integer
env.rand_seed = seed
print(f"Random seed set to: {seed}")

# Read topology and parameter libraries
env.libs.topology.read(file='$(LIB)/top_heav.lib')
env.libs.parameters.read(file='$(LIB)/par.lib')

# Add missing atoms
mdl = complete_pdb(env, intermediate_pdb)

# Set up restraints to ensure stereochemical quality
mdl.restraints.make(mdl, restraint_type='stereo', spline_on_site=False)

# Add harmonic positional restraints to keep Cα atoms close to their original positions
rsr = mdl.restraints
at = mdl.atoms

# Apply positional restraints to all Cα atoms
for atom in at:
    if atom.name == 'CA':
        rsr.add(forms.Gaussian(group=physical.xy_distance,
                               feature=features.XCoordinate(atom),
                               mean=atom.x, stdev=1))
        rsr.add(forms.Gaussian(group=physical.xy_distance,
                               feature=features.YCoordinate(atom),
                               mean=atom.y, stdev=1))
        rsr.add(forms.Gaussian(group=physical.xy_distance,
                               feature=features.ZCoordinate(atom),
                               mean=atom.z, stdev=1))

# Apply trans bond restraints on the omega angle
residues = mdl.residues
for i in range(1, len(mdl.residues)):
    chain_id = residues[i].chain.name
    r_prev = residues[i-1].num
    r_curr = residues[i].num

    omega = features.Dihedral(at['CA:%s:%s' % (r_prev, chain_id)], at['C:%s:%s' % (r_prev, chain_id)], at['N:%s:%s' % (r_curr, chain_id)], at['CA:%s:%s' % (r_curr, chain_id)])
    # Enforce trans configuration for all residues
    rsr.add(forms.Gaussian(group=physical.phi_psi_dihedral,
                           feature=omega, mean=math.pi, stdev=1))


md = MolecularDynamics(output='REPORT')
cg = ConjugateGradients()

# initial minimization
cg.optimize(mdl, max_iterations=200, min_atom_shift=0.001)

#correct the cis
ncis = correct_cis_peptides(mdl, move_N=False, move_C=True, optimization=True)

# initial optimization
#sched = autosched.loop.make_for_model(mdl)
#optimize(Selection(mdl), sched)
md.optimize(mdl, temperature=300, max_iterations=1000, equilibrate=50)
cg.optimize(mdl, max_iterations=200, min_atom_shift=0.001)


# below we do optimization until no cis bonds are found but we do more to make sure
# Initialize counter for additional iterations
additional_iterations = 0
max_additional_iterations = 20

ncis = correct_cis_peptides(mdl, move_N=False, move_C=True, optimization=True)
while ncis > 0 or additional_iterations < max_additional_iterations:
    md.optimize(mdl, temperature=300, max_iterations=100, equilibrate=50)
    cg.optimize(mdl, max_iterations=200, min_atom_shift=0.001)
    ncis = correct_cis_peptides(mdl, move_N=False, move_C=True, optimization=True)
    
    if ncis == 0:
        additional_iterations += 1
    else:
        additional_iterations = 0  # Reset if more cis bonds are found


###
print("Checking CIS")
ncis = correct_cis_peptides(mdl,  move_N=False, move_C=True, optimization=False)
print("Number of CIS = ", str(ncis))
####

# Remove all restraints except stereochemical ones and minimize
mdl.restraints.clear() 
mdl.restraints.make(mdl, restraint_type='stereo', spline_on_site=False)

cg.optimize(mdl, max_iterations=200, min_atom_shift=0.001)

###
print("Checking CIS again")
ncis = correct_cis_peptides(mdl, move_N=False, move_C=False, optimization=False)
print("Number of CIS = ", str(ncis))
####


# Write out the complete PDB structure
mdl.write(file=output_pdb)

print(f"Complete structure has been written to '{output_pdb}'")

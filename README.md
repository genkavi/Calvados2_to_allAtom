# Calvados2_to_allAtom

This repository contains tools for converting **Calvados2 simulation** for **Intrinsically Disordered Regions (IDRs)** outputs into all-atom structures. 
## Included Scripts

### `complete_structure_modeller.py`

This script is designed to convert **Calvados2** coarse-grained simulation outputs into refined all-atom structures. It adds missing atoms, corrects cis-peptides, and refines the structure to ensure stereochemical accuracy and structural stability.

#### Key Features
- Converts coarse-grained amino acid representations (single-bead) into detailed all-atom structures.
- Detects cis-peptides by analyzing omega dihedral angles.
- Corrects cis-peptides by reflecting atomic coordinates of affected residues.
- Adds missing heavy atoms to ensure the structural integrity of the all-atom model.
- Applies stereochemical restraints and ensures proper geometry using molecular dynamics optimization and conjugate gradient minimization.
- Iteratively optimizes the structure until no cis-peptides remain.

## Usage

```bash
python complete_structure_modeller.py input.pdb output.pdb
```

---

### `check_cis_modeller.py`

This script is designed to identify and correct cis-peptides in protein structures, particularly in **Intrinsically Disordered Regions (IDRs)**. It analyzes the omega dihedral angles between residues and applies corrections to cis-peptides to ensure proper stereochemistry.

#### Features

- Detects cis-peptides by analyzing omega dihedral angles.
- Corrects cis-peptides by adjusting atomic coordinates using reflection across defined axes.
- Supports local optimization of corrected residues using `ConjugateGradients` for stereochemical accuracy.
- Configurable to move the nitrogen atom, carbon atom, or both for corrections.

#### Usage

```bash
python check_cis_modeller.py input.pdb
```

---

## Dependencies

The following dependencies are required to run `complete_structure_modeller.py`:

- Python (3.x)
- Modeller
- Numpy


---

### Use Case: Processing a Trajectory with VMD and `complete_structure_modeller.py`

This example demonstrates how to process a Calvados2 trajectory file (`traj.dcd`) and convert selected frames from a coarse-grained simulation into all-atom representations using `complete_structure_modeller.py`. The workflow uses VMD to extract frames and rename residues, followed by Modeller for structure refinement.

#### Workflow:

1. **Extract and Process Frames with VMD**  
   The following script selects residues, renames them, and saves specific frames as PDB files:
   ```bash
   vmd -dispdev text test/top.pdb test/traj.dcd <<EOF

   \# Calvados2 labels the first and the last residues as X and Z.  
   \# We can either remove them or assign them their correct residue names. 

   set a [ atomselect top "resname X"]
   \$a set resname R

   set a [ atomselect top "resname Z"]
   \$a set resname Q

   set all [ atomselect top all ]

   foreach frame "500 600 700 800 900 1000" {
       \$all frame \$frame
       \$all writepdb \$frame.pdb
   }

   exit
   EOF

   for f in 500 600 700 800 900 1000
   do 
        python complete_structure_modeller.py \$f.pdb \$f.aa.pdb
   done

```
---


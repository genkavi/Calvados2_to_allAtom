# Calvados2_to_allAtom

This repository contains tools for converting **Calvados2 simulation** for **Intrinsically Disordered Regions (IDRs)** outputs into all-atom structures. 
## Included Scripts

### `check_cis_modeller.py`
- **Purpose**: Identifies and corrects cis-peptides in protein structures, particularly in IDRs.
- **Features**:
  - Detects cis-peptides by analyzing omega dihedral angles.
  - Corrects cis-peptides by adjusting atomic coordinates and reflecting relevant atoms.
  - Optimizes modified residues to maintain stereochemical integrity.
- **Usage**:
  ```bash
  python check_cis_modeller.py input.pdb
```

### complete_structure_modeller.py

This script is designed to convert **Calvados2** coarse-grained simulation outputs into refined all-atom structures. It adds missing atoms, corrects cis-peptides, and refines the structure to ensure stereochemical accuracy and structural stability.
---

#### Key Features

1. **Conversion from Single-Bead to All-Atom Structures**:
   - Converts coarse-grained amino acid representations (single-bead) into detailed all-atom structures.

2. **Correction of Cis-Peptides**:
   - Detects cis-peptides by analyzing omega dihedral angles.
   - Corrects cis-peptides by reflecting atomic coordinates of affected residues.

3. **Addition of Missing Atoms**:
   - Adds missing heavy atoms to ensure the structural integrity of the all-atom model.

4. **Structure Refinement**:
   - Applies stereochemical restraints and ensures proper geometry using molecular dynamics optimization and conjugate gradient minimization.

5. **Iterative Correction**:
   - Iteratively optimizes the structure until no cis-peptides remain.

## Usage

```bash
python complete_structure_modeller.py input.pdb output.pdb
```

---

### check_cis_modeller.py

This script is designed to identify and correct cis-peptides in protein structures, particularly in **Intrinsically Disordered Regions (IDRs)**. It analyzes the omega dihedral angles between residues and applies corrections to cis-peptides to ensure proper stereochemistry.

#### Features

- Detects cis-peptides by analyzing omega dihedral angles.
- Corrects cis-peptides by adjusting atomic coordinates using reflection across defined axes.
- Supports local optimization of corrected residues using `ConjugateGradients` for stereochemical accuracy.
- Configurable to move the nitrogen atom, carbon atom, or both for corrections.

#### Usage

```bash
python check_cis_modeller.py input.pdb

## Dependencies

The following dependencies are required to run `complete_structure_modeller.py`:

- Python (3.x)
- Modeller
- Numpy





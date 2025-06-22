# RDF Analysis

This tool processes output files from various simulation software,
calculates the radial distribution function (RDF) to characterize structural properties and interactions,
and generates analysis-ready data files (.dat) along with corresponding RDF plots (.png).
---

## Features

- Supports multiple molecular simulation output formats  
- Calculates radial distribution functions (RDF) accurately  
- Outputs both raw data files (.dat) and corresponding figures 

---

## Installation

It is recommended to create a dedicated Conda environment:

```bash
conda create --name rdf python=3.10
conda activate rdf
pip install -r requirements.txt
```

## Usage
Run the program with your topology and trajectory files as input:

```bash
python main.py <topology_file> <trajectory_file>
```

## input
- VASP: (POSCAR + XDATCAR)
- LAMMPS: (.data + dump trajectory)
- GROMACS: (.tpr + .xtc/.trr)

## Project Structure
- parsers/: Modules to read input files for POSCAR/XDATCAR (VASP), LAMMPS .data, and GROMACS .tpr + trajectory files
- rdf/: Functions to compute and plot RDF
- utils/: Helper functions for user interaction
- outputs/: Directory where analysis results (plots and data) are saved
- data/: (optional) Sample/test input data

## Note
- Interactive mode allows selection of atom pairs for RDF calculation, or calculate all pairs at once
- Output files are automatically named with prefixes according to input type (vasp, lmp, gmx)
- Ensure all dependencies in requirements.txt are installed before running
from pymatgen.io.vasp import Xdatcar, Poscar
from ase import Atoms
import numpy as np
from ase.io import write
import MDAnalysis as mda
from pathlib import Path

#去除Selective Dynamics 标签设置的warning
import warnings
from pymatgen.io.vasp.inputs import BadPoscarWarning
warnings.filterwarnings("ignore", category=BadPoscarWarning)

def read_poscar_box(poscar_path: Path):
    poscar = Poscar.from_file(str(poscar_path))
    lattice = poscar.structure.lattice
    a, b, c = lattice.abc
    alpha, beta, gamma = lattice.angles
    return [a, b, c, alpha, beta, gamma]


def process_vasp(xdatcar_path: Path, poscar_path: Path):
    output_dir = Path("outputs")
    output_dir.mkdir(parents=True, exist_ok=True)

    xdat = Xdatcar(str(xdatcar_path))
    ase_atoms = []
    atom_types = []
    for s in xdat.structures:
        symbols = [site.specie.symbol for site in s]
        positions = s.cart_coords
        cell = s.lattice.matrix
        atoms = Atoms(symbols=symbols, positions=positions, cell=cell, pbc=True)
        ase_atoms.append(atoms)
        atom_types.extend(symbols)

    traj_xyz = output_dir / "vasp_traj.xyz"
    write(traj_xyz, ase_atoms)

    box = read_poscar_box(poscar_path)
    u = mda.Universe(str(traj_xyz))
    for ts in u.trajectory:
        ts.dimensions = np.array(box)

    return u, list(set(atom_types))

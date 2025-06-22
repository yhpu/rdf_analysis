import MDAnalysis as mda
from pathlib import Path

def process_gmx(traj_path: Path, tpr_path: Path):
    # 检查轨迹文件扩展名
    if traj_path.suffix not in [".xtc", ".trr"]:
        raise ValueError(f"Unsupported GROMACS trajectory file: {traj_path.name}. Expected .xtc or .trr")

    try:
        u = mda.Universe(str(tpr_path), str(traj_path))
    except Exception as e:
        raise ValueError(f"Failed to load GROMACS files: {e}")

    print(f"Loaded system with {len(u.atoms)} atoms and {len(u.trajectory)} frames")
    atom_types = sorted(set(u.atoms.types))

    return u, atom_types

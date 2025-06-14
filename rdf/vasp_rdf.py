import numpy as np
from pathlib import Path
from typing import List, Tuple, Dict
from .utils import save_rdf_results
from .utils import save_all_rdf_plot

def read_xdatcar(filename: Path) -> Tuple[List[np.ndarray], np.ndarray, List[str], List[int]]:
    """读取XDATCAR文件"""
    with open(filename, 'r') as f:
        lines = [line.strip() for line in f if line.strip()]
    
    scale = float(lines[1])
    lattice = np.array([list(map(float, lines[i].split())) for i in range(2, 5)]) * scale
    elements = lines[5].split()
    natoms = list(map(int, lines[6].split()))
    total_atoms = sum(natoms)
    
    traj = []
    current_line = 7
    while current_line < len(lines):
        if "Direct configuration" in lines[current_line]:
            current_line += 1
            frame = [list(map(float, lines[current_line + i].split())) for i in range(total_atoms)]
            traj.append(np.array(frame))
            current_line += total_atoms
        else:
            current_line += 1
    return traj, lattice, elements, natoms

def calculate_rdf(traj: List[np.ndarray], 
                  lattice: np.ndarray,
                  elements: List[str],
                  natoms: List[int],
                  pair: Tuple[str, str],
                  r_max: float = 10.0,
                  dr: float = 0.05) -> Tuple[np.ndarray, np.ndarray]:
    """计算RDF"""
    from collections import defaultdict
    from scipy.spatial.distance import cdist
    
    elem_indices = defaultdict(list)
    start = 0
    for elem, n in zip(elements, natoms):
        elem_indices[elem].extend(range(start, start + n))
        start += n
    
    elem1, elem2 = pair
    indices1 = elem_indices[elem1]
    indices2 = elem_indices[elem2] if elem1 != elem2 else elem_indices[elem1]
    
    bins = np.arange(0, r_max + dr, dr)
    rdf = np.zeros(len(bins) - 1)
    n_frames = len(traj)
    volume = np.abs(np.linalg.det(lattice))
    rho = len(indices2) / volume
    
    for frame in traj:
        cart_coords = np.dot(frame, lattice)
        coords1 = cart_coords[indices1]
        coords2 = cart_coords[indices2]
        
        dists = cdist(coords1, coords2)
        if elem1 == elem2:
            dists = dists[dists > 1e-8]
        
        hist, _ = np.histogram(dists[dists < r_max], bins=bins)
        rdf += hist

    r = 0.5 * (bins[1:] + bins[:-1])
    shell_vol = 4 * np.pi * r**2 * dr
    rdf /= (n_frames * len(indices1) * rho * shell_vol)
    
    return r, rdf

def process(xdatcar_path: Path, 
            output_dir: Path = Path("."), 
            r_max: float = 10.0,
            dr: float = 0.05) -> None:
    """处理XDATCAR文件"""
    traj, lattice, elements, natoms = read_xdatcar(xdatcar_path)
    print(f"Loaded {len(traj)} frames with elements: {dict(zip(elements, natoms))}")
    
    unique_elements = list(set(elements))
    all_rdfs = {}
    prefix = "vasp"  # <-- 添加前缀控制
    
    for elem1 in unique_elements:
        for elem2 in unique_elements:
            if elements.index(elem1) > elements.index(elem2):
                continue  # 避免重复计算
            
            pair = (elem1, elem2)
            r, g_r = calculate_rdf(traj, lattice, elements, natoms, pair, r_max, dr)
            pair_name = f"{elem1}-{elem2}"
            save_rdf_results(r, g_r, pair_name, output_dir, prefix=prefix)
            all_rdfs[pair_name] = (r, g_r)

    # 绘制所有 RDF 总图
    save_all_rdf_plot(all_rdfs, output_dir, prefix=prefix)

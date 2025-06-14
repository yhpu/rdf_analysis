import argparse
import MDAnalysis as mda
from MDAnalysis.analysis.rdf import InterRDF
from pathlib import Path
from typing import Optional
from .utils import save_rdf_results

def process(dump_path: Path, 
           output_dir: Path = Path("."),
           r_max: float = 10.0,
           nbins: int = 200) -> None:
    """处理LAMMPS dump文件"""
    try:
        u = mda.Universe(str(dump_path), format="LAMMPSDUMP")
    except Exception as e:
        raise ValueError(f"Failed to load LAMMPS dump: {e}")
    
    print(f"Loaded system with {len(u.atoms)} atoms")
    atom_types = sorted(set(u.atoms.types))
    print("Available atom types:", atom_types)
    
    all_rdfs = {}  # 存储所有pair结果
    
    for i, type1 in enumerate(atom_types):
        for type2 in atom_types[i:]:  # 避免重复计算
            group1 = u.select_atoms(f"type {type1}")
            group2 = u.select_atoms(f"type {type2}")
            
            rdf = InterRDF(group1, group2, range=(0.0, r_max), nbins=nbins)
            rdf.run()
            
            pair_name = f"type{type1}-type{type2}"
            save_rdf_results(
                rdf.results.bins,
                rdf.results.rdf,
                pair_name,
                output_dir,
                prefix="lmp"
            )
            all_rdfs[pair_name] = (rdf.results.bins, rdf.results.rdf)
    
    # 画一张所有rdf叠加的总图，传入前缀方便命名
    from .utils import save_all_rdf_plot
    save_all_rdf_plot(all_rdfs, output_dir, prefix="lmp")
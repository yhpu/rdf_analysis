import MDAnalysis as mda
from MDAnalysis.analysis.rdf import InterRDF
from pathlib import Path
from typing import Optional
from .utils import save_rdf_results, save_all_rdf_plot

def process(xtc_path: Path, 
           tpr_path: Path,
           output_dir: Path = Path("."),
           r_max: float = 10.0,
           nbins: int = 100) -> None:
    """处理GROMACS轨迹"""
    try:
        u = mda.Universe(str(tpr_path), str(xtc_path))
    except Exception as e:
        raise ValueError(f"Failed to load GROMACS files: {e}")
    
    print(f"Loaded system with {len(u.atoms)} atoms and {len(u.trajectory)} frames")
    atom_types = sorted(set(u.atoms.types))
    print("Available atom types:", atom_types)
    
    prefix = "gmx"  # 文件名前缀
    all_rdfs = {}
    # 计算所有类型组合的RDF
    for i, type1 in enumerate(atom_types):
        for type2 in atom_types[i:]:  # 避免重复计算
            group1 = u.select_atoms(f"type {type1}")
            group2 = u.select_atoms(f"type {type2}")
            
            rdf = InterRDF(group1, group2, 
                          range=(0.0, r_max), 
                          nbins=nbins,
                          exclusion_block=(1, 1))  # 排除相邻原子
            rdf.run()
            
            pair_name = f"{type1}-{type2}"
            save_rdf_results(
                rdf.results.bins,
                rdf.results.rdf,
                pair_name,
                output_dir,
                prefix=prefix  # 传入前缀
            )
            all_rdfs[pair_name] = (rdf.results.bins, rdf.results.rdf)
    # 绘制所有RDF总图
    save_all_rdf_plot(all_rdfs, output_dir, prefix=prefix)          
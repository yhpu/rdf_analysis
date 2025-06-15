import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from typing import Tuple, Dict

def save_rdf_results(r: np.ndarray, 
                    g_r: np.ndarray, 
                    pair_name: str, 
                    output_dir: Path, 
                    prefix: str = "rdf") -> None:
    """保存RDF结果到文件"""
    # 保存数据
    data_file = output_dir / f"{prefix}_{pair_name}.dat"
    np.savetxt(data_file, np.column_stack((r, g_r)), 
               header=f"# {pair_name} RDF\n# Distance(A) g(r)")
    
    # 保存图像
    plot_file = output_dir / f"{prefix}_{pair_name}.png"
    plt.figure(figsize=(10, 6))
    plt.plot(r, g_r, linewidth=2)
    plt.xlabel("r (Å)")
    plt.ylabel("g(r)")
    plt.xlim(0.2, 10)
    plt.ylim(-0.5, 80)
    plt.title(f"RDF for {pair_name}")
    plt.grid(alpha=0.3)
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    plt.close()

def save_all_rdf_plot(all_rdfs: Dict[str, Tuple[np.ndarray, np.ndarray]], 
                      output_dir: Path,
                      prefix: str = "rdf") -> None:
    """将所有 RDF 曲线画在一张图里"""
    plt.figure(figsize=(10, 6))
    
    for pair_name, (r, g_r) in all_rdfs.items():
        plt.plot(r, g_r, label=pair_name)
    
    plt.xlabel("r (Å)")
    plt.ylabel("g(r)")
    plt.xlim(0.2, 10)
    plt.ylim(bottom=0)
    plt.title("Radial Distribution Function")
    plt.grid(alpha=0.3)
    plt.legend()
    
    combined_plot_file = output_dir /  f"{prefix}_all.png"
    plt.savefig(combined_plot_file, dpi=300, bbox_inches='tight')
    plt.close()
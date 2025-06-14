import sys
from pathlib import Path
from rdf import vasp_rdf, lmp_rdf, gmx_rdf

def detect_format(file_path: Path) -> str:
    """检测轨迹文件格式"""
    if "XDATCAR" in file_path.name:
        return "vasp"
    elif file_path.suffix in (".dump", ".lammpstrj"):
        return "lmp"
    elif file_path.suffix in (".xtc", ".tpr"):
        return "gmx"
    else:
        raise ValueError(f"Unsupported file format: {file_path}")

def main(file_path: str):
    """主处理函数"""
    file_path = Path(file_path)
    if not file_path.exists():
        raise FileNotFoundError(f"File not found: {file_path}")
    
    fmt = detect_format(file_path)
    output_dir = Path("outputs")
    output_dir.mkdir(exist_ok=True)
    
    try:
        if fmt == "vasp":
            vasp_rdf.process(file_path, output_dir=output_dir)
        elif fmt == "lmp":
            lmp_rdf.process(file_path, output_dir=output_dir) 
        elif fmt == "gmx":
            tpr_file = file_path.with_suffix(".tpr")
            if not tpr_file.exists():
                raise FileNotFoundError(f"Required TPR file not found: {tpr_file}")
            gmx_rdf.process(file_path, tpr_file, output_dir=output_dir)
        
        print(f"\nAnalysis completed. Results saved to {output_dir}/")
    except Exception as e:
        print(f"\nError processing file: {e}")
        sys.exit(1)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python main.py <trajectory_file>")
        sys.exit(1)
    main(sys.argv[1])
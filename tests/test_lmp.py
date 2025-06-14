from pathlib import Path
from rdf import lmp_rdf

def test_lmp_rdf_runs(lmp_sample_path, output_dir):
    """测试LAMMPS RDF处理能否正常运行"""
    lmp_rdf.process(lmp_sample_path, output_dir=output_dir)
    
    output_files = list(output_dir.glob("*.png")) + list(output_dir.glob("*.dat"))
    assert len(output_files) > 0
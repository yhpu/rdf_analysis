from pathlib import Path
from rdf import gmx_rdf

def test_gmx_rdf_runs(gmx_sample_path, gmx_topology_path, output_dir):
    """测试GROMACS RDF处理能否正常运行"""
    gmx_rdf.process(gmx_sample_path, gmx_topology_path, output_dir=output_dir)
    
    output_files = list(output_dir.glob("*.png")) + list(output_dir.glob("*.dat"))
    assert len(output_files) > 0
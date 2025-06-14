from pathlib import Path
from rdf import vasp_rdf

def test_vasp_rdf_runs(vasp_sample_path, output_dir):
    """测试VASP RDF处理能否正常运行"""
    vasp_rdf.process(vasp_sample_path, output_dir=output_dir)
    
    # 检查输出文件
    output_files = list(output_dir.glob("*.png")) + list(output_dir.glob("*.dat"))
    assert len(output_files) > 0
    for f in output_files:
        assert f.stat().st_size > 0
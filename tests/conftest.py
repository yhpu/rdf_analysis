import pytest
from pathlib import Path
import shutil

@pytest.fixture
def output_dir(tmp_path):
    """为每个测试创建独立的输出目录"""
    dir_path = tmp_path / "outputs"
    dir_path.mkdir()
    return dir_path

@pytest.fixture
def vasp_sample_path():
    return Path("data/vasp/XDATCAR")

@pytest.fixture
def lmp_sample_path():
    return Path("data/lmp/dump_Al.lammpstrj")

@pytest.fixture
def gmx_sample_path():
    return Path("data/gmx/md_h2o.xtc")

@pytest.fixture
def gmx_topology_path():
    return Path("data/gmx/md_h2o.tpr")
import MDAnalysis as mda
from pymatgen.core.periodic_table import Element
import warnings

warnings.filterwarnings("ignore", message="Reader has no dt information, set to 1.0 ps")

# 检查 dump 文件字段
def check_dump_columns(traj_path):
    required_fields = {"id", "type", "x", "y", "z"}  # 用集合，不关心顺序
    with open(traj_path, 'r') as f:
        for line in f:
            if line.startswith("ITEM: ATOMS"):
                actual_fields = line.strip().split()[2:]
                break
        else:
            raise ValueError("未找到 'ITEM: ATOMS' 行，无法检查字段")

    actual_fields_set = set(actual_fields)
    if not required_fields.issubset(actual_fields_set):
        raise ValueError(f"dump 文件字段缺失，必须包含 {required_fields}，实际为 {actual_fields}")
    
#读取.data文件中的原子质量
def parse_type_mass_map(data_file):
    type_mass_map = {}
    with open(data_file) as f:
        lines = f.readlines()

    in_masses = False
    for line in lines:
        stripped = line.strip()
        if stripped.lower() == "masses":
            in_masses = True
            continue
        if in_masses:
            if stripped.lower().startswith("atoms"):
                break
            if stripped == "":
                continue
            parts = stripped.split()
            if len(parts) >= 2:
                try:
                    type_id = int(parts[0])
                    mass = float(parts[1])
                    type_mass_map[type_id] = mass
                except:
                    pass
    return type_mass_map

#将原子质量映射到元素类型
def map_mass_to_element_type(type_mass_map, tol=0.1):
    mapped = {}
    for t, mass in type_mass_map.items():
        matched = None
        for el in Element:
            if abs(el.atomic_mass - mass) < tol:
                matched = el
                break
        mapped[t] = matched
    return mapped


def process_lmp(data_path, traj_path):
    check_dump_columns(traj_path)# 检查 ump文件字段是不是符合要求
    type_mass_map = parse_type_mass_map(data_path) #解析原子质量
    type_info = map_mass_to_element_type(type_mass_map)#映射到元素类型

    print("Detected atom types:")
    symbol_map = {}
    for t, el in type_info.items():
        if el is None:
            print(f"Type {t}: Unknown mass = {type_mass_map[t]:.3f}")
            symbol_map[t] = f"type{t}"
        else:
            print(f"Type {t}: {el.symbol} (Z={el.Z}, Mass={el.atomic_mass:.3f})")
            symbol_map[t] = el.symbol

    u = mda.Universe(data_path, traj_path, format='LAMMPSDUMP', atom_style='id type x y z')
    for atom in u.atoms:
        atom.mass = type_mass_map.get(int(atom.type), 0)
    atom_types = sorted(symbol_map.keys())
    return u, atom_types, symbol_map
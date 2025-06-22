import argparse
from pathlib import Path
from parsers import vasp, lmp, gmx
from rdf import rdf_calc
from utils.io import get_user_selection, yes_or_no
import sys

def main():
    parser = argparse.ArgumentParser(description="RDF analysis for VASP, LAMMPS, and GROMACS outputs")
    #必须输入且只能输入两个文件
    parser.add_argument(
        "inputs",
        nargs=2,
        help="Input file paths: 2 for VASP (POSCAR XDATCAR), 2 for LAMMPS (.data .lammpstrj), 2 for GROMACS (.tpr .xtc/.trr)"
    )
    args = parser.parse_args()
    input1, input2 = Path(args.inputs[0]), Path(args.inputs[1])

    prefix = "rdf"  # 默认文件名前缀
    symbol_map = None  # 用于LAMMPS元素符号映射，其他格式为None

    # 判断输入文件类型
    # === VASP ===
    if "POSCAR" in input1.name.upper() or "POSCAR" in input2.name.upper():
        print("Detected VASP input.")
        poscar = input1 if "POSCAR" in input1.name.upper() else input2
        xdatcar = input2 if poscar == input1 else input1
        u, atom_types = vasp.process_vasp(str(xdatcar), str(poscar))
        prefix = "vasp"
    # === LAMMPS ===
    elif input1.suffix == ".data" or input2.suffix == ".data":
        print("Detected LAMMPS input.")
        data = input1 if input1.suffix == ".data" else input2
        traj = input2 if data == input1 else input1
        u, atom_types, symbol_map = lmp.process_lmp(str(data), str(traj))
        prefix = "lmp"
    # === GROMACS ===
    elif input1.suffix == ".tpr" or input2.suffix == ".tpr":
        print("Detected GROMACS input.")
        tpr = input1 if input1.suffix == ".tpr" else input2
        traj = input2 if tpr == input1 else input1
        if traj.suffix not in [".xtc", ".trr"]:
            #print("GROMACS trajectory must be .xtc or .trr")
            sys.exit(1)
        u, atom_types = gmx.process_gmx(traj, tpr)
        prefix = "gmx"

    else:
        print("Cannot determine input type. Please check file extensions or names.")
        sys.exit(1)

     # === 计算 r_max ===
    box_lengths = u.dimensions[:3]
    min_length = min(box_lengths)
    if min_length > 30:
        r_max = None  # 让MDAnalysis用默认值(15)
    else:
        r_max = min_length * 0.5 * 0.95
    #print({"min_box_lengths": min(box_lengths), "r_max": f"{r_max:.3f}"})

    # 可用原子类型
    print("Atom types:", sorted(set(atom_types)))
    rdf_data = {}  

    #交互式选择计算所有原子对/特定原子对
    if yes_or_no("Calculate RDF for all atom pairs?"):
        # LAMMPS用 type id，其它用 type label
        type_ids = list(symbol_map.keys()) if symbol_map else sorted(set(atom_types))

        for i in range(len(type_ids)):
            for j in range(i, len(type_ids)):
                t1, t2 = type_ids[i], type_ids[j]
                g1 = u.select_atoms(f"type {t1}")
                g2 = u.select_atoms(f"type {t2}")
                pair_name, bins, rdf = rdf_calc.compute_rdf(
                    g1, g2, r_max=r_max, type_labels=symbol_map
                )
                if pair_name is None or bins is None or rdf is None:
                    print(f"Skipping empty or invalid group pair: {t1}-{t2}")
                    continue
                rdf_data[pair_name] = (bins, rdf)

        rdf_calc.compute_multiple_rdfs(
            rdf_data, output_file=f"outputs/{prefix}_rdf_all.png"
        )

    else:
        # 对 LAMMPS，把 atom_types（int）转成字符串，传给选择函数
        # VASP/GMX 本来就是字符串，无需改变
        if symbol_map:
            atom_type_options = list(map(str, atom_types))
        else:
            atom_type_options = sorted(set(atom_types))

        g1_type = get_user_selection("Select first atom type", atom_type_options)
        g2_type = get_user_selection("Select second atom type", atom_type_options)
        
        g1 = u.select_atoms(f"type {g1_type}")
        g2 = u.select_atoms(f"type {g2_type}")
        pair_name, bins, rdf = rdf_calc.compute_rdf(
        g1, g2, r_max=r_max, type_labels=symbol_map
    )

        rdf_calc.compute_multiple_rdfs(
            {pair_name: (bins, rdf)},
            output_file=f"outputs/{prefix}_rdf_{g1_type}_{g2_type}.png",
    )


if __name__ == '__main__':
    main()

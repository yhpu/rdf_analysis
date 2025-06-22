from MDAnalysis.analysis.rdf import InterRDF
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np
from MDAnalysis.analysis.rdf import InterRDF

def compute_rdf(group1, group2, r_max=None, pair_name=None, type_labels=None):
    # 跳过空 group
    if len(group1) == 0 or len(group2) == 0:
        return None, None, None
    
    if r_max is not None:
        rdf = InterRDF(group1, group2, range=(0.0, r_max))
    else:
        rdf = InterRDF(group1, group2)
    rdf.run()

    if pair_name is None:
        try:
            t1 = int(group1.types[0])
            t2 = int(group2.types[0])
            if type_labels:
                label1 = type_labels.get(t1, str(t1))
                label2 = type_labels.get(t2, str(t2))
                pair_name = f"{label1}-{label2}"
            else:
                pair_name = f"{t1}-{t2}"
        except:
            pair_name = f"{group1.types[0]}-{group2.types[0]}"

    return pair_name, rdf.results.bins, rdf.results.rdf


def compute_multiple_rdfs(rdf_dict, output_file, r_max=None):
    plt.figure()
    output_path = Path(output_file)
    output_dir = output_path.parent
    output_dir.mkdir(parents=True, exist_ok=True)

    plt.figure()
    for pair_name, (bins, rdf) in rdf_dict.items():
        plt.plot(bins, rdf, label=pair_name)

    # 保存dat文件，名字和png同名，只改后缀
    data_file = output_path.with_suffix(".dat")
    np.savetxt(data_file, np.column_stack((bins, rdf)),
            header=f"# {pair_name} RDF\n# Distance(A) g(r)")
    plt.xlabel("r (Angstrom)")
    plt.ylabel("g(r)")
    plt.xlim(0.3, 10)
    plt.ylim(0, 60)
    plt.title("Radial Distribution Function")
    plt.legend()
    plt.savefig(output_file)
    plt.close()
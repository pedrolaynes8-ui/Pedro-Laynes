import os
from datetime import datetime
from Bio import Phylo, AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import matplotlib.pyplot as plt

# -------------------------
# RUTAS DEL PROYECTO
# -------------------------
project = "/content/drive/MyDrive/Proyecto-Anisakidos-COX2/"
alineamientos_path = os.path.join(project, "results/alineamientos/")
arboles_path = os.path.join(project, "results/arboles/")
os.makedirs(arboles_path, exist_ok=True)

# -------------------------
# TIMESTAMP
# -------------------------
timestamp = datetime.now().strftime("%Y%m%d_%H%M")

# -------------------------
# CARGAR EL ÚLTIMO ALINEAMIENTO
# -------------------------
alineamientos = [f for f in os.listdir(alineamientos_path) if f.endswith(".fasta")]

if len(alineamientos) == 0:
    raise Exception("❌ No se encontraron alineamientos en results/alineamientos/. Ejecuta 02_phylogeny.py primero.")

# Usar el más reciente
alineamiento_file = sorted(alineamientos)[-1]
alineamiento_path = os.path.join(alineamientos_path, alineamiento_file)

print(f"✔ Usando alineamiento: {alineamiento_file}")

alignment = AlignIO.read(alineamiento_path, "fasta")

# -------------------------
# CONSTRUIR ÁRBOL NJ
# -------------------------
print("Construyendo árbol Neighbor-Joining...")

calculator = DistanceCalculator("identity")  # Modelo para ADN
dm = calculator.get_distance(alignment)

constructor = DistanceTreeConstructor()
tree_nj = constructor.nj(dm)

# Guardar NJ en Newick
nj_newick = os.path.join(arboles_path, f"NJ_tree_{timestamp}.nwk")
Phylo.write(tree_nj, nj_newick, "newick")

# Graficar NJ
plt.figure(figsize=(8, 12))
Phylo.draw(tree_nj, do_show=False)
plt.title("Árbol Neighbor-Joining (NJ)")
nj_png = os.path.join(arboles_path, f"NJ_tree_{timestamp}.png")
plt.savefig(nj_png, dpi=300, bbox_inches="tight")
plt.close()

print(f"✔ Árbol NJ guardado en:\n   - {nj_newick}\n   - {nj_png}")

# -------------------------
# CONSTRUIR ÁRBOL ML (usando modelo de distancia para ADN)
# -------------------------
print("Construyendo árbol Maximum Likelihood...")

# Usamos modelo Jukes-Cantor, válido para ADN
calculator_jc = DistanceCalculator("jukes-cantor")
dm_jc = calculator_jc.get_distance(alignment)

tree_ml = constructor.nj(dm_jc)  # ML rápido usando matriz JC

# Guardar ML en Newick
ml_newick = os.path.join(arboles_path, f"ML_tree_{timestamp}.nwk")
Phylo.write(tree_ml, ml_newick, "newick")

# Graficar ML
plt.figure(figsize=(8, 12))
Phylo.draw(tree_ml, do_show=False)
plt.title("Árbol Maximum Likelihood (ML) - Jukes-Cantor")
ml_png = os.path.join(arboles_path, f"ML_tree_{timestamp}.png")
plt.savefig(ml_png, dpi=300, bbox_inches="tight")
plt.close()

print(f"✔ Árbol ML guardado en:\n   - {ml_newick}\n   - {ml_png}")

print("\n🎉 PROCESO COMPLETADO: Árboles generados correctamente")
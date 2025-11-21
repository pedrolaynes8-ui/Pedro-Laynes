import os
from Bio import SeqIO
from Bio.Align.Applications import MafftCommandline
from datetime import datetime

# --- Rutas ---
project = "/content/drive/MyDrive/Proyecto-Anisakidos-COX2/"
processed_path = os.path.join(project, "data/processed/")
alineamientos_path = os.path.join(project, "results/alineamientos/")
os.makedirs(alineamientos_path, exist_ok=True)

# --- Timestamp para archivos ---
timestamp = datetime.now().strftime("%Y%m%d_%H%M")

# --- Archivos FASTA a alinear ---
archivos = [os.path.join(processed_path, f) for f in os.listdir(processed_path) if f.endswith(".fasta")]

# --- Concatenar en un solo archivo ---
concatenado = os.path.join(alineamientos_path, f"secuencias_a_alinear_{timestamp}.fasta")
with open(concatenado, "w") as f:
    for archivo in archivos:
        seq = SeqIO.read(archivo, "fasta")
        SeqIO.write(seq, f, "fasta")

# --- Ejecutar MAFFT ---
mafft_cline = MafftCommandline(input=concatenado)
print("Ejecutando MAFFT...")
stdout, stderr = mafft_cline()
alineado_file = os.path.join(alineamientos_path, f"alineamiento_COX2_{timestamp}.fasta")
with open(alineado_file, "w") as f:
    f.write(stdout)
print(f"✔ Alineamiento generado: {alineado_file}")

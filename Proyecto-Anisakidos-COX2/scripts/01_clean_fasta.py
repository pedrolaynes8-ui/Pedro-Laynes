from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import PairwiseAligner
import os
from datetime import datetime

# --- Rutas del proyecto ---
project = "/content/drive/MyDrive/Proyecto-Anisakidos-COX2/"
raw_path = os.path.join(project, "data/raw/")
processed_path = os.path.join(project, "data/processed/")
os.makedirs(processed_path, exist_ok=True)

# --- Funciones ---
def leer_ab1(ruta):
    registro = SeqIO.read(ruta, "abi")
    return registro.seq

def generar_consenso(seqF, seqR):
    aligner = PairwiseAligner()
    aligner.mode = "global"
    seqR_rc = seqR.reverse_complement()
    aln = aligner.align(seqF, seqR_rc)[0]
    cons = []
    for bF, bR in zip(aln.seqA, aln.seqB):
        if bF == bR:
            cons.append(bF)
        else:
            cons.append("N")
    return Seq("".join(cons))

def limpiar_seq(seq):
    return Seq(str(seq).strip("N"))

# --- Timestamp para evitar sobreescritura ---
timestamp = datetime.now().strftime("%Y%m%d_%H%M")

# --- Procesamiento de muestras ---
muestras = ["453","454","455","456"]

for code in muestras:
    fwd = os.path.join(raw_path, f"{code}_Cox2-211F.ab1")
    rev = os.path.join(raw_path, f"{code}_Cox2-211R.ab1")
    seqF = leer_ab1(fwd)
    seqR = leer_ab1(rev)
    consenso = generar_consenso(seqF, seqR)
    limpio = limpiar_seq(consenso)
    out_fasta = os.path.join(processed_path, f"{code}_COX2_clean_{timestamp}.fasta")
    SeqIO.write(SeqRecord(limpio, id=code+"_clean", description="COX2 consensus"), out_fasta, "fasta")
    print(f"✔ {out_fasta} generado")

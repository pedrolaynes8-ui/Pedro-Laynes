import os
import time
import csv
from datetime import datetime
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO

# -------------------------
# CONFIGURACIÓN
# -------------------------
project = "/content/drive/MyDrive/Proyecto-Anisakidos-COX2/"
processed_path = os.path.join(project, "data/processed/")
blast_results_path = os.path.join(project, "results/blast/")
os.makedirs(blast_results_path, exist_ok=True)

email = "tu_correo@dominio.com"  # <<<< CAMBIA ESTO (necesario para NCBI)

timestamp = datetime.now().strftime("%Y%m%d_%H%M")
output_csv = os.path.join(blast_results_path, f"BLAST_results_{timestamp}.csv")

# -------------------------
# PREPARAR CSV DE SALIDA
# -------------------------
header = ["Codigo", "Accesion", "Especie", "Porc_identidad", "Cobertura", "Evalue", "Score"]
csvfile = open(output_csv, "w", newline="")
writer = csv.writer(csvfile)
writer.writerow(header)

# -------------------------
# PROCESAR CADA FASTA LIMPIO
# -------------------------
fastas = [f for f in os.listdir(processed_path) if f.endswith(".fasta")]

if len(fastas) == 0:
    raise Exception("❌ No se encontraron FASTA en data/processed/. Ejecuta 01_clean_fasta.py primero.")

print(f"🔍 Ejecutando BLAST para {len(fastas)} secuencias...\n")

for fasta in fastas:
    codigo = fasta.split("_")[0]  # Ej: 453
    fasta_path = os.path.join(processed_path, fasta)
    seq = SeqIO.read(fasta_path, "fasta")
    
    print(f"🚀 Ejecutando BLAST para {codigo} ...")
    
    # Ejecutar BLAST
    result_handle = NCBIWWW.qblast(
        program="blastn",
        database="nt",
        sequence=seq.seq,
        entrez_query="Anisakidae[Organism]",
        hitlist_size=10,
        format_type="XML",
        expect=0.001
    )
    
    # Guardar XML para referencia
    xml_output = os.path.join(blast_results_path, f"{codigo}_BLAST_{timestamp}.xml")
    with open(xml_output, "w") as xml_file:
        xml_file.write(result_handle.read())
    
    result_handle.close()

    # Analizar resultados
    with open(xml_output) as result:
        blast_record = NCBIXML.read(result)

    if len(blast_record.alignments) == 0:
        print(f"⚠ No se encontraron hits para {codigo}")
        writer.writerow([codigo, "-", "-", "-", "-", "-", "-"])
        continue

    # Mejor hit
    best_hit = blast_record.alignments[0]
    hsp = best_hit.hsps[0]

    accesion = best_hit.accession
    descripcion = best_hit.hit_def
    especie = descripcion.split(" ", 1)[1] if " " in descripcion else descripcion

    identidad = round((hsp.identities / hsp.align_length) * 100, 2)
    cobertura = round((hsp.align_length / len(seq)) * 100, 2)
    evalue = hsp.expect
    score = hsp.score

    # Guardar en CSV
    row = [codigo, accesion, especie, identidad, cobertura, evalue, score]
    writer.writerow(row)

    print(f"✔ {codigo}: {especie} ({identidad}% identidad)")

    # Pausa para no saturar NCBI
    time.sleep(2)

csvfile.close()

print("\n🎉 BLAST COMPLETADO")
print(f"📄 Resultados guardados en: {output_csv}")
print(f"📁 Archivos XML en: {blast_results_path}")

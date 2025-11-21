
# Caracterización molecular de larvas de anisákidos (Nematoda: Anisakidae) en peces comerciales de Loreto - Perú

## Descripción
Este repositorio contiene el flujo bioinformático aplicado a larvas de nemátodos de la familia Anisakidae obtenidas de peces comerciales en Loreto, Perú.
El análisis se basa en la secuenciación Sanger del gen mitocondrial COX2.

## Hipótesis
**H1:** El uso de herramientas moleculares permite identificar e inferir la diversidad genética de nemátodos Anisakidae presentes en peces comerciales en Loreto.

## Objetivos
- Caracterización molecular de larvas de Anisakidae mediante secuenciación del gen COX2.
- Identificación taxonómica mediante análisis filogenético.
- Evaluación de diversidad genética intra e interespecífica.

## Muestras y secuencias
- 100 peces analizados, 20 especies diferentes.
- 4 larvas procesadas: 453, 454, 455, 456.
- Secuencias Forward y Reverse, procesadas a FASTA limpio.

## Scripts y pipeline
- `/scripts/01_clean_fasta.py` : limpieza y generación de secuencias consenso
- `/scripts/02_phylogeny.py` : alineamiento múltiple y preparación de árboles
- Flujo resumido:
  1. Lectura de archivos AB1
  2. Limpieza y generación de consenso
  3. BLAST preliminar
  4. Alineamiento múltiple (MAFFT/MUSCLE)
  5. Árboles filogenéticos (NJ, ML, Bayesian)
  6. Análisis de diversidad genética

## Estructura del repositorio
/data/ --> datos crudos y procesados
/scripts/ --> scripts de análisis
/results/ --> alineamientos y árboles
/docs/ --> informes y pipeline documentado
/notebooks/ --> notebooks de Google Colab

# 🧬 Fastq-Architect v7.0 (Clinical & Swarm Edition)

**Orquestador de Eficiencia Computacional y Hardware-Aware (OE-HA)** para el procesamiento masivo de datos genómicos con trazabilidad integral.

## 🌟 Propósito
Este sistema automatiza la purificación de archivos FASTQ, detectando de forma autónoma la topología de la secuencia (Short-Read vs Long-Read) y gestionando la carga de trabajo sobre arquitecturas multi-hilo (optimizado para Ryzen Zen 4). La versión 7.0 introduce estándares clínicos de registro y normalización.

## 🚀 Capacidades de Ingeniería (v7.0)
- **Trazabilidad Maestra (Logging):** Generación automática de `orquestador_master.log`. Un registro inmutable de cada éxito, fallo y decisión del sistema para auditorías posteriores.
- **Normalización de IDs por Regex:** Uso de expresiones regulares para limpiar sufijos técnicos de máquinas Illumina (ej. `_L001_R1_001`) y obtener el ID biológico puro del paciente.
- **Modo Enjambre (Batch):** Escaneo recursivo de directorios para el procesamiento desatendido de cohortes completas de pacientes.
- **Auditoría de Almacenamiento Pre-vuelo:** Suma el peso de la cola de trabajo y verifica el espacio libre en disco para prevenir colapsos por I/O Overflow.
- **Conciencia de Hardware:** Ajuste dinámico de hilos (threads) basado en la CPU del sistema anfitrión.

## 🛠️ Requisitos del Ecosistema
- `fastp` (Motor C++ para Short-Reads)
- `chopper` (Motor Rust para Long-Reads)
- `python` (>= 3.6)

```bash
# Instalación recomendada
micromamba install fastp chopper -c conda-forge -c bioconda -c defaults

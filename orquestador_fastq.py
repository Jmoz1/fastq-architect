import sys
import os
import shutil
import subprocess
import json
import gzip
import re
import logging
from pathlib import Path
from collections import defaultdict

# 1. SISTEMA DE TRAZABILIDAD (Caja Negra)
logging.basicConfig(
    filename='orquestador_master.log', 
    level=logging.INFO,
    format='%(asctime)s - [%(levelname)s] - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)

def auditar_dependencias():
    faltantes = []
    if shutil.which("fastp") is None: faltantes.append("fastp")
    if shutil.which("chopper") is None: faltantes.append("chopper")
    if faltantes:
        logging.critical("Motores C++/Rust no detectados.")
        print("\n[FATAL] Motores C++/Rust no detectados. Ejecuta:")
        print("micromamba install fastp chopper -c conda-forge -c bioconda -c defaults\n")
        sys.exit(1)

def obtener_hilos_optimos():
    return str(min(16, os.cpu_count() or 4))

def auditar_almacenamiento(archivos):
    peso_total_bytes = sum(f.stat().st_size for f in archivos)
    espacio_libre = shutil.disk_usage(Path.cwd()).free
    espacio_requerido = peso_total_bytes * 3
    
    print("\n[AUDITORÍA I/O] Verificando fricción de estado sólido...")
    print(f" - Peso crudo: {peso_total_bytes / (1024**3):.2f} GB")
    print(f" - Espacio libre: {espacio_libre / (1024**3):.2f} GB")
    
    if espacio_libre < espacio_requerido:
        logging.critical(f"Riesgo I/O Overflow. Espacio insuficiente ({espacio_libre} bytes).")
        print(f"[FATAL] Riesgo de colapso de disco (I/O Overflow).")
        sys.exit(1)
    logging.info(f"Auditoría I/O superada. Peso total: {peso_total_bytes / (1024**3):.2f} GB.")
    print("[OK] Espacio en disco aprobado.")

def extraer_id_puro(nombre_archivo):
    """
    2. FILTRO MATEMÁTICO (Expresiones Regulares)
    Amputa todas las colas de las máquinas Illumina (ej. _L001_R1_001.fastq.gz)
    o de Nanopore para obtener el ID puro del paciente.
    """
    # Elimina cualquier cosa desde _R1, _R2, _1, _2 o extensiones fastq
    patron = r'(_S\d+)?(_L\d+)?(_R?[12])?(_\d+)?\.f(ast)?q(\.gz)?$'
    return re.sub(patron, '', nombre_archivo)

def agrupar_pacientes(directorio):
    archivos = list(directorio.glob("*.fastq*")) + list(directorio.glob("*.fq*"))
    archivos = [f for f in archivos if f.is_file() and not f.name.startswith("limpio_")]
    
    if not archivos:
        logging.warning("No se detectaron archivos FASTQ en la ejecución de lote.")
        print("[AVISO] No se detectaron archivos FASTQ en el directorio actual.")
        sys.exit(0)
        
    pacientes = defaultdict(list)
    for f in archivos:
        base_limpia = extraer_id_puro(f.name)
        pacientes[base_limpia].append(f)
        
    for p in pacientes: pacientes[p].sort()
    return pacientes, archivos

def detectar_topologia_secuencia(ruta_archivo):
    abrir = gzip.open if str(ruta_archivo).endswith('.gz') else open
    modo = 'rt' if str(ruta_archivo).endswith('.gz') else 'r'
    longitudes = []
    try:
        with abrir(ruta_archivo, modo) as f:
            for _ in range(100):
                id_seq = f.readline()
                if not id_seq: break
                dna = f.readline().strip()
                f.readline()
                f.readline()
                if dna: longitudes.append(len(dna))
    except Exception as e:
        logging.error(f"Fallo en sensor topológico para {ruta_archivo}: {e}")
        sys.exit(1)
    media = sum(longitudes) / len(longitudes) if longitudes else 0
    return "LONG-READ" if media > 500 else "SHORT-READ"

def parsear_json_fastp(ruta_json):
    try:
        with open(ruta_json, 'r') as f: data = json.load(f)
        antes = data['summary']['before_filtering']['total_reads']
        despues = data['summary']['after_filtering']['total_reads']
        q30 = data['summary']['after_filtering']['q30_rate'] * 100
        return antes, despues, q30
    except: return 0, 0, 0

def procesar_paciente(id_paciente, archivos_paciente, hilos):
    in1 = archivos_paciente[0]
    topologia = detectar_topologia_secuencia(in1)
    dir_salida = Path.cwd() / id_paciente
    dir_salida.mkdir(exist_ok=True)
    
    msg_inicio = f"INICIANDO: {id_paciente} | {topologia} | {len(archivos_paciente)} archivo(s)"
    logging.info(msg_inicio)
    print("\n" + "-"*50)
    print(f"🚀 {msg_inicio}")
    
    if topologia == "SHORT-READ":
        out1 = dir_salida / f"limpio_{in1.name}"
        html_rep, json_rep = dir_salida / f"{id_paciente}.html", dir_salida / f"{id_paciente}.json"
        
        if len(archivos_paciente) == 2:
            in2 = archivos_paciente[1]
            out2 = dir_salida / f"limpio_{in2.name}"
            cmd = ["fastp", "-i", str(in1), "-I", str(in2), "-o", str(out1), "-O", str(out2), "--detect_adapter_for_pe", "--thread", hilos, "-h", str(html_rep), "-j", str(json_rep)]
        else:
            cmd = ["fastp", "-i", str(in1), "-o", str(out1), "--thread", hilos, "-h", str(html_rep), "-j", str(json_rep)]
            
        proceso = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if proceso.returncode != 0: 
            logging.error(f"Fallo en motor C++ para {id_paciente}.")
            print(f"[ERROR] Fallo en {id_paciente}. Saltando al siguiente.")
            return
            
        antes, despues, q30 = parsear_json_fastp(json_rep)
        resultado = f"Supervivencia: {despues:,} / {antes:,} | Q30: {q30:.2f}%"
        logging.info(f"Completado {id_paciente}. {resultado}")
        print(f"   ↳ {resultado}")
        
    elif topologia == "LONG-READ":
        out1 = dir_salida / f"limpio_{id_paciente}.fastq.gz"
        comando_leer = "zcat" if str(in1).endswith('.gz') else "cat"
        cmd_bash = f"{comando_leer} '{in1}' | chopper -q 10 -l 500 --threads {hilos} | gzip > '{out1}'"
        
        proceso = subprocess.run(cmd_bash, shell=True, executable='/bin/bash', stderr=subprocess.PIPE)
        if proceso.returncode != 0: 
            logging.error(f"Fallo en motor Rust para {id_paciente}.")
            print(f"[ERROR] Fallo en {id_paciente}. Saltando al siguiente.")
            return
        logging.info(f"Completado {id_paciente} [Nanopore].")
        print(f"   ↳ Curación Nanopore (Q>10, L>500) completada.")
        
    print(f"✅ Resultados en: {dir_salida.name}/")

if __name__ == "__main__":
    logging.info("--- INICIO DE SESIÓN DE ORQUESTADOR ---")
    auditar_dependencias()
    hilos_activos = obtener_hilos_optimos()
    
    if len(sys.argv) > 1:
        archivos_input = [Path(f) for f in sys.argv[1:]]
        auditar_almacenamiento(archivos_input)
        id_manual = extraer_id_puro(archivos_input[0].name)
        procesar_paciente(id_manual, archivos_input, hilos_activos)
        
    else:
        print("="*50)
        print("🤖 MODO ENJAMBRE ACTIVADO")
        print("="*50)
        pacientes_dict, todos_los_archivos = agrupar_pacientes(Path.cwd())
        
        print(f"[!] Detectados {len(pacientes_dict)} pacientes ({len(todos_los_archivos)} archivos).")
        logging.info(f"Modo Enjambre: {len(pacientes_dict)} pacientes en cola.")
        auditar_almacenamiento(todos_los_archivos)
        
        for i, (id_paciente, archivos_paciente) in enumerate(pacientes_dict.items(), 1):
            print(f"\n[PROCESO EN COLA: {i}/{len(pacientes_dict)}]")
            procesar_paciente(id_paciente, archivos_paciente, hilos_activos)
            
        print("\n" + "="*50)
        print("🏆 MODO ENJAMBRE FINALIZADO.")
        print("="*50)
    logging.info("--- FIN DE SESIÓN ---")

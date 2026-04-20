import sys
import os
import shutil
import subprocess
import json
import gzip
import re
import logging
import concurrent.futures
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
    faltantes = [cmd for cmd in ["fastp", "chopper"] if shutil.which(cmd) is None]
    if faltantes:
        logging.critical(f"Motores no detectados: {faltantes}")
        sys.stderr.write("\n[FATAL] Motores C++/Rust no detectados. Ejecuta:\n")
        sys.stderr.write("micromamba install fastp chopper -c conda-forge -c bioconda -c defaults\n\n")
        sys.exit(1)

def obtener_hilos_optimos(total_tareas=1):
    cpu_cores = os.cpu_count() or 4
    if total_tareas > 1:
        # En modo enjambre, limita los hilos por subproceso para correr múltiples pacientes en paralelo
        return str(max(2, cpu_cores // min(total_tareas, cpu_cores // 2)))
    return str(min(16, cpu_cores))

def auditar_almacenamiento(archivos):
    peso_total_bytes = sum(f.stat().st_size for f in archivos)
    espacio_libre = shutil.disk_usage(Path.cwd()).free
    espacio_requerido = peso_total_bytes * 3
    
    sys.stdout.write("\n[AUDITORÍA I/O] Verificando fricción de estado sólido...\n")
    sys.stdout.write(f" - Peso crudo: {peso_total_bytes / (1024**3):.2f} GB\n")
    sys.stdout.write(f" - Espacio libre: {espacio_libre / (1024**3):.2f} GB\n")
    
    if espacio_libre < espacio_requerido:
        logging.critical(f"Riesgo I/O Overflow. Espacio insuficiente ({espacio_libre} bytes).")
        sys.stderr.write(f"[FATAL] Riesgo de colapso de disco (I/O Overflow).\n")
        sys.exit(1)
    
    logging.info(f"Auditoría I/O superada. Peso total: {peso_total_bytes / (1024**3):.2f} GB.")
    sys.stdout.write("[OK] Espacio en disco aprobado.\n")

def extraer_id_puro(nombre_archivo):
    patron = r'(_S\d+)?(_L\d+)?(_R?[12])?(_\d+)?\.f(ast)?q(\.gz)?$'
    return re.sub(patron, '', nombre_archivo)

def agrupar_pacientes(directorio):
    archivos = [f for f in directorio.iterdir() if f.is_file() and (f.suffix in ['.gz', '.fastq', '.fq']) and not f.name.startswith("limpio_")]
    
    if not archivos:
        logging.warning("No se detectaron archivos FASTQ en la ejecución de lote.")
        sys.stdout.write("[AVISO] No se detectaron archivos FASTQ en el directorio actual.\n")
        sys.exit(0)
        
    pacientes = defaultdict(list)
    for f in archivos:
        base_limpia = extraer_id_puro(f.name)
        pacientes[base_limpia].append(f)
        
    for p in pacientes: 
        pacientes[p].sort()
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
        with open(ruta_json, 'r') as f: 
            data = json.load(f)
        antes = data['summary']['before_filtering']['total_reads']
        despues = data['summary']['after_filtering']['total_reads']
        q30 = data['summary']['after_filtering']['q30_rate'] * 100
        return antes, despues, q30
    except (json.JSONDecodeError, FileNotFoundError, KeyError) as e:
        logging.error(f"Corrupción de I/O en lectura de JSON de reporte: {e}")
        return -1, -1, 0.0

def procesar_paciente(id_paciente, archivos_paciente, hilos):
    in1 = archivos_paciente[0]
    topologia = detectar_topologia_secuencia(in1)
    dir_salida = Path.cwd() / id_paciente
    dir_salida.mkdir(exist_ok=True)
    
    msg_inicio = f"INICIANDO: {id_paciente} | {topologia} | {len(archivos_paciente)} archivo(s)"
    logging.info(msg_inicio)
    
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
            logging.error(f"Fallo en motor C++ para {id_paciente}. STDERR: {proceso.stderr}")
            return f"[ERROR] Fallo en {id_paciente}."
            
        antes, despues, q30 = parsear_json_fastp(json_rep)
        if antes == -1:
            return f"[ERROR] Fallo I/O métricas en {id_paciente}."
            
        resultado = f"Supervivencia: {despues:,} / {antes:,} | Q30: {q30:.2f}%"
        logging.info(f"Completado {id_paciente}. {resultado}")
        return f"✅ {id_paciente}: {resultado}"
        
    elif topologia == "LONG-READ":
        out1 = dir_salida / f"limpio_{id_paciente}.fastq.gz"
        comando_leer = "zcat" if str(in1).endswith('.gz') else "cat"
        cmd_bash = f"{comando_leer} '{in1}' | chopper -q 10 -l 500 --threads {hilos} | gzip > '{out1}'"
        
        proceso = subprocess.run(cmd_bash, shell=True, executable='/bin/bash', stderr=subprocess.PIPE)
        if proceso.returncode != 0: 
            logging.error(f"Fallo en motor Rust para {id_paciente}.")
            return f"[ERROR] Fallo en {id_paciente}."
            
        logging.info(f"Completado {id_paciente} [Nanopore].")
        return f"✅ {id_paciente}: Curación Nanopore completada."

if __name__ == "__main__":
    logging.info("--- INICIO DE SESIÓN DE ORQUESTADOR ---")
    auditar_dependencias()
    
    if len(sys.argv) > 1:
        archivos_input = [Path(f) for f in sys.argv[1:]]
        auditar_almacenamiento(archivos_input)
        id_manual = extraer_id_puro(archivos_input[0].name)
        hilos_activos = obtener_hilos_optimos(total_tareas=1)
        res = procesar_paciente(id_manual, archivos_input, hilos_activos)
        sys.stdout.write(f"{res}\n")
        
    else:
        sys.stdout.write("="*50 + "\n🤖 MODO ENJAMBRE ACTIVADO\n" + "="*50 + "\n")
        pacientes_dict, todos_los_archivos = agrupar_pacientes(Path.cwd())
        total_pacientes = len(pacientes_dict)
        
        sys.stdout.write(f"[!] Detectados {total_pacientes} pacientes ({len(todos_los_archivos)} archivos).\n")
        logging.info(f"Modo Enjambre: {total_pacientes} pacientes en cola.")
        auditar_almacenamiento(todos_los_archivos)
        
        hilos_por_proceso = obtener_hilos_optimos(total_tareas=total_pacientes)
        max_workers = max(1, (os.cpu_count() or 4) // int(hilos_por_proceso))
        
        sys.stdout.write(f"[*] Aceleración por hardware: {max_workers} pacientes simultáneos ({hilos_por_proceso} hilos c/u).\n\n")
        
        with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
            futuros = {
                executor.submit(procesar_paciente, p_id, p_archivos, hilos_por_proceso): p_id 
                for p_id, p_archivos in pacientes_dict.items()
            }
            
            for index, futuro in enumerate(concurrent.futures.as_completed(futuros), 1):
                try:
                    resultado = futuro.result()
                    sys.stdout.write(f"[{index}/{total_pacientes}] {resultado}\n")
                except Exception as exc:
                    id_paciente = futuros[futuro]
                    sys.stderr.write(f"[{index}/{total_pacientes}] Excepción crítica en {id_paciente}: {exc}\n")
            
        sys.stdout.write("\n" + "="*50 + "\n🏆 MODO ENJAMBRE FINALIZADO.\n" + "="*50 + "\n")
    logging.info("--- FIN DE SESIÓN ---")

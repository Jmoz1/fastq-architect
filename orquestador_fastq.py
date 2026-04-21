import sys
import os
import shutil
import subprocess
import json
import gzip
import logging
import concurrent.futures
from pathlib import Path
from collections import defaultdict
from typing import Dict, List, Tuple, Optional

# --- CONFIGURACIÓN DE ALTO RENDIMIENTO (8845HS) ---
logging.basicConfig(
    filename='clinical_orchestrator.log',
    level=logging.INFO,
    format='%(asctime)s - [%(levelname)s] - %(message)s'
)

class FastqArchitect:
    """
    Motor de pre-procesamiento genómico optimizado para entornos de 
    diagnóstico clínico y análisis de variantes.
    """
    
    def __init__(self, threads_max: int = 16):
        self.threads_max = threads_max
        self.check_dependencies(["fastp", "chopper", "samtools"])

    @staticmethod
    def check_dependencies(cmds: List[str]):
        missing = [c for c in cmds if shutil.which(c) is None]
        if missing:
            logging.critical(f"Missing engines: {missing}")
            raise RuntimeError(f"Faltan dependencias: {missing}")

    def get_topology(self, path: Path) -> str:
        """Sensor topológico con manejo de flujo binario."""
        try:
            opener = gzip.open(path, 'rb') if path.suffix == '.gz' else open(path, 'rb')
            with opener as f:
                lengths = []
                for _ in range(400): # 100 reads (4 líneas por read)
                    _ = f.readline() # ID
                    seq = f.readline().strip()
                    _ = f.readline() # +
                    _ = f.readline() # Qual
                    if seq: lengths.append(len(seq))
                avg = sum(lengths) / len(lengths) if lengths else 0
                return "LONG-READ" if avg > 500 else "SHORT-READ"
        except Exception as e:
            logging.error(f"Error en sensor de {path.name}: {e}")
            return "UNKNOWN"

    def run_fastp(self, patient_id: str, files: List[Path], out_dir: Path, threads: int) -> str:
        """Ejecución de Short-Reads (Illumina) con lógica de Variant Scientist."""
        html = out_dir / f"{patient_id}_qc.html"
        json_rep = out_dir / f"{patient_id}_qc.json"
        
        cmd = ["fastp", "--thread", str(threads), "-h", str(html), "-j", str(json_rep)]
        
        # Detección de PE (Paired-End) vs SE (Single-End)
        if len(files) == 2:
            o1, o2 = out_dir / f"clean_{files[0].name}", out_dir / f"clean_{files[1].name}"
            cmd += ["-i", str(files[0]), "-I", str(files[1]), "-o", str(o1), "-O", str(o2), "--detect_adapter_for_pe"]
        else:
            o1 = out_dir / f"clean_{files[0].name}"
            cmd += ["-i", str(files[0]), "-o", str(o1)]

        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            return f"ERR: {result.stderr[:50]}"
        
        # Extracción de métricas Q30 para firma clínica
        with open(json_rep) as f:
            data = json.load(f)
            q30 = data['summary']['after_filtering']['q30_rate'] * 100
            reads = data['summary']['after_filtering']['total_reads']
        return f"PASS | Q30: {q30:.2f}% | Reads: {reads:,}"

    def run_chopper(self, patient_id: str, file: Path, out_dir: Path, threads: int) -> str:
        """Pipeline de Long-Reads (Nanopore) para Clinical Variant Science."""
        out_file = out_dir / f"clean_{patient_id}.fastq.gz"
        
        # Eliminamos shell=True mediante tuberías de subprocess
        cat_cmd = ["zcat" if file.suffix == '.gz' else "cat", str(file)]
        chop_cmd = ["chopper", "-q", "10", "-l", "500", "--threads", str(threads)]
        gzip_cmd = ["gzip"]

        with open(out_file, "wb") as out_f:
            p1 = subprocess.Popen(cat_cmd, stdout=subprocess.PIPE)
            p2 = subprocess.Popen(chop_cmd, stdin=p1.stdout, stdout=subprocess.PIPE)
            p3 = subprocess.Popen(gzip_cmd, stdin=p2.stdout, stdout=out_f)
            p3.wait()
            
        return "PASS | Nanopore Curated"

def process_node(p_id: str, p_files: List[Path], threads: int):
    """Nodo de ejecución para ProcessPoolExecutor."""
    arch = FastqArchitect()
    topology = arch.get_topology(p_files[0])
    out_dir = Path.cwd() / p_id
    out_dir.mkdir(exist_ok=True)
    
    if topology == "SHORT-READ":
        return p_id, arch.run_fastp(p_id, p_files, out_dir, threads)
    elif topology == "LONG-READ":
        return p_id, arch.run_chopper(p_id, p_files[0], out_dir, threads)
    return p_id, "ERR: Unknown Topology"

# --- BLOQUE PRINCIPAL (ENJAMBRE) ---
if __name__ == "__main__":
    # 1. Auditoría I/O Inicial
    # [Mantener tu lógica de auditar_almacenamiento aquí]
    
    # 2. Agrupamiento de Pacientes
    # [Mantener tu lógica de agrupar_pacientes aquí]
    pacientes_dict = {} # Mock para el ejemplo, usar tu función agrupar_pacientes
    
    total_p = len(pacientes_dict)
    if total_p == 0: sys.exit(0)

    # Optimización Zen4 (16 hilos)
    hilos_per_p = max(2, 16 // total_p) if total_p < 8 else 2
    workers = 16 // hilos_per_p

    print(f"[*] Ejecución: {workers} procesos en paralelo | {hilos_per_p} hilos/paciente")

    with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
        futures = [executor.submit(process_node, pid, files, hilos_per_p) 
                   for pid, files in pacientes_dict.items()]
        
        for f in concurrent.futures.as_completed(futures):
            pid, status = f.result()
            print(f"[REPORT] {pid}: {status}")

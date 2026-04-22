# 🧬 Fastq-Architect

**Hardware-aware FASTQ preprocessing orchestrator for short-read and long-read sequencing data.**

Fastq-Architect is a lightweight orchestration framework designed to transform raw sequencing data (FASTQ) into quality-controlled datasets suitable for downstream variant analysis workflows. It integrates established bioinformatics tools (`fastp`, `chopper`) into a unified, automated pipeline with a focus on reproducibility, data integrity, and efficient resource utilization.

---

## 🚀 Core Features

### 🔹 Automatic Read Topology Detection
- Infers sequencing type (short-read vs long-read) based on sequence length distribution.
- Enables dynamic dispatch of appropriate preprocessing engines.

### 🔹 Dual Engine Processing
- **Short-Reads (Illumina)** → `fastp`
  - Adapter trimming
  - Quality filtering
  - QC report generation (HTML + JSON)
  - Extraction of Q30 and read survival metrics

- **Long-Reads (Oxford Nanopore)** → `chopper`
  - Length filtering
  - Quality thresholding (`Q > 10`)
  - Stream-based processing for large files

---

### 🔹 Parallel Processing (Process-Level)
- Multi-sample execution using `ProcessPoolExecutor`
- Dynamic allocation of threads per sample
- Designed for multi-core CPU environments

---

### 🔹 QC-Oriented Output
- Structured metrics extracted from preprocessing:
  - Q30 rate
  - Total reads after filtering
- Designed to support downstream validation workflows

---

### 🔹 Environment Validation
- Pre-flight dependency checks (`fastp`, `chopper`, `samtools`)
- Ensures required tools are available before execution

---

## 🛠️ Technology Stack

| Component | Role | Language |
|----------|------|----------|
| fastp | Short-read preprocessing | C++ |
| chopper | Long-read filtering | Rust |
| Python 3.12+ | Orchestration layer | Python |
| concurrent.futures | Parallel execution | Standard Library |

---

## 📦 Installation

Recommended environment setup using **micromamba**:

```bash
micromamba create -n fastq_architect \
    fastp chopper python=3.12 \
    -c bioconda -c conda-forge

micromamba activate fastq_architect

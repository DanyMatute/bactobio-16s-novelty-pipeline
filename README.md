# 16S rRNA Novelty Screening Pipeline

A reproducible, production-oriented **Snakemake pipeline** for screening 16S rRNA gene sequencing data to identify **potentially novel bacterial taxa**.

This pipeline is designed to turn raw amplicon sequencing data into **clear, auditable novelty calls** (`KNOWN`, `POTENTIALLY_NOVEL`, `INCONCLUSIVE`), with explicit handling of ambiguity, controls, and biological limitations.

---

## 🎯 Project Goals

- Provide a **robust screening workflow** for detecting candidate novel taxa from 16S data  
- Emphasize **traceability, reproducibility, and conservative decision-making**
- Demonstrate production-quality pipeline design (not just analysis scripts)

This project is intentionally framed as **novelty screening**, not formal species discovery.

---

## 🧠 Design Principles

- **ASV-based inference (DADA2)** for maximum biological resolution and reproducibility
- **Explicit decision logic** with documented thresholds and reasons
- **Strong I/O contracts**: each step produces predictable, versioned artifacts
- **Deterministic inputs** via a manifest file
- **Separation of concerns**: workflow orchestration (Snakemake) vs logic (Python/R)
- **Production-ready structure**: scalable, debuggable, safe to re-run

---

## 🧬 Pipeline Overview

Raw FASTQ

↓

Primer trimming & QC (cutadapt, FastQC)

↓

Denoising & ASV inference (DADA2)

↓

Taxonomic assignment (SILVA trainset)

↓

Novelty decision logic

↓

Per-sample call + metrics


### Core Outputs (per sample)
- `call.json` — final classification + reasons + provenance
- `metrics.tsv` — summary metrics for downstream analysis
- ASV tables and representative sequences
- Logs for each pipeline stage

---

## 📂 Repository Structure
.
├── Snakefile

├── config/

│ └── config.yaml

├── data/

│ └── manifest.tsv

├── workflow/

│ ├── envs/ # Conda environments

│ ├── scripts/ # Python/R logic

│ └── rules/ # Snakemake rules

├── resources/

│ └── db/README.md # Reference DB instructions

└── README.md


**Note:**  
Raw sequencing data, reference databases, logs, and results are intentionally **excluded from version control** and regenerated deterministically.

---

## 📄 Input Manifest (Deterministic Inputs)

All samples are defined explicitly in a tab-delimited manifest:

```tsv
sample        read1                     read2                     sample_type
TEST_F3D146   data/...R1.fastq.gz       data/...R2.fastq.gz       test
NEG_BLANK     data/...R1.fastq.gz       data/...R2.fastq.gz       negative
POS_MOCK      data/...R1.fastq.gz       data/...R2.fastq.gz       positive
```
---
**🧪 Denoising & ASVs (Why ASVs > OTUs)**
* This pipeline uses Amplicon Sequence Variants (ASVs) rather than OTUs.
* Why:
  - ASVs represent exact biological sequences
  - Reproducible across runs and datasets
  - Preserve subtle sequence differences important for novelty detection
  - Avoid arbitrary clustering thresholds
* Denoising models sequencing error to separate true biological variation from noise.

**🧠 Novelty Decision Logic (v0)**

* The pipeline makes a conservative screening call:
* Decision categories
  - KNOWN → At least one ASV confidently assigned to a known genus
  - POTENTIALLY_NOVEL → Sufficient signal, but no genus-level assignment
  - INCONCLUSIVE → Insufficient data, failed QC, or control issues
* Key principles
  - Negative controls are evaluated first
  - Low read depth or too few ASVs → no over-interpretation
  - Novelty indicates candidates for follow-up, not formal species claims
* All decisions include explicit reasons and supporting metrics.

** 🔁 Reproducibility & Provenance**
* Each run is reproducible and auditable via:
  - manifest-defined inputs
  - versioned pipeline code
  - recorded thresholds
  - reference database versioning
  - per-sample provenance metadata
* Outputs can be safely regenerated or compared across runs.

🚀 Running the Pipeline
* Requirements
  - Snakemake
  - Conda / Mamba
  - Linux environment (local, HPC, or cloud)
* Example
```snakemake -j 16```
* The same workflow can scale to clusters or cloud executors using Snakemake profiles.

⚠️ Limitations (By Design)
* 16S cannot resolve all species (e.g. closely related taxa)
* Novelty is relative to reference databases
* Final species confirmation requires:
  - full-length 16S, and/or
  - isolate sequencing, and/or
  - whole-genome data
* This pipeline is intentionally conservative to avoid false novelty claims.

📌 Future Extensions
* Phylogenetic placement for borderline cases
* Cross-database taxonomy comparison
* Integration with isolate/WGS workflows
* Persistent metadata storage (e.g. SQLite/Postgres)
* Automated reprocessing on DB updates

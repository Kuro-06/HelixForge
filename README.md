
---

# HelixForge 🧬

> **A High-Performance, Dependency-Free Genomic Alignment Engine**

(tldr available at the end)

HelixForge is a custom-built, systems-level bioinformatics engine designed for ultra-fast genomic querying and sequence alignment. Built entirely from scratch in **C++17** with zero external bioinformatics dependencies, this engine implements a highly optimized **"Seed-and-Extend"** architecture.

It allows for instantaneous querying of massive reference genomes (e.g., *E. coli*) using a mathematically compressed index, followed by rigorous dynamic programming for precision local alignment and statistical significance calculation.

---

## 🧠 Core Architecture & Algorithmic Pipeline

HelixForge is designed around four mathematically rigorous computational phases, shifting seamlessly from global search heuristics to local precision.

### 1. BWT & FM-Index Compression (Search Space Reduction)

Instead of linearly scanning a multi-million base-pair genome, the reference sequence is compressed using the **Burrows-Wheeler Transform (BWT)** via an in-place Prefix Doubling Suffix Array algorithm O(N log² N).

* **LF-Mapping Engine:** The BWT allows the engine to utilize Last-to-First (LF) mapping, reducing exact string-matching searches to **O(P)** time complexity (where P is the length of the query), rendering the search time entirely independent of the reference genome's massive size.

### 2. Seed-and-Extend Clustering (The Heuristic)

The engine extracts k-mers (seeds) from the query and locates all absolute coordinates in the reference using the FM-Index.

* **1D Collinear Chaining:** The mathematical coordinates are sorted and clustered to identify the highest-density hit regions, successfully filtering out billions of base pairs of noise and isolating the exact biological window where the mutation/alignment likely exists.

### 3. Smith-Waterman Local Alignment (Dynamic Programming)

Once the optimal genomic window is extracted, it is passed into an O(M × N) dynamic programming matrix.

* **Affine Gap Penalties:** Implements sophisticated gap-open and gap-extend penalties (-10/-1) to biologically model major indels (insertions/deletions) without prematurely terminating the alignment.
* **Dual-Alphabet Matrix Routing:** Automatically dynamically routes scoring to **BLOSUM62** for Protein amino acids, and **+1/-1 matrices** for DNA nucleotides.

### 4. Karlin-Altschul Statistics

A raw alignment score is scientifically meaningless without statistical context. The engine calculates the true **E-Value (Expectation Value)** using Karlin-Altschul mathematical constants (λ and K), verifying whether an alignment is a biological homology or mere random background noise.

---

## ⚡ Systems-Level Optimizations (C++)

* **Dynamic Alphabet Reduction:** The engine dynamically analyzes the raw string data and allocates Occurrence Matrix memory *only* for actively used ASCII characters. This drops RAM overhead by **>90%** compared to standard 256-character ASCII allocation arrays.
* **Zero-Copy Substring Slicing:** Utilizes `std::string_view` throughout the Seed-and-Extend pipeline to pass reference windows between functions without triggering costly heap-memory allocations.
* **Binary Serialization:** Mathematical tables (Suffix Arrays, BWTs) are serialized directly into binary `.idx` formats. This allows the software to bypass the O(N log² N) build phase on subsequent runs via direct memory-mapping.
* **REPL State Machine:** The frontend operates as a persistent Read-Eval-Print Loop. Massive genomes are loaded into RAM exactly once, allowing for hundreds of consecutive user queries without I/O disk bottlenecks.

---

## 📊 Performance Benchmarks

*Tested on an Apple Silicon Unified Memory Architecture against a 4.64 Million Base-Pair E. coli Reference Genome.*

> **Architect's Note:** All benchmarks below were recorded on a strict **single-threaded execution model**. These execution times are projected to decrease exponentially upon the integration of parallel multi-threading across modern multi-core architectures.

| Computational Phase | Execution Time | Complexity |
| --- | --- | --- |
| **Initial Index Build** (Prefix Doubling) | ~3.0 seconds | O(N log² N) |
| **Binary `.idx` Load to RAM** | **< 0.02 seconds** | O(1) File I/O |
| **FM-Index K-mer Extraction** (15-mers) | **< 0.01 seconds** | O(P) |
| **Full Pipeline Execution** (Search + Align + Stats) | **< 0.05 seconds** | O(M × N) Local |

---

## 🚀 Quick Start & Usage

### Prerequisites

* C++17 Compliant Compiler (GCC, Clang, Apple Clang)
* CMake 3.10+

### Compilation

Clone the repository and build using CMake:

```bash
git clone https://github.com/yourusername/HelixForge.git
cd HelixForge
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build

```
(must build in release mode,this reduces runtime by ~4 times)
### The Interactive Workbench

HelixForge is operated via an interactive, persistent REPL shell.

```text
user@macbook % ./build/HelixForge

=========================================================
      🧬 HELIXFORGE INTERACTIVE WORKBENCH v1.0
=========================================================
Type 'help' for commands or 'exit' to quit.

# 1. Switch to Protein Mode (Defaults to DNA)
HelixForge [DNA] > mode protein

# 2. Load Reference (Builds .idx if not found)
HelixForge [Protein] > load ref sp_human_p53

# 3. Load Query
HelixForge [Protein] > load query sp_mouse_p53

# 4. Execute the mathematical pipeline
HelixForge [Protein] > align

```

---

## 🗺️ Roadmap & Future Architecture

HelixForge is under active development. The current roadmap focuses on expanding the engine's algorithmic complexity, parallel processing capabilities, and memory safety:

* **Multi-threading Computation Engine:** Implementing parallel execution via `std::thread` pools to partition the FM-Index search space and DP matrices simultaneously. This will exponentially decrease compute times and enable massive batch-processing of thousands of `.fastq` sequence reads.
* **Advanced Memory Guardrails:** Integrating stricter memory allocation parameters, custom smart pointers, and bounds-checking guardlines to provide absolute segmentation fault safety, ensuring zero memory leaks during heavy multithreaded loads.
* **O(N) Suffix Array Construction:** Upgrading the BWT construction from Prefix-Doubling to the industry-standard SA-IS (Suffix Array Induced Sorting) algorithm to handle multi-gigabyte genomes (like the Human Genome) in pure linear time.
* **Dynamic Programming Clustering:** Replacing the current 1D distance heuristic with a 2D Collinear Chaining algorithm with dynamic max-hit values to calculate mathematically optimal, biologically ordered seed coordinates.
* **Custom Matrices & User Parameters:** Allowing researchers to dynamically ingest custom PAM/BLOSUM scoring matrices and fine-tune gap-penalty environments via CLI variables.
* **GUI Integration:** Wrapping the highly optimized C++ backend in a modern, lightweight graphical user interface.

## TLDR
In layman terms without biological context, the program engine can search a query of any size in 0(1) time withing a refference file by building its BWT file once.
(the files for a 4.6 million characters long bacteria DNA was produced in ~3.0 second)
same is used in biology in lot of places (matching dna).
the tool is optimised specifically for DNA and Protiens fasta format files with aditional tools for their matching,scoring and calculating statistical significance.
(build entirely from scratch)


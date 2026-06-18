# HelixForge 🧬

> **A High-Performance, Dependency-Free Genomic Alignment Engine written in C++**

HelixForge is a custom-built, systems-level bioinformatics engine designed for ultra-fast genomic querying and sequence alignment. It implements a highly optimized **"Seed-and-Extend"** architecture, allowing for instantaneous searching of massive reference genomes (such as *E. coli*) using a compressed mathematical index, followed by precision local alignment.

Built entirely from scratch without the use of external bioinformatics libraries, HelixForge serves as a demonstration of low-level systems architecture, memory management, and advanced algorithmic design.

---

## 🏗️ Core Architecture & Algorithms

### 1. Burrows-Wheeler Transform (BWT) & Suffix Arrays
Compresses the reference genome using an in-place **Prefix Doubling** algorithm (O(N log² N)). 
* *Engineering Choice:* This approach was deliberately chosen over O(N log N) Radix/Counting implementations to avoid RAM-thrashing and bypass the "Memory Wall." By using an in-place comparative sort, the engine prioritizes CPU L1/L2 cache locality, making it highly stable across different hardware profiles when processing multi-million-element arrays.

### 2. FM-Index & LF-Mapping Search Engine
Implements a mathematically pure Last-to-First (LF) mapping engine for exact substring matching without scanning the reference string. 
* **Dynamic Alphabet Reduction:** The engine dynamically analyzes the reference sequence to allocate Occurrence Matrix memory *only* for active characters. This prevents the generation of billions of useless zero-values for unused ASCII characters, reducing RAM overhead by >90% and dropping query search times to O(1) relative to genome size.

### 3. Smith-Waterman Local Alignment
A dynamic programming matrix used for the "Extend" phase of the pipeline.
* Features **Affine Gap Penalties** to accurately model biological realities, smoothly bridging large deletion mutations without prematurely terminating the alignment score.

### 4. Automated Sequence Inference
Dynamically scans query inputs to distinguish between DNA nucleotides and Protein amino acids. It automatically adjusts substitution matrices and mismatch penalties on the fly, ensuring biologically accurate scoring without manual user configuration.

---

## 📂 Directory Structure

    HelixForge/
    ├── CMakeLists.txt           # Build configuration
    ├── README.md                # Project documentation
    ├── include/                 # Header files (.hpp)
    │   ├── AlignmentEngine.hpp
    │   ├── AlignmentWorkflow.hpp
    │   ├── Genome.hpp
    │   └── SequenceUtils.hpp
    ├── src/                     # Source files (.cpp)
    │   ├── AlignmentEngine.cpp
    │   ├── AlignmentWorkflow.cpp
    │   ├── Genome.cpp
    │   ├── SequenceUtils.cpp
    │   └── main.cpp
    └── database/                # Local data storage
        └── dna/                 # Place .fasta files here
            ├── sampleA.fasta    # Target Reference Genome (e.g., E. coli)
            └── sampleB.fasta    # Query Sequence

---

## 🚀 Build & Execution

HelixForge uses CMake for cross-platform compilation.

### Prerequisites
* C++17 Compiler (GCC, Clang, or MSVC)
* CMake (3.10+)

### Compilation
Open your terminal at the root of the project and run:
    
    cmake -B build
    cmake --build build

### Usage
Ensure your target sequences (`sampleA.fasta` and `sampleB.fasta`) are placed in the `database/dna/` directory.

    ./build/HelixForge

---

## ⚡ Database Serialization & Performance
*Benchmarks tested on Apple Silicon (M3 / 8GB Unified Memory) against a 4.64 million base pair E. coli genome.*

To minimize repeated computation, HelixForge serializes its mathematical tables. On the initial run with a new reference FASTA, the engine will compute the BWT and Suffix Array. This database is then serialized into a binary `.idx` file. 

Subsequent executions will automatically bypass the calculation phase and load the binary index directly into RAM.

| Operation | Execution Time |
| :--- | :--- |
| **Index Build Time** (First Run) | ~12.0 - 14.0 seconds |
| **Index Load Time** (From `.idx`) | **< 0.05 seconds** |
| **335-bp Query Search & Alignment** | **< 0.10 seconds** |

---

## 🗺️ Future Roadmap
- [ ] **Batch Processing:** Implement reading from multi-query `.fastq` files.
- [ ] **Multi-threading:** Utilize `std::thread` to process multiple alignments concurrently across different CPU cores.
- [ ] **TUI (Terminal User Interface):** Build an interactive command-line menu for dynamic file loading and parameter tweaking.
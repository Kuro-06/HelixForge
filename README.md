# HelixForge 🧬

HelixForge is a high-performance, C++ based bioinformatics alignment pipeline. It implements a two-stage **"Seed-and-Extend" heuristic** designed to rapidly search massive genomic reference sequences and perform highly accurate, statistically robust local alignments.

By decoupling the high-speed search phase from the computationally heavy alignment phase, HelixForge drastically reduces the standard $O(nm)$ time complexity of dynamic programming algorithms while preserving biological accuracy.

---

## 🚀 Architecture & Core Pipeline

The pipeline operates in three distinct phases:

1. **Phase 1: The Search Engine (`KmerIndex`)** The reference genome is chopped into $k$-mers and indexed using an $O(1)$ Hash Map. When a query is submitted, it is parsed into $k$-mers to rapidly locate exact matches ("seeds") in the reference.
2. **Phase 2: Seed Chaining & Extraction**
   Raw hits are clustered using a dynamic **Seed Chaining algorithm**. A Bounding Box is calculated around the densest cluster, and a tiny local neighborhood is mathematically "scissored" out of the massive reference sequence.
3. **Phase 3: High-Precision Alignment (`AlignmentEngine`)**
   The extracted sub-sequence and the query are passed to a highly optimized Smith-Waterman engine, which computes the final alignment score using Affine Gap penalties.

---

## 🧠 Key Algorithmic Optimizations (Interview Highlights)

### 1. Hash-Based $k$-mer Indexing
Instead of running a dynamic programming matrix over an entire genome, the reference is indexed once using `std::unordered_map<std::string, std::vector<int>>`. This allows $O(1)$ average time complexity for lookup operations, bypassing 99.9% of the reference sequence mathematically.

### 2. Dynamic Seed Chaining (Overcoming the "Hula Hoop" Flaw)
A naive sliding window approach fails when biological sequences contain massive insertions (like introns or transposons). HelixForge implements a **Seed Chaining Algorithm** that dynamically links $k$-mer hits based on proximity.
* **Auto-Scaling Gap Threshold:** The `maxGap` allowed between seeds dynamically scales relative to the query length (Query Length / 2), clamping within biological bounds (50 to 3000 bases). This prevents short queries from chaining random noise, while allowing long queries to successfully bridge massive biological insertions.

### 3. Memory-Safe Subsequence Extraction
Extracting data from large memory structures often leads to out-of-bounds segmentation faults. The Bounding Box calculation utilizes `std::max` and `std::min` to enforce strict boundary safety:
* `startCut = std::max(0, bestHit - padding);` ensures we never access negative memory addresses.
* `safeLength = std::min(desiredLength, referenceLength - startCut);` ensures we never bleed past the end of the chromosome.

### 4. Affine Gap Penalty Model
The Smith-Waterman engine does not use a linear gap penalty. It tracks gap openings versus gap extensions independently. This requires tracking three separate matrices simultaneously but yields a biologically accurate alignment that accurately models the reality of grouped mutational events.

---

## 🛠️ API & End-User Functions

### `AlignmentEngine` Class
Handles the heavy matrix mathematics and pathfinding.
* `AlignmentEngine(std::string seqA, std::string seqB)`: Constructor initializing the exact strings to be aligned.
* `void scorePara(int customMatch, int customMismatch, int customGapOpen, int customGapExtend)`: Allows the user to override default scoring parameters.
* `int Waterman()`: Executes the Smith-Waterman algorithm and returns the optimal local alignment score.
* `void Traceback()`: Reconstructs and prints the visual alignment of the two sequences.

### `KmerIndex` Class
Handles the fast-search database and spatial clustering.
* `KmerIndex(int kSize)`: Constructor to set the $k$-mer length (default $k=11$ for DNA).
* `void buildIndex(const std::string& referenceSeq)`: Uses a sliding window to populate the hash map with exact coordinates.
* `std::vector<int> getHits(const std::string& queryKmer)`: Performs the $O(1)$ database lookup.
* `int findBestCluster(const std::string& querySeq)`: Executes the Seed Chaining algorithm, sorts coordinates, and returns the starting integer index of the most biologically relevant alignment zone.

---

## 🗺️ Future Roadmap
* **`Genome` Class Wrapper:** Object-oriented management to cleanly handle biological sequences.
* **`.fasta` File Ingestion:** A custom file parser to ingest real-world genomic data.
* **`.idx` Binary Serialization:** Writing the `std::unordered_map` directly to disk so genomes only need to be indexed once, allowing instant loading into RAM on future runs.

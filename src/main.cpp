#include <iostream>
#include "AlignmentEngine.hpp"
#include "Genome.hpp"
#include "AlignmentWorkflow.hpp"



void fetchFromDatabase(Genome& genomeObj, const std::string& filename) {
    std::string fastaPath = "database/dna/" + filename + ".fasta";
    std::string idxPath   = "database/dna/" + filename + ".idx";

    if (!genomeObj.loadIndex(idxPath)) {
        genomeObj.loadFromFasta(fastaPath);
        genomeObj.buildBWT();
        genomeObj.saveIndex(idxPath);
    }
}

int main() {
    std::cout << "... HelixForge Workbench Initialized ..." << std::endl;

    Genome myReference;
    fetchFromDatabase(myReference, "sampleA"); 

    Genome myQuery;
    std::cout << "\n... Fetching Query 'sampleB' ..." << std::endl;

    std::string queryPath = "database/dna/sampleB.fasta";
    std::cout << "Target FASTA : " << queryPath << std::endl;
    
   myQuery.loadFromFasta("database/dna/sampleB.fasta");
    std::string querySeq = myQuery.getOriginalSequence();

    if (querySeq.empty()) {
        std::cout << "[ERROR] Query sequence is empty or missing from " << queryPath << "\n";
        return 1;
    }

    int kSize = 15; // make seperate functions for this in future
    std::pair<int, int> clusterBounds = AlignmentWorkflow::findBestCluster(myReference, querySeq, kSize);

    if (clusterBounds.first == -1) {
        std::cout << "[RESULT] No viable matches found for alignment.\n";
        return 0;
    }

    std::string fullRefDNA = myReference.getOriginalSequence();
    
    int start = clusterBounds.first;
    int length = clusterBounds.second - clusterBounds.first;

    std::string targetWindow = fullRefDNA.substr(start, length);


    std::cout << "Absolute Coordinates : [" << start << " to " << clusterBounds.second << "] in Reference\n";
    std::cout << "\n--- Final Alignment Window ---\n";
    std::cout << "Reference Window : " << targetWindow << "\n";
    
    std::cout << "Query Sequence   : " << querySeq << "\n";

 
    AlignmentEngine aligner(targetWindow, querySeq); // alignment engine object;
    
    int localScore = aligner.Waterman();
    
    std::cout << ">>> FINAL LOCAL ALIGNMENT SCORE: " << localScore << " <<<\n\n";
    std::cout << "--- Visual Traceback ---\n";
    
    aligner.TraceBack();

    std::cout << "\n(Note: The alignment above occurred within the absolute genomic bounds of [" 
              << start << ", " << clusterBounds.second << "])\n";

    return 0;
}
#include <iostream>
#include "AlignmentEngine.hpp"
#include "FastaReader.hpp"

int main() {
    std::cout << "HelixForge Workbench Initialized!" << std::endl;

    std::string seqA = extractSequenceFromFasta("data/sampleA.fasta");
    std::string seqB = extractSequenceFromFasta("data/sampleB.fasta");
    
    // 2. Pass the real file data into the engine!
    AlignmentEngine engine(seqA,seqB);
    std::cout << engine.Needleman() <<std::endl;
    engine.TraceBack();

    std::cout << engine.Waterman() <<std::endl;
    engine.TraceBack();
    return 0;
}
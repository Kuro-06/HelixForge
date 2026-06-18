#ifndef ALIGNMENT_WORKFLOW_HPP
#define ALIGNMENT_WORKFLOW_HPP

#include <string>
#include <utility>
class Genome;

namespace AlignmentWorkflow {
    
    std::pair<int, int> findBestCluster(Genome& targetGenome, const std::string& querySeq, int k);
}

#endif
#pragma once

#include <string>
#include <utility>
#include <string_view>
class Genome;

namespace AlignmentWorkflow {
    struct ClusterResult {
        bool success;
        int start;
        int end;
        int maxDensity;
        int totalHits;
    };
    
    ClusterResult findBestCluster(Genome& targetGenome, const std::string_view querySeq, int k);
}

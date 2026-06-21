#include "AlignmentWorkflow.hpp"
#include "Genome.hpp"
#include <iostream>
#include <string>
#include <vector>
#include <algorithm> 
#include <utility> 


namespace AlignmentWorkflow {

    ClusterResult findBestCluster(Genome& targetGenome, const std::string_view querySeq, int k) {
        
        ClusterResult result = {false, -1, -1, 0, 0};
        //query too smalll
        if (querySeq.length() < static_cast<size_t>(k)) return result;

        std::vector<int> allHits;
        int limit = querySeq.length() - k;
        const size_t MAX_HITS_PER_SEED = 65;

        for (int i = 0; i <= limit; i++) {
            std::string_view currentKmer= querySeq.substr(i, k);
            std::vector<int> hits = targetGenome.search(currentKmer);

            if (hits.size() > MAX_HITS_PER_SEED) continue;
            allHits.insert(allHits.end(), hits.begin(), hits.end());
        }

        if (allHits.empty()) return result;

        std::sort(allHits.begin(), allHits.end());
        result.totalHits = allHits.size();

        int queryLen = querySeq.length();
        int dynamicGap = queryLen / 2; 

        if (dynamicGap < 50) dynamicGap = 50;       
        if (dynamicGap > 3000) dynamicGap = 3000;   

        // 4. Linear Clustering Algorithm
        int bestChainStart = allHits[0]; 
        int bestChainEnd = allHits[0]; 
        int maxDensity = 0;

        int currentChainStart = allHits[0];
        int currentChainEnd = allHits[0];
        int currentDensity = 1;

        for (size_t i = 1; i < allHits.size(); i++) {
            int currentHit = allHits[i]; 
            
            if (currentHit - currentChainEnd <= dynamicGap) {
                currentChainEnd = currentHit;
                currentDensity++;
            } else {
                if(currentDensity > maxDensity) {
                    maxDensity = currentDensity;
                    bestChainStart = currentChainStart;
                    bestChainEnd = currentChainEnd; 
                }
                currentChainStart = currentHit;
                currentChainEnd = currentHit;
                currentDensity = 1;
            }
        }

       
        if (currentDensity > maxDensity) {
            maxDensity = currentDensity;
            bestChainStart = currentChainStart;
            bestChainEnd = currentChainEnd;
        }

        result.success = true;
        result.start = bestChainStart;
        result.end = bestChainEnd + k;
        result.maxDensity = maxDensity;
        
        return result;
    }

}
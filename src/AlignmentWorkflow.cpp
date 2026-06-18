#include "AlignmentWorkflow.hpp"
#include "Genome.hpp"
#include <iostream>
#include <string>
#include <vector>
#include <algorithm> 
#include <utility> 

namespace AlignmentWorkflow {

    std::pair<int, int> findBestCluster(Genome& targetGenome, const std::string& querySeq, int k) {
        //query too smalll
        if (querySeq.length() < k) return {-1, -1};

        std::vector<int> allHits;
        int limit = querySeq.length() - k;

        std::cout << "[SYSTEM] Extracting " << k << "-mer seeds and querying FM-Index...\n";

        for (int i = 0; i <= limit; i++) {
            std::string currentKmer = querySeq.substr(i, k);
            std::vector<int> hits = targetGenome.search(currentKmer);
            allHits.insert(allHits.end(), hits.begin(), hits.end());
        }

        if (allHits.empty()) {
            return {-1, -1}; 
        }

        std::sort(allHits.begin(), allHits.end());

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

        std::cout << "[RESULT] Densest cluster found with " << maxDensity << " overlapping seeds.\n";
        return {bestChainStart, bestChainEnd + k};
    }

}
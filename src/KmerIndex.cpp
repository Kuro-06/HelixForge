#include "KmerIndex.hpp"
#include <iostream>
#include <algorithm>

// Constructor
KmerIndex::KmerIndex(int kSize) {
    k = kSize;
}

void KmerIndex::buildIndex(const std::string& referenceSeq) {
    indexMap.clear();

    if (referenceSeq.length() < k) {
        std::cerr << "Warning: Reference sequence is shorter than k-mer size." << std::endl;
        return;
    }

    int limit = referenceSeq.length() - k;
    
    for (int i = 0; i <= limit; i++) {
        std::string currentKmer = referenceSeq.substr(i, k);
        indexMap[currentKmer].push_back(i);
    }
}

// The Query Engine
std::vector<int> KmerIndex::getHits(const std::string& queryKmer) {

    if (indexMap.find(queryKmer) != indexMap.end()) {
        return indexMap[queryKmer]; 
    }
    
    return std::vector<int>();
}

// Utility Function
int KmerIndex::getIndexSize() {
    return indexMap.size();
}

//finding the cluster that is the most dense (chaining with dynamic gap length)
int KmerIndex::findBestCluster(const std::string& querySeq) {
    if (querySeq.length() < k) return -1; // Query too small!

    std::vector<int> allHits;

    // finding hits for each possible query and appending to master list
    int limit = querySeq.length() - k;
    for (int i = 0; i <= limit; i++) {
        std::string currentKmer = querySeq.substr(i, k);
        std::vector<int> hits = getHits(currentKmer);
        allHits.insert(allHits.end(), hits.begin(), hits.end());
    }

    // if none found
    if (allHits.empty()) {
        return -1; 
    }

    //sorting hits to detect clusters
    std::sort(allHits.begin(), allHits.end());

   // dynamic gap calculation (50% of the query length)
    int queryLen = querySeq.length();
    int dynamicGap = queryLen / 2; 

    
    if (dynamicGap < 50) dynamicGap = 50;       // Minimum allowed gap
    if (dynamicGap > 3000) dynamicGap = 3000;   // Maximum allowed gap

    //seed chaining
    int bestChainStart = allHits[0];
    int maxDensity = 0;

    //tracking current chain (patch finding algo)
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
            }
            //starting of new chain after breaking of last one
            currentChainStart = currentHit;
            currentChainEnd = currentHit;
            currentDensity = 1;
        }
    }

    //final chain check
    if (currentDensity > maxDensity) {
        maxDensity = currentDensity;
        bestChainStart = currentChainStart;
    }

    return bestChainStart;
}
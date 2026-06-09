#ifndef KMERINDEX_HPP
#define KMERINDEX_HPP

#include <string>
#include <vector>
#include <unordered_map>

class KmerIndex {
private:
    int k;
    std::unordered_map<std::string, std::vector<int>> indexMap;

public:
   
    KmerIndex(int kSize);

    // Build kmer index hash maps
    void buildIndex(const std::string& referenceSeq);

    // Seeding to get hits
    std::vector<int> getHits(const std::string& queryKmer);
    
    // Optional utility to see how big the index got
    int getIndexSize();
    
    // cluster finding
    int findBestCluster(const std::string& querySeq);
};

#endif
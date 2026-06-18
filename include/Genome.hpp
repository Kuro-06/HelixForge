#ifndef GENOME_HPP
#define GENOME_HPP

#include <string>
#include <vector>

class Genome {
private:
    std::string header;
    std::string originalSequence;
    std::string bwtSequence;
    std::vector<int> suffixArray;

    // Pre-calculated FM-Index Tables
    std::vector<int> counts;                     
    std::vector<std::vector<int>> occ;          

    // Internal helper to build the tables
    void generateFMTables();

public:
    Genome() = default;

    // Core Functions
    void loadFromFasta(const std::string& filepath);
    void buildBWT();
    bool saveIndex(const std::string& filepath);
    bool loadIndex(const std::string& filepath);

    std::vector<int> search(const std::string& query);
    

    std::string getHeader() const { return header; }
    std::string getOriginalSequence() const { return originalSequence; }
    std::string getBWTSequence() const { return bwtSequence; }
};

#endif
#pragma once

#include <string>
#include <vector>
#include <string_view>

class Genome {
private:
    std::string header;
    std::string originalSequence;
    std::string bwtSequence;
    std::vector<int> suffixArray;

    bool isProtein;
    int alphabetSize;

    // Pre-calculated FM-Index Tables
    std::vector<int> counts;                     
    std::vector<std::vector<int>> occ;  
    
    uint8_t getBaseIndex(char c) const;

    // Internal helper to build the tables
    void generateFMTables();

public:
    Genome(bool proteinMode = false);

    // Core Functions
    bool loadFromFasta(const std::string& filepath);
    void buildBWT();
    bool saveIndex(const std::string& filepath);
    bool loadIndex(const std::string& filepath);

    std::vector<int> search(const std::string_view query);
    

    std::string_view getHeader() const { return header; }
    std::string_view getOriginalSequence() const { return originalSequence; }
    std::string_view getBWTSequence() const { return bwtSequence; }

    bool getIsProtein() const { return isProtein; }
};


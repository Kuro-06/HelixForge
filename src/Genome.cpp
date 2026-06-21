#include "Genome.hpp"
#include <fstream>
#include <iostream>
#include <algorithm>
#include <cstdint>
#include <cctype>


Genome::Genome(bool proteinMode) : isProtein(proteinMode) {
    alphabetSize = isProtein ? 25 : 6; 
}

uint8_t Genome::getBaseIndex(char c) const {
    if (!isProtein) {
        switch (c) {
            case '$': return 0;
            case 'A': return 1;
            case 'C': return 2;
            case 'G': return 3;
            case 'T': return 4;
            default:  return 5; 
        }
    } else {
        switch(c) {
            case '$': return 0;
            case '*': return 1;
            case 'A': return 2;
            case 'B': return 3;
            case 'C': return 4;
            case 'D': return 5;
            case 'E': return 6;
            case 'F': return 7;
            case 'G': return 8;
            case 'H': return 9;
            case 'I': return 10;
            case 'K': return 11;
            case 'L': return 12;
            case 'M': return 13;
            case 'N': return 14;
            case 'P': return 15;
            case 'Q': return 16;
            case 'R': return 17;
            case 'S': return 18;
            case 'T': return 19;
            case 'V': return 20;
            case 'W': return 21;
            case 'X': return 22;
            case 'Y': return 23;
            case 'Z': return 24;
            default: return 22; // Map unknowns to 'X' (22)
        }
    }
}

// info extractor for fasta fiile
bool Genome::loadFromFasta(const std::string& filename) {
    std::ifstream file(filename);

    if (!file.is_open()) return false;

    std::string line;
    std::string tempSequence = "";

    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '>') {
            continue;
        }
        
        
        if (line.back() == '\r') {
            line.pop_back();
        }
        line.erase(std::remove_if(line.begin(), line.end(), ::isspace), line.end());
        tempSequence += line;
    }
    
    file.close();

    
    if (tempSequence.empty()) {
        return false;
    }
    

    
    originalSequence = tempSequence;
    return true; //load sucess
}

//struct forprefix doubling
struct Suffix {
    int index;     
    int rank[2];
};
void Genome::buildBWT() {
    if (originalSequence.empty()) return;
    

    std::string seq = originalSequence + "$";
    int n = seq.length();

    std::vector<Suffix> suffixes(n);
    
    for (int i = 0; i < n; i++) {
        suffixes[i].index = i;
        suffixes[i].rank[0] = seq[i];
        suffixes[i].rank[1] = ((i + 1) < n) ? (seq[i + 1]) : -1;   
    }

    //sorting based on 1st 2 chars
    auto compareSuffixes = [](const Suffix& a, const Suffix& b) {
        return (a.rank[0] == b.rank[0]) ? (a.rank[1] < b.rank[1]) : (a.rank[0] < b.rank[0]);
    };
    std::sort(suffixes.begin(), suffixes.end(), compareSuffixes);

    std::vector<int> ind(n);
    
    for (int k = 4; k < 2 * n; k = k * 2) {
        int rank = 0;
        int prev_rank = suffixes[0].rank[0];
        suffixes[0].rank[0] = rank;          
        ind[suffixes[0].index] = 0;          

        for (int i = 1; i < n; i++) {
            
            if (suffixes[i].rank[0] == prev_rank && suffixes[i].rank[1] == suffixes[i - 1].rank[1]) {
                prev_rank = suffixes[i].rank[0]; 
                suffixes[i].rank[0] = rank;      
            } else {
                prev_rank = suffixes[i].rank[0]; 
                suffixes[i].rank[0] = ++rank;   
            }
            ind[suffixes[i].index] = i;          
        }

        for (int i = 0; i < n; i++) {
            int nextindex = suffixes[i].index + k / 2;            
            suffixes[i].rank[1]= (nextindex < n) ? suffixes[ind[nextindex]].rank[0] : -1; 
        }

        std::sort(suffixes.begin(), suffixes.end(), compareSuffixes);
    }

    
    suffixArray.resize(n);
    
    for (int i = 0; i < n; i++) {
        suffixArray[i] = suffixes[i].index; 
    }

    // extract BWT
    bwtSequence.resize(n); 
    
    for (int i = 0; i < n; i++) {
        if (suffixArray[i] == 0) {       
            bwtSequence[i] = seq[n - 1]; 
        } else {
            bwtSequence[i] = seq[suffixArray[i] - 1]; 
        }
    }
    generateFMTables(); //for query searching
}

void Genome::generateFMTables() {
    int n = bwtSequence.length();

    //count occ
    counts.assign(alphabetSize, 0); 
    for (int i = 0; i < n; i++) {
        counts[getBaseIndex(bwtSequence[i])]++; 
    }
    
    //running totals
    int runningTotal = 0;
    for (int i = 0; i < alphabetSize; i++) {
        if (counts[i] > 0) {
            int temp = counts[i];     
            counts[i] = runningTotal; 
            runningTotal += temp;
        }
    }

    //occ matrix
    occ.assign(alphabetSize, std::vector<int>(n + 1, 0));
    for (int i = 0; i < n; i++) {
        for (int c = 0; c < alphabetSize; c++) {
            occ[c][i + 1] = occ[c][i]; 
        }
        occ[getBaseIndex(bwtSequence[i])][i + 1]++; 
    }
}

bool Genome::saveIndex(const std::string& filepath) {
    std::ofstream out(filepath, std::ios::binary);
    if (!out.is_open()) return false;

    uint8_t proteinFlag = isProtein ? 1 : 0;
    out.write(reinterpret_cast<const char*>(&proteinFlag), sizeof(proteinFlag));

    // Binary String writing Lamda
    auto writeString = [&out](const std::string& str) {
        uint64_t len = str.length();
        out.write(reinterpret_cast<const char*>(&len), sizeof(len)); // syntax : out.write(adress , no of bytes)
        // reinterpret cast : writes directly binary bits to file 
        out.write(str.data(), len);                                  
    };

    // Binary Vector<int> writing Lamda
    auto writeIntVector = [&out](const std::vector<int>& vec) {
        uint64_t len = vec.size();
        out.write(reinterpret_cast<const char*>(&len), sizeof(len)); 
        if (len > 0) {
            out.write(reinterpret_cast<const char*>(vec.data()), len * sizeof(int));
        }
    };

    // Save variables in exact order
    writeString(header);
    writeString(originalSequence);
    writeString(bwtSequence); 
    writeIntVector(suffixArray);

    out.close();
    return true;
}

bool Genome::loadIndex(const std::string& filepath) {
    std::ifstream in(filepath, std::ios::binary);
    if (!in.is_open()) return false;

    uint8_t proteinFlag = 0;
    in.read(reinterpret_cast<char*>(&proteinFlag), sizeof(proteinFlag));
    isProtein = (proteinFlag == 1);
    alphabetSize = isProtein ? 25 : 6;

   // Binary string reader (for the exact format we encoded the file)
    auto readString = [&in](std::string& str) {
        uint64_t len = 0;
        in.read(reinterpret_cast<char*>(&len), sizeof(len)); 
        str.resize(len);                                     
        in.read(&str[0], len);                               
    };

    //Binary vector reader (same as writing format)
    auto readIntVector = [&in](std::vector<int>& vec) {
        uint64_t len = 0;
        in.read(reinterpret_cast<char*>(&len), sizeof(len)); 
        vec.resize(len);
        if (len > 0) {
            in.read(reinterpret_cast<char*>(vec.data()), len * sizeof(int));
        }
    };

    //loading data form file into the class 
    readString(header);
    readString(originalSequence);
    readString(bwtSequence);
    readIntVector(suffixArray);


    generateFMTables(); //mush faster to create than again than storing and loading

    in.close();
    return true;
}

std::vector<int> Genome::search(const std::string_view query) {
    std::vector<int> matches; 
    if (bwtSequence.empty() || query.empty()) return matches;
    
    int top = 0;
    int bottom = bwtSequence.length() - 1;

    for (int i = query.length() - 1; i >= 0; i--) {
        
        uint8_t c_idx = getBaseIndex(query[i]); 
        
        
        if (c_idx >= alphabetSize) {
            return matches; 
        }

       
        top = counts[c_idx] + occ[c_idx][top]; 
        bottom = counts[c_idx] + occ[c_idx][bottom + 1] - 1; 
        
        if (top > bottom) return matches; 
    }

    for (int i = top; i <= bottom; i++) {
        matches.push_back(suffixArray[i]); 
    }
    
    return matches;
}
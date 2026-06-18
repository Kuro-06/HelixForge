#include "Genome.hpp"
#include <fstream>
#include <iostream>
#include <algorithm>

// info extractor for fasta fiile
void Genome::loadFromFasta(const std::string& filepath) {
    std::ifstream file(filepath);
    
    //check file path
    if (!file.is_open()) {
        std::cerr << "ERROR: Could not open FASTA file at " << filepath << std::endl;
        return; 
    }

    std::string line;
    bool isFirstHeaderFound = false;

    //resrting variables for the class
    header = "";
    originalSequence = "";

    std::cout << "Processing FASTA file: " << filepath << "..." << std::endl;

    while (std::getline(file, line)) {
        //remove empty char
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        // remove empty lines
        if (line.empty()) continue;

        // header existence check
        if (line[0] == '>') {
            if (isFirstHeaderFound) break;
            header = line.substr(1);
            isFirstHeaderFound = true;
        } 
        //else its a sequence file
        else {
            originalSequence += line;
        }
    }

    file.close();

    std::cout << "Loaded sequence: '" << header << "'" << std::endl;
    std::cout << "Total Base Pairs: " << originalSequence.length() << std::endl;
}

//helper struct for prefix doubling
struct Suffix {
    int index;     
    int rank[2];
};
void Genome::buildBWT() {
    if (originalSequence.empty()) {
        std::cout << "[ERROR] No sequence loaded. Cannot build BWT." << std::endl;
        return;
    }

    std::string seq = originalSequence + "$";
    int n = seq.length();

    std::vector<Suffix> suffixes(n);
    
    for (int i = 0; i < n; i++) {
        suffixes[i].index = i;
        suffixes[i].rank[0] = seq[i];
        suffixes[i].rank[1] = ((i + 1) < n) ? (seq[i + 1]) : -1;   
    }

    // 2. Sort based on the first two characters
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

    std::cout << "[SUCCESS] High-Speed BWT Compressed Engine Online!" << std::endl;
}

void Genome::generateFMTables() {
    int n = bwtSequence.length();

    //Count occurrences of each character
    counts.assign(256, 0); 
    for (int i = 0; i < n; i++) {
        counts[(unsigned char)bwtSequence[i]]++; 
    }


    //scan the counts to know which chars exist
    std::vector<unsigned char> activeAlphabet;
    for (int i = 0; i < 256; i++) {
        if (counts[i] > 0) {
            activeAlphabet.push_back((unsigned char)i);
        }
    }

    //running cumilative total
    int runningTotal = 0;
    for (int i = 0; i < 256; i++) {
        if (counts[i] > 0) {
            int temp = counts[i];     
            counts[i] = runningTotal; 
            runningTotal += temp;
        }
    }

    //occourance matrix
    occ.assign(256, std::vector<int>(n + 1, 0));
    for (int i = 0; i < n; i++) {
        // dosen't loop 256 times (at max 5-6 DNA / 21 Protein )
        for (unsigned char c : activeAlphabet) {
            occ[c][i + 1] = occ[c][i]; 
        }
        occ[(unsigned char)bwtSequence[i]][i + 1]++; 
    }
}

bool Genome::saveIndex(const std::string& filepath) {
    std::ofstream out(filepath, std::ios::binary);
    if (!out.is_open()) return false;

    // Binary String writing Lamda
    auto writeString = [&out](const std::string& str) {
        size_t len = str.length();
        out.write(reinterpret_cast<const char*>(&len), sizeof(len)); // syntax : out.write(adress , no of bytes)
        // reinterpret cast : writes directly binary bits to file 
        out.write(str.data(), len);                                  
    };

    // Binary Vector<int> writing Lamda
    auto writeIntVector = [&out](const std::vector<int>& vec) {
        size_t len = vec.size();
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

   // Binary string reader (for the exact format we encoded the file)
    auto readString = [&in](std::string& str) {
        size_t len = 0;
        in.read(reinterpret_cast<char*>(&len), sizeof(len)); 
        str.resize(len);                                     
        in.read(&str[0], len);                               
    };

    //Binary vector reader (same as writing format)
    auto readIntVector = [&in](std::vector<int>& vec) {
        size_t len = 0;
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

std::vector<int> Genome::search(const std::string& query) {
    std::vector<int> matches; 
    if (bwtSequence.empty() || query.empty()) return matches;
    
    int top = 0;
    int bottom = bwtSequence.length() - 1;

    for (int i = query.length() - 1; i >= 0; i--) {
        unsigned char c = query[i]; 
        top = counts[c] + occ[c][top]; 
        bottom = counts[c] + occ[c][bottom + 1] - 1; 
        if (top > bottom) return matches; 
    }

    for (int i = top; i <= bottom; i++) matches.push_back(suffixArray[i]); 
    return matches;
}
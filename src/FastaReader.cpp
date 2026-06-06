#include "FastaReader.hpp"
#include <fstream>
#include <iostream>

std::string extractSequenceFromFasta(std::string filePath) {
    // 1. Open the file
    std::ifstream file(filePath); 
    
    // 2. Professional Safety Check: Did the file actually open successfully?
    if (!file.is_open()) {
        std::cout << "Error: Could not open file at " << filePath << std::endl;
        return ""; // Return an empty string if we failed
    }

    std::string sequence = "";
    std::string currentLine;

    // 3. Read the file one line at a time until we reach the end
    while (std::getline(file, currentLine)) {
        
        // Skip any blank lines to be safe
        if (currentLine.empty()) {
            continue; 
        }

        // Your Logic: If the line starts with '>', it's a name, not DNA. Skip it!
        if (currentLine[0] == '>') {
            continue; 
        }

        // If it's not a blank line and doesn't start with '>', it must be our sequence data.
        // We append it to our master sequence string.
        sequence += currentLine;
    }

    // 4. Always close the file when you are done!
    file.close(); 
    
    return sequence;
}
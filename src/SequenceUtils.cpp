#include "SequenceUtils.hpp"
#include <algorithm>
#include <cctype>

namespace SequenceUtils {

    std::string detectSequenceType(const std::string& sequence) {
        int checkLength = std::min(100, (int)sequence.length());
        
        for (int i = 0; i < checkLength; i++) {
            char c = std::toupper(sequence[i]);
            // If we find letters exclusive to proteins
            if (c == 'M' || c == 'F' || c == 'P' || c == 'Z' || c == 'Q' || c == 'E') {
                return "protein";
            }
        }
        return "dna"; // Default fallback
    }

}
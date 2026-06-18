#ifndef SEQUENCE_UTILS_HPP
#define SEQUENCE_UTILS_HPP

#include <string>

namespace SequenceUtils {
    // Analyzes a sequence and returns "dna" or "protein"
    std::string detectSequenceType(const std::string& sequence);
}

#endif
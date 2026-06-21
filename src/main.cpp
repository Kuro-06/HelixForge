#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <memory>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include "AlignmentEngine.hpp"
#include "Genome.hpp"
#include "AlignmentWorkflow.hpp"


// state machine
struct HelixState {
    bool isProtein = false;
    int kSize = 15;
    std::string refName = "None";
    std::string queryName = "None";
    std::unique_ptr<Genome> refGenome = nullptr;
    std::unique_ptr<Genome> queryGenome = nullptr;
};

//utility
inline bool fileExists(const std::string& name) {
    std::ifstream f(name.c_str());
    return f.good();
}

std::pair<std::string, std::string> resolvePaths(const std::string& input, bool isProtein) {
    std::string fastaPath, idxPath;
    //OS path check
    if (input.find('/') != std::string::npos || input.find('\\') != std::string::npos || input.find('~') != std::string::npos) {
        fastaPath = input;
        size_t dotPos = fastaPath.find_last_of('.');
        if (dotPos != std::string::npos) {
            idxPath = fastaPath.substr(0, dotPos) + ".idx";
        } else {
            idxPath = fastaPath + ".idx";
        }
    } else { // Database lookup
        std::string folder = isProtein ? "database/protein/" : "database/dna/";
        fastaPath = folder + input + ".fasta";
        idxPath = folder + input + ".idx";
    }
    return {fastaPath, idxPath};
}

// ui/ux
class UIController {
public:
    static void printHeader() {
        std::cout << "\n=========================================================\n";
        std::cout << "      🧬 HELIXFORGE INTERACTIVE WORKBENCH v1.0\n";
        std::cout << "=========================================================\n";
        std::cout << "Type 'help' for commands or 'exit' to quit.\n";
    }

    static void printHelp() {
        std::cout << "\n--- HelixForge Commands ---\n";
        std::cout << "  mode <dna|protein>       : Switch scoring and sequence alphabets\n";
        std::cout << "  load ref <name_or_path>  : Load reference (e.g., 'load ref ecoli' or 'load ref ~/seq.fasta')\n";
        std::cout << "  load query <name_or_path>: Load target query\n";
        std::cout << "  set kmer <int>           : Change Seed Size (Default: DNA=15, Prot=4)\n";
        std::cout << "  status                   : View loaded memory status\n";
        std::cout << "  align                    : Execute Smith-Waterman Seed-and-Extend pipeline\n";
        std::cout << "  clear                    : Wipe terminal\n";
        std::cout << "  exit                     : Safely dump memory and quit\n";
    }

    static void printChunkedTraceback(const std::string& ref, const std::string& qry, int refStart, int qryStart) {
        int currentRefPos = refStart;
        int currentQryPos = qryStart;

        for (size_t i = 0; i < ref.length(); i += 60) {
            std::string rChunk = ref.substr(i, 60);
            std::string qChunk = qry.substr(i, 60);
            std::string matchLine = "";
            
            int rAdv = 0, qAdv = 0;
            for(size_t j = 0; j < rChunk.length(); j++){
                if(rChunk[j] == qChunk[j]) matchLine += "|";
                else matchLine += " ";
                if(rChunk[j] != '-') rAdv++;
                if(qChunk[j] != '-') qAdv++;
            }

            std::cout << "Ref: [" << currentRefPos << "]\t" << rChunk << "  [" << currentRefPos + rAdv - 1 << "]\n";
            std::cout << "     \t\t" << matchLine << "\n";
            std::cout << "Qry: [" << currentQryPos << "]\t" << qChunk << "  [" << currentQryPos + qAdv - 1 << "]\n\n";
            
            currentRefPos += rAdv;
            currentQryPos += qAdv;
        }
    }

    static void printReport(const HelixState& state, const AlignmentWorkflow::ClusterResult& cluster, int score, double evalue, const TracebackResult& tb,int paddedStart) {
        std::cout << "=========================================================\n";
        std::cout << "[🧬] HELIXFORGE ALIGNMENT REPORT\n";
        std::cout << "=========================================================\n";
        std::cout << "Sequence Mode : " << (state.isProtein ? "Protein (BLOSUM62)" : "DNA (+1/-1)") << "\n";
        std::cout << "Reference     : " << state.refName << " (" << state.refGenome->getOriginalSequence().length() << " bp)\n";
        std::cout << "Query         : " << state.queryName << " (" << state.queryGenome->getOriginalSequence().length() << " bp)\n";
        std::cout << "Algorithm     : Smith-Waterman (Local)\n";
        std::cout << "---------------------------------------------------------\n";
        std::cout << " SEED & EXTEND MAPPING (k=" << state.kSize << ")\n";
        std::cout << "---------------------------------------------------------\n";
        std::cout << "Best Cluster  : [Ref Bounds: " << cluster.start << " - " << cluster.end << "]\n";
        std::cout << "Seed Density  : " << cluster.maxDensity << " overlapping segments (from " << cluster.totalHits << " hits)\n";
        std::cout << "---------------------------------------------------------\n";
        std::cout << " STATISTICAL SIGNIFICANCE\n";
        std::cout << "---------------------------------------------------------\n";
        std::cout << "Raw Score (S) : " << score << "\n";
        std::cout << "E-Value       : " << std::scientific << evalue << std::defaultfloat << " \n";
        std::cout << "---------------------------------------------------------\n";
        std::cout << " ALIGNMENT TRACEBACK\n";
        std::cout << "---------------------------------------------------------\n";

        int trueRefStart = paddedStart + tb.refStartOffset + 1;
        int trueQryStart = tb.queryStartOffset + 1;
        
        printChunkedTraceback(tb.alignedRef, tb.alignedQuery, trueRefStart, trueQryStart);
        
        std::cout << "=========================================================\n";
    }
};

//repl shell
int main() {
    HelixState state;
    UIController::printHeader();

    std::string input;
    while (true) {
        std::cout << "\nHelixForge [" << (state.isProtein ? "Protein" : "DNA") << "] > ";
        if (!std::getline(std::cin, input)) break;
        
        if (input.empty()) continue;

        std::vector<std::string> tokens;
        std::string token;
        std::istringstream tokenStream(input);
        while (tokenStream >> token) tokens.push_back(token);

        std::string cmd = tokens[0];

        //utility (pipeline)
        if (cmd == "exit" || cmd == "quit") {
            std::cout << "[SYSTEM] Dumping RAM... Goodbye.\n";
            break;
        } 
        else if (cmd == "clear") {
            #if defined(_WIN32)
                system("cls");
            #else
                system("clear");
            #endif
            UIController::printHeader();
        }
        else if (cmd == "help") {
            UIController::printHelp();
        }
        else if (cmd == "status") {
            std::cout << "\n[SYSTEM MEMORY STATUS]\n";
            std::cout << "Mode  : " << (state.isProtein ? "Protein" : "DNA") << "\n";
            std::cout << "K-mer : " << state.kSize << "\n";
            std::cout << "Ref   : " << state.refName << (state.refGenome ? " [LOCKED IN RAM]" : " [EMPTY]") << "\n";
            std::cout << "Query : " << state.queryName << (state.queryGenome ? " [LOCKED IN RAM]" : " [EMPTY]") << "\n";
        }
        
        else if (cmd == "set" && tokens.size() >= 3 && tokens[1] == "kmer") {
            state.kSize = std::stoi(tokens[2]);
            std::cout << "[SYSTEM] Seed size set to " << state.kSize << ".\n";
        }
        else if (cmd == "mode" && tokens.size() >= 2) {
            if (tokens[1] == "protein") {
                state.isProtein = true;
                state.kSize = 4;
            } else if (tokens[1] == "dna") {
                state.isProtein = false;
                state.kSize = 15;
            } else {
                std::cout << "[ERROR] Invalid mode. Use 'dna' or 'protein'.\n";
                continue;
            }
            state.refGenome.reset();
            state.queryGenome.reset();
            state.refName = "None";
            state.queryName = "None";
            std::cout << "[SYSTEM] Engine switched to " << (state.isProtein ? "Protein" : "DNA") << " Mode. Memory cleared.\n";
        }
        
        else if (cmd == "load" && tokens.size() >= 3) {
            std::string type = tokens[1];
            std::string pathInput = tokens[2];

            auto paths = resolvePaths(pathInput, state.isProtein);
            
            if (type == "ref") {
                state.refGenome = std::make_unique<Genome>(state.isProtein);
                
                if (fileExists(paths.second)) {
                    std::cout << "[SYSTEM] Found .idx. Instant RAM load initiated...\n";
                    state.refGenome->loadIndex(paths.second);
                } else {
                    std::cout << "[SYSTEM] No .idx found. Building FM-Index (This may take a moment)...\n";
                    if (!state.refGenome->loadFromFasta(paths.first)) {
                        std::cout << "[ERROR] Failed to load " << paths.first << "\n";
                        state.refGenome.reset();
                        continue;
                    }
                    state.refGenome->buildBWT();
                    state.refGenome->saveIndex(paths.second);
                }
                state.refName = pathInput;
                std::cout << "[SUCCESS] Reference locked in RAM.\n";
            } 
            else if (type == "query") {
                state.queryGenome = std::make_unique<Genome>(state.isProtein);
                if (state.queryGenome->loadFromFasta(paths.first)) {
                    state.queryName = pathInput;
                    std::cout << "[SUCCESS] Query locked in RAM.\n";
                } else {
                    std::cout << "[ERROR] Failed to load Query.\n";
                    state.queryGenome.reset();
                }
            }
        }
        
        else if (cmd == "align") {
            if (!state.refGenome || !state.queryGenome) {
                std::cout << "[ERROR] Must load both reference and query before aligning.\n";
                continue;
            }

            
            std::string_view refSeq = state.refGenome->getOriginalSequence();
            std::string_view qrySeq = state.queryGenome->getOriginalSequence();

            if (qrySeq.length() > refSeq.length()) {
                std::cout << "[WARNING] Query sequence is larger than Reference.\n";
                std::cout << "This breaks the FM-Index bounds. Swap them? (y/n): ";
                std::string ans;
                std::getline(std::cin, ans);
                if (ans == "y" || ans == "Y") {
                    std::swap(state.refGenome, state.queryGenome);
                    std::swap(state.refName, state.queryName);
                    refSeq = state.refGenome->getOriginalSequence();
                    qrySeq = state.queryGenome->getOriginalSequence();
                    std::cout << "[SUCCESS] Sequences swapped.\n";
                } else {
                    std::cout << "[ABORTED] Alignment cancelled.\n";
                    continue;
                }
            }

            auto cluster = AlignmentWorkflow::findBestCluster(*state.refGenome, qrySeq, state.kSize);
            
            if (!cluster.success) {
                std::cout << "[RESULT] No viable matches found. Try lowering k-mer size.\n";
                continue;
            }

            int padding = state.isProtein ? 20 : 50;
            int paddedStart = std::max(0, cluster.start - padding);
            int paddedLength = std::min(static_cast<int>(refSeq.length() - paddedStart), (cluster.end - cluster.start) + (padding * 2));
            std::string_view targetWindow = refSeq.substr(paddedStart, paddedLength);

            AlignmentEngine aligner(targetWindow, qrySeq, state.isProtein);
            int rawScore = aligner.Waterman();
            double evalue = aligner.EValue(rawScore);
            auto tb = aligner.TraceBack();

            UIController::printReport(state, cluster, rawScore, evalue, tb,paddedStart);
        } 
        else {
            std::cout << "[ERROR] Unknown command. Type 'help'.\n";
        }
    }

    return 0;
}
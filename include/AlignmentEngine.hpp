#pragma once
#include <string>
#include <vector>
#include <string_view>

struct TracebackResult {
    std::string alignedRef;
    std::string alignedQuery;
    int refStartOffset;
    int queryStartOffset;
};


class AlignmentEngine {
private:
    std::string_view seqA;
    std::string_view seqB;
    bool isProtein;

    int match;
    int mismatch;
    int gapOpen;
    int gapExtend;

    double lambda;
    double k_stat;

    std::vector<std::vector<int>> matrixM;
    std::vector<std::vector<int>> matrixX;
    std::vector<std::vector<int>> matrixY;
    std::vector<std::vector<char>> dirMatrix;

    bool isLocal;
    int starti;
    int startj;

    int getBlosumScore(char a, char b) const;

public:
    // constructor
    AlignmentEngine(std::string_view inputSeqA, std::string_view inputSeqB, bool proteinMode = false);

    int Needleman();
    int Waterman();
    TracebackResult TraceBack();
//    void scorePara(int a,int b,int c,int d);  (future feature)
    double EValue(int rawScore) const;
};
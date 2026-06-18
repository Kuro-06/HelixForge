#pragma once
#include <string>
#include <vector>

class AlignmentEngine {
private:
    std::string seqA;
    std::string seqB;

    int match;
    int mismatch;
    int gapOpen;
    int gapExtend;

    std::vector<std::vector<int>> matrixM;
    std::vector<std::vector<int>> matrixX;
    std::vector<std::vector<int>> matrixY;
    std::vector<std::vector<char>> dirMatrix;

    bool isLocal;
    int starti;
    int startj;

public:
    // constructor
    AlignmentEngine(std::string inputSeqA, std::string inputSeqB);

    // Action Buttons
    int Needleman();
    int Waterman();
    void TraceBack();
    void scorePara(int a,int b,int c,int d);
};
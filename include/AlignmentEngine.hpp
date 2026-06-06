#pragma once
#include <string>
#include <vector>

// THIS is the only place the word "class" goes!
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
    // The Constructor (Boots up the machine)
    AlignmentEngine(std::string inputSeqA, std::string inputSeqB);

    // The Action Buttons
    int Needleman();
    int Waterman();
    void TraceBack();
    void scorePara(int a,int b,int c,int d);
};
#include "AlignmentEngine.hpp"
#include <iostream> 
#include <string>
#include <vector>   
#include <algorithm> // Required for std::max()

AlignmentEngine::AlignmentEngine(std::string inputSeqA, std::string inputSeqB) {
    seqA = inputSeqA;
    seqB = inputSeqB;
    match = 1;
    mismatch = -1;
    gapOpen = -10;
    gapExtend = -1;
}

void AlignmentEngine::scorePara(int customMatch, int customMismatch, int customGapOpen,int customGapExtend) {
    match = customMatch;
    mismatch = customMismatch;
    gapOpen = customGapOpen;
    gapExtend = customGapExtend;
}

int AlignmentEngine::Needleman() {

    isLocal=false;

    int row = seqA.size()+1;
    int col = seqB.size()+1;

    int NEG_INF = -1000000;

    matrixM.assign(row, std::vector<int>(col, NEG_INF));
    matrixX.assign(row, std::vector<int>(col, NEG_INF));
    matrixY.assign(row, std::vector<int>(col, NEG_INF));
    dirMatrix.assign(row, std::vector<char>(col, ' '));

    matrixM[0][0] = 0;

    for(int i=1;i<col;i++){
        matrixY[0][i]= gapOpen + i*gapExtend;
        dirMatrix[0][i]='L';
    }
    for(int i=1;i<row;i++){
        matrixX[i][0]= gapOpen + i*gapExtend;
        dirMatrix[i][0]='U';
    }

    for(int i=1;i<row;i++){
        for(int j=1;j<col;j++){
            int dig=0;
            seqA[i-1]==seqB[j-1]?dig=match:dig=mismatch;

            matrixM[i][j] = std::max({
                matrixM[i-1][j-1], 
                matrixX[i-1][j-1], 
                matrixY[i-1][j-1]
            }) + dig;

            matrixX[i][j] = std::max({
                matrixM[i-1][j] + gapOpen + gapExtend,
                matrixX[i-1][j] + gapExtend
            });

            matrixY[i][j] = std::max({
                matrixM[i][j-1] + gapOpen + gapExtend,
                matrixY[i][j-1] + gapExtend
            });

            int best = std::max({matrixM[i][j], matrixX[i][j], matrixY[i][j]});
            if (best == matrixM[i][j]) {
                dirMatrix[i][j] = 'D';
            } else if (best == matrixX[i][j]) {
                dirMatrix[i][j] = 'U';
            } else {
                dirMatrix[i][j] = 'L';
            }
        }
    }

    return std::max({matrixM[row-1][col-1],matrixX[row-1][col-1], matrixY[row-1][col-1]});
}

int AlignmentEngine::Waterman() {

    isLocal=true;

    int row = seqA.size()+1;
    int col = seqB.size()+1;

    matrixM.assign(row, std::vector<int>(col, 0));
    matrixX.assign(row, std::vector<int>(col, 0));
    matrixY.assign(row, std::vector<int>(col, 0));
    dirMatrix.assign(row, std::vector<char>(col, ' '));

    matrixM[0][0] = 0;
    int maxScore=0;
    int maxi=0;
    int maxj=0;

    for(int i=1;i<row;i++){
        for(int j=1;j<col;j++){
            int dig=0;
            seqA[i-1]==seqB[j-1]?dig=match:dig=mismatch;

            matrixM[i][j] = std::max({
                matrixM[i-1][j-1]+ dig, 
                matrixX[i-1][j-1]+ dig, 
                matrixY[i-1][j-1]+ dig,
                0
            });

            matrixX[i][j] = std::max({
                matrixM[i-1][j] + gapOpen + gapExtend,
                matrixX[i-1][j] + gapExtend,
                0
            });

            matrixY[i][j] = std::max({
                matrixM[i][j-1] + gapOpen + gapExtend,
                matrixY[i][j-1] + gapExtend,
                0
            });

            int best = std::max({matrixM[i][j], matrixX[i][j], matrixY[i][j]});
            if(best>maxScore){
                maxScore=best;
                maxi=i;
                maxj=j;
            }

            if (best == matrixM[i][j]) {
                dirMatrix[i][j] = 'D';
            } else if (best == matrixX[i][j]) {
                dirMatrix[i][j] = 'U';
            } else {
                dirMatrix[i][j] = 'L';
            }
        }
    }

    starti=maxi;
    startj=maxj;
    return maxScore;
}

void AlignmentEngine::TraceBack() {
    std::string alignedA = "";
    std::string alignedB = "";

    int i = isLocal ? starti : seqA.size();
    int j = isLocal ? startj : seqB.size();

    while (i > 0 || j > 0) {

        if (isLocal && (matrixM[i][j] <= 0 && matrixX[i][j] <= 0 && matrixY[i][j] <= 0)) {
            break;
        }

        char direction = dirMatrix[i][j];

        if (direction == 'D') {
            alignedA += seqA[i - 1];
            alignedB += seqB[j - 1];
            i--; j--;
        } 
        else if (direction == 'U') {
            alignedA += seqA[i - 1];
            alignedB += '-';
            i--;
        } 
        else if (direction == 'L') {
            alignedA += '-';
            alignedB += seqB[j - 1];
            j--;
        }
        else {
            break;
        }
    }

    std::reverse(alignedA.begin(), alignedA.end());
    std::reverse(alignedB.begin(), alignedB.end());

    std::cout << "Aligned Sequence A: " << alignedA << std::endl;
    std::cout << "Aligned Sequence B: " << alignedB << std::endl;
}
#include "AlignmentEngine.hpp"
#include <iostream> 
#include <string>
#include <vector>   
#include <algorithm>
#include <cmath> // Required for std::max()


//blosum62
namespace {
    
    int getAAIndex(char c) {
        switch(c) {
            case 'A': return 0;  case 'R': return 1;  case 'N': return 2;  case 'D': return 3;
            case 'C': return 4;  case 'Q': return 5;  case 'E': return 6;  case 'G': return 7;
            case 'H': return 8;  case 'I': return 9;  case 'L': return 10; case 'K': return 11;
            case 'M': return 12; case 'F': return 13; case 'P': return 14; case 'S': return 15;
            case 'T': return 16; case 'W': return 17; case 'Y': return 18; case 'V': return 19;
            case 'B': return 20; case 'Z': return 21; case 'X': return 22; case '*': return 23;
            default: return 22; 
        }
    }

   
    const int BLOSUM62[24][24] = {
     // A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X 
/*A*/ { 4,-1,-2,-2, 0,-1,-1, 0,-2,-1,-1,-1,-1,-2,-1, 1, 0,-3,-2, 0,-2,-1, 0,-4},
/*R*/ {-1, 5, 0,-2,-3, 1, 0,-2, 0,-3,-2, 2,-1,-3,-2,-1,-1,-3,-2,-3,-1, 0,-1,-4},
/*N*/ {-2, 0, 6, 1,-3, 0, 0, 0, 1,-3,-3, 0,-2,-3,-2, 1, 0,-4,-2,-3, 3, 0,-1,-4},
/*D*/ {-2,-2, 1, 6,-3, 0, 2,-1,-1,-3,-4,-1,-3,-3,-1, 0,-1,-4,-3,-3, 4, 1,-1,-4},
/*C*/ { 0,-3,-3,-3, 9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1,-3,-3,-2,-4},
/*Q*/ {-1, 1, 0, 0,-3, 5, 2,-2, 0,-3,-2, 1, 0,-3,-1, 0,-1,-2,-1,-2, 0, 3,-1,-4},
/*E*/ {-1, 0, 0, 2,-4, 2, 5,-2, 0,-3,-3, 1,-2,-3,-1, 0,-1,-3,-2,-2, 1, 4,-1,-4},
/*G*/ { 0,-2, 0,-1,-3,-2,-2, 6,-2,-4,-4,-2,-3,-3,-2, 0,-2,-2,-3,-3,-1,-2,-1,-4},
/*H*/ {-2, 0, 1,-1,-3, 0, 0,-2, 8,-3,-3,-1,-2,-1,-2,-1,-2,-2, 2,-3, 0, 0,-1,-4},
/*I*/ {-1,-3,-3,-3,-1,-3,-3,-4,-3, 4, 2,-3, 1, 0,-3,-2,-1,-3,-1, 3,-3,-3,-1,-4},
/*L*/ {-1,-2,-3,-4,-1,-2,-3,-4,-3, 2, 4,-2, 2, 0,-3,-2,-1,-2,-1, 1,-4,-3,-1,-4},
/*K*/ {-1, 2, 0,-1,-3, 1, 1,-2,-1,-3,-2, 5,-1,-3,-1, 0,-1,-3,-2,-2, 0, 1,-1,-4},
/*M*/ {-1,-1,-2,-3,-1, 0,-2,-3,-2, 1, 2,-1, 5, 0,-2,-1,-1,-1,-1, 1,-3,-1,-1,-4},
/*F*/ {-2,-3,-3,-3,-2,-3,-3,-3,-1, 0, 0,-3, 0, 6,-4,-2,-2, 1, 3,-1,-3,-3,-1,-4},
/*P*/ {-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4, 7,-1,-1,-4,-3,-2,-2,-1,-2,-4},
/*S*/ { 1,-1, 1, 0,-1, 0, 0, 0,-1,-2,-2, 0,-1,-2,-1, 4, 1,-3,-2,-2, 0, 0, 0,-4},
/*T*/ { 0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 1, 5,-2,-2, 0,-1,-1, 0,-4},
/*W*/ {-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1, 1,-4,-3,-2,11, 2,-3,-4,-3,-2,-4},
/*Y*/ {-2,-2,-2,-3,-2,-1,-2,-3, 2,-1,-1,-2,-1, 3,-3,-2,-2, 2, 7,-1,-3,-2,-1,-4},
/*V*/ { 0,-3,-3,-3,-1,-2,-2,-3,-3, 3, 1,-2, 1,-1,-2,-2, 0,-3,-1, 4,-3,-2,-1,-4},
/*B*/ {-2,-1, 3, 4,-3, 0, 1,-1, 0,-3,-4, 0,-3,-3,-2, 0,-1,-4,-3,-3, 4, 1,-1,-4},
/*Z*/ {-1, 0, 0, 1,-3, 3, 4,-2, 0,-3,-3, 1,-1,-3,-1, 0,-1,-3,-2,-2, 1, 4,-1,-4},
/*X*/ { 0,-1,-1,-1,-2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-2, 0, 0,-2,-1,-1,-1,-1,-1,-4},
/***/ {-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4, 1}
    };
}



AlignmentEngine::AlignmentEngine(std::string_view inputSeqA, std::string_view inputSeqB, bool proteinMode) {
    seqA = inputSeqA;
    seqB = inputSeqB;
    isProtein = proteinMode;
    match = 1;
    mismatch = -1;
    gapOpen = -10;
    gapExtend = -1;

    if (isProtein) {
       
        lambda = 0.267;
        k_stat = 0.041;
    } else {
        
        lambda = 1.37;
        k_stat = 0.711;
    }

}

int AlignmentEngine::getBlosumScore(char a, char b) const {
    return BLOSUM62[getAAIndex(a)][getAAIndex(b)];
}

// void AlignmentEngine::scorePara(int customMatch, int customMismatch, int customGapOpen,int customGapExtend) {
//     match = customMatch;
//     mismatch = customMismatch;
//     gapOpen = customGapOpen;
//     gapExtend = customGapExtend;
// }  Future feature

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
            dig = isProtein ? getBlosumScore(seqA[i-1], seqB[j-1]) : (seqA[i-1] == seqB[j-1] ? match : mismatch);

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
            dig = isProtein ? getBlosumScore(seqA[i-1], seqB[j-1]) : (seqA[i-1] == seqB[j-1] ? match : mismatch);

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

TracebackResult AlignmentEngine::TraceBack() {
    TracebackResult result;
    result.alignedRef = "";
    result.alignedQuery = "";

    int i = isLocal ? starti : seqA.size();
    int j = isLocal ? startj : seqB.size();

    while (i > 0 || j > 0) {

        if (isLocal && (matrixM[i][j] <= 0 && matrixX[i][j] <= 0 && matrixY[i][j] <= 0)) {
            break;
        }

        char direction = dirMatrix[i][j];

        if (direction == 'D') {
            result.alignedRef += seqA[i - 1];
            result.alignedQuery += seqB[j - 1];
            i--; j--;
        } 
        else if (direction == 'U') {
            result.alignedRef += seqA[i - 1];
            result.alignedQuery += '-';
            i--;
        } 
        else if (direction == 'L') {
            result.alignedRef += '-';
            result.alignedQuery += seqB[j - 1];
            j--;
        }
        else {
            break;
        }
    }

    result.refStartOffset = i;
    result.queryStartOffset = j;

    std::reverse(result.alignedRef.begin(), result.alignedRef.end());
    std::reverse(result.alignedQuery.begin(), result.alignedQuery.end());

    return result;
}

double AlignmentEngine::EValue(int rawScore) const {
    
    double m = static_cast<double>(seqA.length());
    double n = static_cast<double>(seqB.length());
    
    // E = K * m * n * e^(-lambda * S)
    double eValue = k_stat * m * n * std::exp(-lambda * static_cast<double>(rawScore));
    
    return eValue;
}
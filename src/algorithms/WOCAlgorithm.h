// WOCAlgorithm.h

#ifndef ESTPGENETICALGORITHMWOC_WOCALGORITHM_H
#define ESTPGENETICALGORITHMWOC_WOCALGORITHM_H

#include "globalData.h"
#include <vector>
#include <unordered_set>
#include <unordered_map>

// Enumeration for selecting the strategy
enum class Strategy {
    MSTBased,            // Minimum Spanning Tree Based Algorithm
    SPH,                 // Shortest Path Heuristic Algorithm
    Probabilistic,       // Probabilistic Edge Selection Algorithm
    IterativeAddition,   // Iterative Edge Addition Algorithm
    Threshold            // Edge Frequency Thresholding Algorithm
};

class WOCAlgorithm {
public:
    // Generates a Steiner Tree based on the selected strategy
    SteinerTree generateSteinerTree(Strategy strategy);

private:
    // Strategy-specific methods
    SteinerTree generateMSTBasedSteinerTree();
    SteinerTree generateSPHSteinerTree();
    SteinerTree generateProbabilisticSteinerTree();
    SteinerTree generateIterativeAdditionSteinerTree();
    SteinerTree generateThresholdSteinerTree();

    // Additional helper functions for strategies
    std::vector<int> findValidPath(const std::unordered_map<int, std::vector<int>>& graph, int start);
};

#endif //ESTPGENETICALGORITHMWOC_WOCALGORITHM_H

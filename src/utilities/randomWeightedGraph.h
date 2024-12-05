// randomWeightedGraph.h
// Created by LuisF

#ifndef ESTPGENETICALGORITHMWOC_RANDOMWEIGHTEDGRAPH_H
#define ESTPGENETICALGORITHMWOC_RANDOMWEIGHTEDGRAPH_H

#include <vector>
#include "globalData.h" // Include the global variables

class randomWeightedGraph {
public:
    randomWeightedGraph();
    void generateGraph();

private:
    void getUserInput();
    void determineTerminalCount();
    void generateVertices();
    void generateEdges();
    bool edgesIntersect(int a1, int a2, int b1, int b2);

};

#endif // ESTPGENETICALGORITHMWOC_RANDOMWEIGHTEDGRAPH_H

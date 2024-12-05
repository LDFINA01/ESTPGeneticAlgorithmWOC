//
// Created by LuisF on 11/14/2024.
//

#ifndef ESTPGENETICALGORITHMWOC_GLOBALDATA_H
#define ESTPGENETICALGORITHMWOC_GLOBALDATA_H


#include <vector>
#include "unordered_set"

struct EdgeInfo {
    int vertexA;
    int vertexB;
    double weight;
    int frequency;
};

// struct to store Steiner tree information
struct SteinerTree {
    std::vector<int> path;       // The sequence of vertices in the tree
    double totalDistance;        // Total distance of the tree
};


extern int n; // Number of vertices
extern int numTerminals;
extern std::vector<std::vector<double>> vertices; // [vertex number, x, y, terminal status]
extern std::vector<EdgeInfo> edges; // List of all edges with weights
extern std::unordered_set<int> terminalVertices; // Indices of terminal vertices
extern std::vector<double> avgDistanceOverGenerations; // The avg distance over a lot of time fr

#endif //ESTPGENETICALGORITHMWOC_GLOBALDATA_H

//
// Created by LuisF on 11/14/2024.
//

#ifndef ESTPGENETICALGORITHMWOC_GLOBALDATA_H
#define ESTPGENETICALGORITHMWOC_GLOBALDATA_H


#include <vector>
#include "unordered_set"

struct Node {
    double x;
    double y;
    bool isTerminal;
    
    Node() : x(0.0), y(0.0), isTerminal(false) {}
    Node(double x, double y, bool isTerminal = false) : x(x), y(y), isTerminal(isTerminal) {}
};

struct EdgeInfo {
    int vertexA;
    int vertexB;
    double weight;
    int frequency;
};

// struct to store Steiner tree information
struct SteinerTree {
    std::vector<int> path;       // The sequence of nodes in the tree
    double totalDistance;        // Total distance of the tree
};


extern int n; // Number of nodes
extern int numTerminals;
extern std::vector<Node> nodes; // Vector of nodes with x, y coordinates and terminal status
extern std::vector<EdgeInfo> edges; // List of all edges with weights
extern std::unordered_set<int> terminalVertices; // Indices of terminal nodes
extern std::vector<double> avgDistanceOverGenerations; // The avg distance over a lot of time fr

#endif //ESTPGENETICALGORITHMWOC_GLOBALDATA_H

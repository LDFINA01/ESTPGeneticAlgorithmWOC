// randomWeightedGraph.cpp
// Created by LuisF on
#include "randomWeightedGraph.h"
#include "globalData.h" // Include the global variables
#include <iostream>
#include <cmath>
#include <random>
#include <algorithm>
#include <set>

randomWeightedGraph::randomWeightedGraph() {
}

void randomWeightedGraph::generateGraph() {
    //getUserInput();
    determineTerminalCount();
    generateVertices();
    generateEdges();
}

void randomWeightedGraph::getUserInput() {
    std::cout << "How many vertices would you like? ";
    std::cin >> n;
    if (n <= 0) {
        std::cerr << "Wrong input." << std::endl;
        exit(1);
    }
}

void randomWeightedGraph::determineTerminalCount() {
    std::random_device rd;
    std::mt19937 gen(rd());
    double minPercentage, maxPercentage;



    minPercentage = 0.10;
    maxPercentage = 0.30;

    std::uniform_real_distribution<> dis(minPercentage, maxPercentage);
    double percentage = dis(gen);
    numTerminals = static_cast<int>(std::round(n * percentage));
}

void randomWeightedGraph::generateVertices() {
    vertices.resize(n, std::vector<double>(4)); // [vertex number, x, y, terminal status]
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> coordDis(0.0, 100.0);
    std::vector<int> indices(n);
    for (int i = 0; i < n; ++i) indices[i] = i;
    std::shuffle(indices.begin(), indices.end(), gen);

    terminalVertices.clear();
    for (int i = 0; i < numTerminals; ++i) {
        terminalVertices.insert(indices[i]);
    }

    for (int i = 0; i < n; ++i) {
        vertices[i][0] = i; // vertex number
        vertices[i][1] = coordDis(gen); // x coordinate
        vertices[i][2] = coordDis(gen); // y coordinate
        vertices[i][3] = 0; // default to non-terminal
    }
    for (int idx : terminalVertices) {
        vertices[idx][3] = 1; // mark as terminal
    }
}

// Helper function to check if two edges intersect
bool randomWeightedGraph::edgesIntersect(int a1, int a2, int b1, int b2) {
    double x1 = vertices[a1][1], y1 = vertices[a1][2];
    double x2 = vertices[a2][1], y2 = vertices[a2][2];
    double x3 = vertices[b1][1], y3 = vertices[b1][2];
    double x4 = vertices[b2][1], y4 = vertices[b2][2];

    // Exclude shared vertices
    if (a1 == b1 || a1 == b2 || a2 == b1 || a2 == b2)
        return false;

    // Compute direction of the lines
    double denom = (y4 - y3)*(x2 - x1) - (x4 - x3)*(y2 - y1);
    if (denom == 0) return false; // Lines are parallel

    double ua = ((x4 - x3)*(y1 - y3) - (y4 - y3)*(x1 - x3)) / denom;
    double ub = ((x2 - x1)*(y1 - y3) - (y2 - y1)*(x1 - x3)) / denom;

    return (ua > 0 && ua < 1 && ub > 0 && ub < 1);
}

void randomWeightedGraph::generateEdges() {
    double additionalEdgePercentage = 2; // Adjust this value to change the percentage of additional edges

    edges.clear();

    std::random_device rd;
    std::mt19937 gen(rd());

    // Initialize Union-Find structure
    std::vector<int> parent(n);
    std::iota(parent.begin(), parent.end(), 0); // parent[i] = i

    auto find = [&](int u) {
        while (parent[u] != u)
            u = parent[u];
        return u;
    };

    auto unite = [&](int u, int v) {
        parent[find(u)] = find(v);
    };

    // Set to store existing edges with ordered pairs
    std::set<std::pair<int, int>> existingEdges;

    // Generate all possible edges without intersections
    std::vector<EdgeInfo> possibleEdges;

    for (int P = 0; P < n; ++P) {
        // Coordinates of P
        double Px = vertices[P][1];
        double Py = vertices[P][2];

        // Compute polar coordinates of other points relative to P
        std::vector<std::tuple<double, double, int>> polarCoords;
        for (int i = 0; i < n; ++i) {
            if (i == P) continue;
            double dx = vertices[i][1] - Px;
            double dy = vertices[i][2] - Py;
            double r = std::sqrt(dx * dx + dy * dy);
            double theta = std::atan2(dy, dx);
            if (theta < 0) theta += 2 * M_PI;
            polarCoords.emplace_back(theta, r, i);
        }

        // Sort points by angle
        std::sort(polarCoords.begin(), polarCoords.end());

        // For each point, check visibility
        for (const auto& [thetaC, rC, C] : polarCoords) {
            bool blocked = false;
            // Check against existing edges
            for (const auto& edge : edges) {
                int A = edge.vertexA;
                int B = edge.vertexB;
                if (A == P || B == P || A == C || B == C)
                    continue;

                if (edgesIntersect(P, C, A, B)) {
                    blocked = true;
                    break;
                }
            }
            if (!blocked) {
                double dx = vertices[P][1] - vertices[C][1];
                double dy = vertices[P][2] - vertices[C][2];
                double dist = std::sqrt(dx * dx + dy * dy);
                possibleEdges.push_back({P, C, dist});
            }
        }
    }

    // Sort possible edges by weight (optional)
    std::sort(possibleEdges.begin(), possibleEdges.end(), [](const EdgeInfo& a, const EdgeInfo& b) {
        return a.weight < b.weight;
    });

    // Kruskal's algorithm to ensure connectivity
    for (const auto& edge : possibleEdges) {
        int u = edge.vertexA;
        int v = edge.vertexB;
        if (find(u) != find(v)) {
            // Check if edge intersects existing edges
            bool intersects = false;
            for (const auto& existingEdge : edges) {
                if (edgesIntersect(u, v, existingEdge.vertexA, existingEdge.vertexB)) {
                    intersects = true;
                    break;
                }
            }
            if (!intersects) {
                edges.push_back(edge);
                unite(u, v);
                // Add edge to existingEdges set
                int a = u, b = v;
                if (a > b) std::swap(a, b);
                existingEdges.insert({a, b});
            }
        }
    }

    // Check if all vertices are connected
    int root = find(0);
    for (int i = 1; i < n; ++i) {
        if (find(i) != root) {
            // Try to connect the disconnected component
            for (int j = 0; j < n; ++j) {
                if (find(j) == root) {
                    double dx = vertices[i][1] - vertices[j][1];
                    double dy = vertices[i][2] - vertices[j][2];
                    double dist = std::sqrt(dx * dx + dy * dy);
                    EdgeInfo newEdge = {i, j, dist};

                    // Check if edge intersects existing edges
                    bool intersects = false;
                    for (const auto& existingEdge : edges) {
                        if (edgesIntersect(i, j, existingEdge.vertexA, existingEdge.vertexB)) {
                            intersects = true;
                            break;
                        }
                    }
                    if (!intersects) {
                        edges.push_back(newEdge);
                        unite(i, j);
                        // Add edge to existingEdges set
                        int a = i, b = j;
                        if (a > b) std::swap(a, b);
                        existingEdges.insert({a, b});
                        break;
                    }
                }
            }
            // Update root
            root = find(0);
        }
    }

    // **Add Additional Edges**

    int numExistingEdges = edges.size();
    int numAdditionalEdges = static_cast<int>(numExistingEdges * additionalEdgePercentage);

// Prepare candidate edges
    std::vector<EdgeInfo> candidateEdges;

// Shuffle possibleEdges for randomness
    std::vector<EdgeInfo> shuffledPossibleEdges = possibleEdges;
    std::shuffle(shuffledPossibleEdges.begin(), shuffledPossibleEdges.end(), gen);

    for (const auto& edge : shuffledPossibleEdges) {
        if (numAdditionalEdges <= 0)
            break;

        int u = edge.vertexA;
        int v = edge.vertexB;
        int a = u, b = v;
        if (a > b) std::swap(a, b);

        if (existingEdges.count({a, b}) == 0) {
            // Directly add the edge without checking for intersection
            edges.push_back(edge);
            existingEdges.insert({a, b});
            numAdditionalEdges--;
        }
    }
}



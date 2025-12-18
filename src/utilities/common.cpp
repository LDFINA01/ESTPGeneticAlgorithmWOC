//
// Created by LuisF
//

#include "common.h"
#include "globalData.h"
#include <iostream>
#include <algorithm>
#include <execution>
#include "cmath"
#include <queue>


Timer::Timer() : running(false) {}

void Timer::start() {
    startTime = std::chrono::high_resolution_clock::now();
    running = true;
}

void Timer::stop() {
    if (running) {
        endTime = std::chrono::high_resolution_clock::now();
        running = false;
    }
}

double Timer::elapsedMilliseconds() const {
    if (running) {
        auto now = std::chrono::high_resolution_clock::now();
        return std::chrono::duration_cast<std::chrono::milliseconds>(now - startTime).count();
    } else {
        return std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();
    }
}

double Timer::elapsedSeconds() const {
    return elapsedMilliseconds() / 1000.0;
}

void printVector(const std::vector<int>& vec) {
    for (std::size_t i = 0; i < vec.size(); ++i) {
        std::cout << vec[i];
        if (i < vec.size() - 1) {
            std::cout << " -> ";
        }
    }
    std::cout << std::endl;
}

// Function to print a SteinerTree
void printSteinerTree(const SteinerTree& tree) {
    // Print the path
    std::cout << "Path: ";
    for (std::size_t i = 0; i < tree.path.size(); ++i) {
        std::cout << tree.path[i];
        if (i < tree.path.size() - 1) {
            std::cout << " -> ";
        }
    }

    // Print the total distance with fixed precision
    std::cout << "\nTotal Distance: " << tree.totalDistance << std::endl;
}

void printEdges() {
    for (const auto& edge : edges) {
        std::cout << "edge: (" << edge.vertexA << ", " << edge.vertexB << "), weight: " << edge.weight << ", frequency: " << edge.frequency << std::endl;
    }
}


double computeTotalDistance(const std::vector<int> &path) {
    double totalDistance = 0.0;
    for (size_t i = 1; i < path.size(); ++i) {
        int u = path[i - 1];
        int v = path[i];
        bool edgeFound = false;

        // Search for the edge between u and v
        for (const auto& edge : edges) {
            if ((edge.vertexA == u && edge.vertexB == v) || (edge.vertexA == v && edge.vertexB == u)) {
                totalDistance += edge.weight;
                edgeFound = true;
                break;
            }
        }

        // If edge not found in the edges list, compute Euclidean distance
        if (!edgeFound) {
            double dx = nodes[u].x - nodes[v].x;
            double dy = nodes[u].y - nodes[v].y;
            double dist = std::sqrt(dx * dx + dy * dy);
            totalDistance += dist;
        }
    }
    return totalDistance;
}

double computeTotalDistanceOfPopulation(const std::vector<SteinerTree>& population) {
    double totalPopulationDistance = 0.0;

    // Iterate over each SteinerTree in the population
    for (const SteinerTree& tree : population) {
        totalPopulationDistance += computeTotalDistance(tree.path); // Use the path from each SteinerTree
    }

    return totalPopulationDistance; // Return the total distance for the population
}

// Function to increment the frequency of edges used in a vector of Steiner Trees
void incrementEdgeFrequencies(const std::vector<SteinerTree>& trees) {
    // Create a map from edge key to index in the edges vector for quick lookup
    std::unordered_map<std::string, int> edgeKeyToIndex;
    for (size_t i = 0; i < edges.size(); ++i) {
        int u = edges[i].vertexA;
        int v = edges[i].vertexB;
        std::string key = (u < v) ? std::to_string(u) + "-" + std::to_string(v)
                                  : std::to_string(v) + "-" + std::to_string(u);
        edgeKeyToIndex[key] = static_cast<int>(i);
    }

    // Iterate over each tree in the vector
    for (const auto& tree : trees) {
        const std::vector<int>& path = tree.path; // The sequence of nodes in the tree

        // Increment frequency for each edge in the tree's path
        for (size_t i = 0; i + 1 < path.size(); ++i) {
            int u = path[i];
            int v = path[i + 1];
            std::string key = (u < v) ? std::to_string(u) + "-" + std::to_string(v)
                                      : std::to_string(v) + "-" + std::to_string(u);

            auto it = edgeKeyToIndex.find(key);
            if (it != edgeKeyToIndex.end()) {
                int idx = it->second;
                edges[idx].frequency += 1; // Increment the frequency of the edge
            }
        }
    }
}

// Helper function to check if an edge exists between two nodes
bool edgeExists(int u, int v) {
    for (const auto& edge : edges) {
        if ((edge.vertexA == u && edge.vertexB == v) || (edge.vertexA == v && edge.vertexB == u)) {
            return true;
        }
    }
    return false;
}

// Function to find a path between two nodes using BFS
std::vector<int> findPath(int start, int end) {
    std::queue<std::vector<int>> q;
    std::unordered_set<int> visited;
    q.push({start});
    visited.insert(start);

    while (!q.empty()) {
        std::vector<int> currentPath = q.front();
        q.pop();
        int currentNode = currentPath.back();

        if (currentNode == end) {
            return currentPath;
        }

        for (const auto& edge : edges) {
            int neighbor = -1;
            if (edge.vertexA == currentNode) {
                neighbor = edge.vertexB;
            } else if (edge.vertexB == currentNode) {
                neighbor = edge.vertexA;
            }
            if (neighbor != -1 && visited.find(neighbor) == visited.end()) {
                visited.insert(neighbor);
                std::vector<int> newPath = currentPath;
                newPath.push_back(neighbor);
                q.push(newPath);
            }
        }
    }
    // No path found
    return {};
}

// Function to find the closest node in the path to the given terminal
int findClosestNode(int terminal, const std::vector<int>& path) {
    double minDistance = std::numeric_limits<double>::max();
    int closestNode = -1;
    for (int node : path) {
        double dx = nodes[node].x - nodes[terminal].x;
        double dy = nodes[node].y - nodes[terminal].y;
        double dist = std::sqrt(dx * dx + dy * dy);
        if (dist < minDistance) {
            minDistance = dist;
            closestNode = node;
        }
    }
    return closestNode;
}

// Repair function to ensure path validity
std::vector<int> repairPath(const std::vector<int>& path) {
    std::unordered_set<int> visited;
    std::vector<int> repairedPath;
    int previousNode = -1;

    for (int node : path) {
        if (visited.find(node) == visited.end()) {
            // Ensure the edge between previousNode and node exists
            if (previousNode != -1 && !edgeExists(previousNode, node)) {
                // Find a path between previousNode and node
                std::vector<int> connectingPath = findPath(previousNode, node);
                if (!connectingPath.empty()) {
                    // Exclude the starting node to avoid duplication
                    repairedPath.insert(repairedPath.end(), connectingPath.begin() + 1, connectingPath.end());
                }
            } else {
                repairedPath.push_back(node);
            }
            visited.insert(node);
            previousNode = node;
        }
    }

    // Ensure all terminal vertices are included
    for (int terminal : terminalVertices) {
        if (visited.find(terminal) == visited.end()) {
            // Find the closest node in repairedPath
            int closestNode = findClosestNode(terminal, repairedPath);
            // Find a path to connect the terminal
            std::vector<int> connectingPath = findPath(closestNode, terminal);
            if (!connectingPath.empty()) {
                // Exclude the starting node to avoid duplication
                repairedPath.insert(repairedPath.end(), connectingPath.begin() + 1, connectingPath.end());
                visited.insert(terminal);
            }
        }
    }

    return repairedPath;
}

// Function to select top k Steiner Trees
std::vector<SteinerTree> selectTopK(std::vector<SteinerTree> &population) {
    size_t populationSize = population.size();
    size_t k = 0;

    // Calculate k based on population size
    if (populationSize > 500) {
        k = static_cast<size_t>(populationSize * 0.05); // 5% for population > 500
    } else if (populationSize > 250) {
        k = static_cast<size_t>(populationSize * 0.10); // 10% for 250 < population <= 500
    } else {
        k = static_cast<size_t>(populationSize * 0.15); // 15% for population <= 250
    }

    // Ensure k is at least 1 to avoid returning an empty vector
    k = std::max<size_t>(1, k);

    if (k >= populationSize) {
        // If k is larger than the population, return the whole population sorted
        std::sort(population.begin(), population.end(),
                  [](const SteinerTree& a, const SteinerTree& b) {
                      return a.totalDistance < b.totalDistance;
                  });
        return population;
    }

    // Partially sort the population to find the kth element
    std::nth_element(population.begin(), population.begin() + k, population.end(),
                     [](const SteinerTree& a, const SteinerTree& b) {
                         return a.totalDistance < b.totalDistance;
                     });

    // Extract the top k elements
    std::vector<SteinerTree> topK(population.begin(), population.begin() + k);

    // Optionally sort the top k elements if order matters
    std::sort(topK.begin(), topK.end(),
              [](const SteinerTree& a, const SteinerTree& b) {
                  return a.totalDistance < b.totalDistance;
              });

    return topK;
}




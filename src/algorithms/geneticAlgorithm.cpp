// geneticAlgorithm.cpp
// Created by LuisF on 11/14/2024.

#include "geneticAlgorithm.h"
#include "globalData.h"
#include "common.h"
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <random>
#include <algorithm>
#include <iostream>




#include <iostream>
#include <vector>
#include <unordered_set>
#include <random>
#include <algorithm>

// Global variables representing the graph
const int MAXN = 100;
std::vector<int> adj[MAXN]; // Adjacency list representation of the graph

// Randomly selects a starting vertex from the terminalVertices
int getRandomStart() {
    // Initialize random number generator
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(1, n);

    // Return a random number between 1 and n
    return dis(gen);
}


std::vector<int> removeConsecutiveDuplicates(const std::vector<int>& path) {
    if (path.empty()) return {};

    std::vector<int> finalPath;
    finalPath.push_back(path[0]);

    for (size_t i = 1; i < path.size(); ++i) {
        // Add the current city if it's not a consecutive duplicate
        if (path[i] != path[i - 1]) {
            finalPath.push_back(path[i]);
        }

        // Check if we can remove the last two cities due to backtracking
        if (i >= 2) {
            // Current city
            int current = path[i];
            // City two steps back
            int twoBack = path[i - 2];
            // City one step back
            int oneBack = path[i - 1];

            // Check if current city equals two steps back
            // and both current and oneBack are non-terminal
            if (current == twoBack &&
                terminalVertices.find(current) == terminalVertices.end() &&
                terminalVertices.find(oneBack) == terminalVertices.end()) {

                // Remove the last city (oneBack) from finalPath
                finalPath.pop_back();

                // Do not add the current city to avoid backtracking
                // Thus, finalPath remains up to twoBack
            }
        }
    }

    return finalPath;
}


// Function to perform DFS to find a valid Steiner tree
SteinerTree randomSteinerTree() {
    const int MAX_ITERATIONS = 100000; // Maximum allowed iterations before restarting

    std::random_device rd;
    std::mt19937 gen(rd());

    while (true) {
        int iterationCount = 0;

        std::unordered_set<int> visitedTerminals; // Track visited terminal vertices
        std::vector<int> visitedCities;           // Track all cities visited (with duplicates)
        std::vector<EdgeInfo> availableEdges;     // Track available edges from the current city
        std::vector<int> finalPath;               // Write out the official path of the Steiner Tree

        // Choose a random starting vertex from terminal vertices
        int currentCity = getRandomStart();

        while (true) {
            iterationCount++;

            // Add the current city to the visited list
            visitedCities.push_back(currentCity);

            // If the current city is a terminal vertex, add it to the visited terminals set
            if (terminalVertices.count(currentCity)) {
                visitedTerminals.insert(currentCity);
            }

            // Gather all available edges connected to the current city
            availableEdges.clear();
            for (const auto& edge : edges) {
                if (edge.vertexA == currentCity || edge.vertexB == currentCity) {
                    availableEdges.push_back(edge);
                }
            }

            // If no more available edges, break the loop
            if (availableEdges.empty()) {
                break;
            }

            // Randomly pick an edge from the available edges
            std::uniform_int_distribution<> dis(0, availableEdges.size() - 1);
            EdgeInfo chosenEdge = availableEdges[dis(gen)];

            // Determine the next city to visit based on the chosen edge
            int nextCity = (chosenEdge.vertexA == currentCity) ? chosenEdge.vertexB : chosenEdge.vertexA;

            currentCity = nextCity; // Move to the next city

            // Stop if all terminal vertices are visited
            if (visitedTerminals.size() == terminalVertices.size()) {
                finalPath = removeConsecutiveDuplicates(visitedCities);

                double totalLength = computeTotalDistance(finalPath);
                SteinerTree tree = {finalPath, totalLength};
                return tree;

            }

            // If maximum iterations reached, restart the function
            if (iterationCount >= MAX_ITERATIONS) {
                std::cerr << "Maximum iterations reached. Restarting the traversal." << std::endl;
                break; // Break out of the inner while loop to restart
            }
        }

        // If we haven't returned by now, we need to restart
        // Continue the outer while loop to restart the traversal
    }
}

// Function to perform crossover between two Steiner trees
SteinerTree crossoverFunction(const SteinerTree &parentA, const SteinerTree &parentB) {
    // Initialize random number generator
    std::random_device rd;
    std::mt19937 gen(rd());

    // Ensure parents have paths
    if (parentA.path.empty() || parentB.path.empty()) {
        // Return the non-empty parent or an empty tree
        return parentA.path.empty() ? parentB : parentA;
    }

    // Step 1: Find valid crossover points where the nodes are connected
    std::vector<size_t> validCrossoverPointsA;
    std::vector<size_t> validCrossoverPointsB;

    // Build sets for faster lookup
    std::unordered_set<int> parentA_nodes(parentA.path.begin(), parentA.path.end());
    std::unordered_set<int> parentB_nodes(parentB.path.begin(), parentB.path.end());

    // Find all nodes common to both parents
    std::vector<int> commonNodes;
    for (int node : parentA_nodes) {
        if (parentB_nodes.count(node)) {
            commonNodes.push_back(node);
        }
    }

    // If there are common nodes, use them as crossover points
    if (!commonNodes.empty()) {
        // Randomly select a common node as the crossover point
        std::uniform_int_distribution<> disCommon(0, commonNodes.size() - 1);
        int crossoverNode = commonNodes[disCommon(gen)];

        // Find indices of the crossover node in both parents
        auto itA = std::find(parentA.path.begin(), parentA.path.end(), crossoverNode);
        auto itB = std::find(parentB.path.begin(), parentB.path.end(), crossoverNode);

        size_t indexA = std::distance(parentA.path.begin(), itA);
        size_t indexB = std::distance(parentB.path.begin(), itB);

        // Create offspring by combining the subpaths at the crossover node
        std::vector<int> offspringPath;

        // Add subpath from parent A up to and including the crossover node
        offspringPath.insert(offspringPath.end(), parentA.path.begin(), parentA.path.begin() + indexA + 1);

        // Add subpath from parent B after the crossover node
        offspringPath.insert(offspringPath.end(), parentB.path.begin() + indexB + 1, parentB.path.end());

        // Repair the offspring path to ensure all terminals are included and edges exist
        offspringPath = repairPath(offspringPath);

        // Compute total distance
        double totalDistance = computeTotalDistance(offspringPath);

        return {offspringPath, totalDistance};
    } else {
        // No common nodes; find valid crossover points where nodes are connected by existing edges
        std::vector<std::pair<size_t, size_t>> connectedCrossoverPoints;

        for (size_t i = 1; i < parentA.path.size(); ++i) {
            int nodeA = parentA.path[i];
            for (size_t j = 1; j < parentB.path.size(); ++j) {
                int nodeB = parentB.path[j];
                if (edgeExists(nodeA, nodeB)) {
                    connectedCrossoverPoints.emplace_back(i, j);
                }
            }
        }

        if (!connectedCrossoverPoints.empty()) {
            // Randomly select a connected crossover point
            std::uniform_int_distribution<> disConnected(0, connectedCrossoverPoints.size() - 1);
            auto [indexA, indexB] = connectedCrossoverPoints[disConnected(gen)];

            // Create offspring by combining subpaths and connecting crossover points
            std::vector<int> offspringPath;

            // Add subpath from parent A up to indexA
            offspringPath.insert(offspringPath.end(), parentA.path.begin(), parentA.path.begin() + indexA + 1);

            // Add node from parent B at indexB to connect
            offspringPath.push_back(parentB.path[indexB]);

            // Add subpath from parent B after indexB
            offspringPath.insert(offspringPath.end(), parentB.path.begin() + indexB + 1, parentB.path.end());

            // Repair the offspring path
            offspringPath = repairPath(offspringPath);

            // Compute total distance
            double totalDistance = computeTotalDistance(offspringPath);

            return {offspringPath, totalDistance};
        } else {
            // If no valid crossover points, return one of the parents
            return (parentA.totalDistance < parentB.totalDistance) ? parentA : parentB;
        }
    }
}

// Function to find the shortest path between two nodes using BFS
std::vector<int> findShortestPath(int start, int end, const std::unordered_map<int, std::vector<int>>& adjacencyMap) {
    std::queue<std::vector<int>> q;
    std::unordered_set<int> visited;
    q.push({start});
    visited.insert(start);

    while (!q.empty()) {
        std::vector<int> path = q.front();
        q.pop();

        int currentNode = path.back();
        if (currentNode == end) {
            return path;
        }

        for (int neighbor : adjacencyMap.at(currentNode)) {
            if (visited.find(neighbor) == visited.end()) {
                visited.insert(neighbor);
                std::vector<int> newPath = path;
                newPath.push_back(neighbor);
                q.push(newPath);
            }
        }
    }
    // No path found
    return {};
}

// Function to perform random mutations on a vector of Steiner Trees
void randomMutation(std::vector<SteinerTree>& population, double mutationRate) {
    // Build adjacency map
    std::unordered_map<int, std::vector<int>> adjacencyMap;
    for (const auto& edge : edges) {
        adjacencyMap[edge.vertexA].push_back(edge.vertexB);
        adjacencyMap[edge.vertexB].push_back(edge.vertexA);
    }

    // Initialize random number generator
    std::random_device rd;
    std::mt19937 gen(rd());

    for (auto& tree : population) {
        std::uniform_real_distribution<> dis(0.5, 1.0);
        if (dis(gen) < mutationRate) {
            // Apply mutation to this tree
            // Decide which mutation to apply
            std::uniform_int_distribution<> mutationTypeDis(1, 3);
            int mutationType = mutationTypeDis(gen);

            if (mutationType == 1) {
                // **Mutation Type 1:** Randomly remove a non-terminal node from the path

                std::vector<size_t> nonTerminalIndices;
                for (size_t i = 1; i < tree.path.size() - 1; ++i) {
                    int node = tree.path[i];
                    if (terminalVertices.find(node) == terminalVertices.end()) {
                        nonTerminalIndices.push_back(i);
                    }
                }

                if (!nonTerminalIndices.empty()) {
                    std::uniform_int_distribution<> indexDis(0, nonTerminalIndices.size() - 1);
                    size_t idxToRemove = nonTerminalIndices[indexDis(gen)];

                    // Remove the node
                    tree.path.erase(tree.path.begin() + idxToRemove);

                    // Check if the two adjacent nodes are connected
                    if (idxToRemove > 0 && idxToRemove < tree.path.size()) {
                        int prevNode = tree.path[idxToRemove - 1];
                        int nextNode = tree.path[idxToRemove];

                        if (edgeExists(prevNode, nextNode)) {
                            // Nodes are connected, path remains valid
                        } else {
                            // Need to find a path between prevNode and nextNode
                            std::vector<int> connectingPath = findShortestPath(prevNode, nextNode, adjacencyMap);
                            if (!connectingPath.empty()) {
                                // Replace the edge between prevNode and nextNode with connectingPath
                                tree.path.erase(tree.path.begin() + idxToRemove, tree.path.begin() + idxToRemove);
                                tree.path.insert(tree.path.begin() + idxToRemove, connectingPath.begin() + 1, connectingPath.end() - 1);
                            } else {
                                // Cannot connect prevNode and nextNode, revert mutation
                                // (Alternatively, could try another mutation)
                                continue;
                            }
                        }
                    }
                }

            } else if (mutationType == 2) {
                // **Mutation Type 2:** Randomly add a node between two nodes in the path, if possible

                if (tree.path.size() >= 2) {
                    std::uniform_int_distribution<> indexDis(0, tree.path.size() - 2);
                    size_t idx = indexDis(gen);

                    int u = tree.path[idx];
                    int v = tree.path[idx + 1];

                    // Find nodes connected to both u and v
                    std::vector<int> candidateNodes;
                    for (int node : adjacencyMap[u]) {
                        if (node != v && std::find(adjacencyMap[v].begin(), adjacencyMap[v].end(), node) != adjacencyMap[v].end()) {
                            candidateNodes.push_back(node);
                        }
                    }

                    if (!candidateNodes.empty()) {
                        std::uniform_int_distribution<> nodeDis(0, candidateNodes.size() - 1);
                        int newNode = candidateNodes[nodeDis(gen)];

                        // Insert newNode between u and v
                        tree.path.insert(tree.path.begin() + idx + 1, newNode);
                    }
                }

            } else if (mutationType == 3) {
                // **Mutation Type 3:** Rewire a portion of the tree

                if (tree.path.size() >= 4) {
                    std::uniform_int_distribution<> indexDis(0, tree.path.size() - 1);
                    size_t idx1 = indexDis(gen);
                    size_t idx2 = indexDis(gen);

                    if (idx1 > idx2) std::swap(idx1, idx2);
                    if (idx2 - idx1 < 2) continue; // Need at least one node in between

                    int u = tree.path[idx1];
                    int v = tree.path[idx2];

                    // Find alternative path between u and v
                    std::vector<int> altPath = findShortestPath(u, v, adjacencyMap);

                    if (!altPath.empty()) {
                        // Replace the subpath between idx1 and idx2 with altPath
                        tree.path.erase(tree.path.begin() + idx1 + 1, tree.path.begin() + idx2);
                        tree.path.insert(tree.path.begin() + idx1 + 1, altPath.begin() + 1, altPath.end() - 1);
                    }
                }
            }

            // After mutation, ensure all terminals are included
            std::unordered_set<int> pathSet(tree.path.begin(), tree.path.end());
            for (int terminal : terminalVertices) {
                if (pathSet.find(terminal) == pathSet.end()) {
                    // Need to connect terminal to the tree
                    int closestNode = -1;
                    double minDistance = std::numeric_limits<double>::max();
                    for (int node : pathSet) {
                        double dx = nodes[node].x - nodes[terminal].x;
                        double dy = nodes[node].y - nodes[terminal].y;
                        double dist = std::sqrt(dx * dx + dy * dy);
                        if (dist < minDistance) {
                            minDistance = dist;
                            closestNode = node;
                        }
                    }

                    if (closestNode != -1) {
                        // Find a path from closestNode to terminal
                        std::vector<int> connectingPath = findShortestPath(closestNode, terminal, adjacencyMap);
                        if (!connectingPath.empty()) {
                            // Insert connectingPath into tree.path
                            tree.path.insert(tree.path.end(), connectingPath.begin() + 1, connectingPath.end());
                            pathSet.insert(connectingPath.begin(), connectingPath.end());
                        }
                    }
                }
            }

            // Remove duplicates and recompute total distance
            tree.path = removeConsecutiveDuplicates(tree.path);
            tree.totalDistance = computeTotalDistance(tree.path);
        }
    }
}

std::vector<SteinerTree> initialPopulation (int size) {
    std::vector<SteinerTree> initialPopulation;

    // Progress milestones
    int milestone25 = size * 0.25;
    int milestone50 = size * 0.50;
    int milestone75 = size * 0.75;
    int milestone90 = size * 0.90;

    for (int i = 0; i < size; i++) {
        initialPopulation.push_back(randomSteinerTree());

//        // Print progress at milestones
//        if (i == milestone25) {
//            std::cout << "25% of the initial population generated." << std::endl;
//        } else if (i == milestone50) {
//            std::cout << "50% of the initial population generated." << std::endl;
//        } else if (i == milestone75) {
//            std::cout << "75% of the initial population generated." << std::endl;
//        }
//    }
//    std::cout << "initial population generated.\n" << std::endl;
    }
    return initialPopulation;
}

std::vector<SteinerTree> runGeneticAlgorithm(int populationSize, int generations) {
    std::vector<SteinerTree> parentPopulation = initialPopulation(populationSize);
    int p = parentPopulation.size();


    double minMutationRate = 0.01; // 1%
    double maxMutationRate = 0.5;  // 50%

    // Progress milestones
    int milestone25 = generations * 0.25;
    int milestone50 = generations * 0.50;
    int milestone75 = generations * 0.75;
    int milestone90 = generations * 0.90;

    int currentGeneration = 0;
    while(currentGeneration < generations) {
        avgDistanceOverGenerations.push_back(computeTotalDistanceOfPopulation(parentPopulation));

        // Adjust mutation rate based on generation
        double mutationRate;
        if (currentGeneration < generations / 2) {
            mutationRate = minMutationRate;
        } else {
            double progress = (currentGeneration - generations / 2) / static_cast<double>(generations / 2);
            mutationRate = minMutationRate + progress * (maxMutationRate - minMutationRate);
        }

        std::vector<SteinerTree> fitPopulation = selectTopK(parentPopulation);
        incrementEdgeFrequencies(fitPopulation);
        std::vector<SteinerTree> childPopulation;

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> dis(0, fitPopulation.size() - 1);

        while (childPopulation.size() < p) {
            int parentAIdx = dis(gen);
            int parentBIdx = dis(gen);

            // Ensure parents are different
            while (parentAIdx == parentBIdx) {
                parentBIdx = dis(gen);
            }

            // Create offspring using crossoverFunction
            SteinerTree offspring = crossoverFunction(fitPopulation[parentAIdx], fitPopulation[parentBIdx]);
            childPopulation.push_back(offspring);
        }

        // Mutation
        randomMutation(childPopulation, mutationRate); // Using a 10% mutation rate

        // Prepare for next generation
        parentPopulation = childPopulation;

//        // Print progress at milestones
//        if (currentGeneration == milestone25) {
//            std::cout << "25% of generations completed." << std::endl;
//        } else if (currentGeneration == milestone50) {
//            std::cout << "50% of generations completed." << std::endl;
//        } else if (currentGeneration == milestone75) {
//            std::cout << "75% of generations completed." << std::endl;
//        }

        currentGeneration++;
    }

//    std::cout << "generations completed.\n" << std::endl;
    return parentPopulation;
}



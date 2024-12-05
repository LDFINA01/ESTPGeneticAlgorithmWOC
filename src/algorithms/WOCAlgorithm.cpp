#include "WOCAlgorithm.h"
#include "common.h"
#include <algorithm>
#include <queue>
#include <stack>
#include <unordered_map>
#include <random>
#include <cmath>
#include <iostream>

// Public method to generate Steiner Tree based on selected strategy
SteinerTree WOCAlgorithm::generateSteinerTree(Strategy strategy) {
    switch (strategy) {
        case Strategy::MSTBased:
            return generateMSTBasedSteinerTree();
        case Strategy::SPH:
            return generateSPHSteinerTree();
        case Strategy::Probabilistic:
            return generateProbabilisticSteinerTree();
        case Strategy::IterativeAddition:
            return generateIterativeAdditionSteinerTree();
        case Strategy::Threshold:
            return generateThresholdSteinerTree();
        default:
            std::cerr << "Invalid strategy selected. Returning empty Steiner Tree." << std::endl;
            return {};
    }
}

/////////////////////////////////////////////////////////////////
// Method 1: Minimum Spanning Tree (MST) Based Algorithm
/////////////////////////////////////////////////////////////////

SteinerTree WOCAlgorithm::generateMSTBasedSteinerTree() {
    // Adjust edge weights based on frequencies
    std::vector<EdgeInfo> adjustedEdges = edges;
    for (auto& edge : adjustedEdges) {
        double freq = edge.frequency > 0 ? static_cast<double>(edge.frequency) : 1.0;
        edge.weight = edge.weight / freq;
    }

    // Sort edges by adjusted weight (ascending)
    std::sort(adjustedEdges.begin(), adjustedEdges.end(), [](const EdgeInfo& a, const EdgeInfo& b) {
        return a.weight < b.weight;
    });

    // Kruskal's algorithm to compute MST
    std::vector<int> parent(n);
    std::iota(parent.begin(), parent.end(), 0);

    auto find_set = [&](int u) -> int {
        while (parent[u] != u)
            u = parent[u];
        return u;
    };

    auto unite_set = [&](int u, int v) {
        int pu = find_set(u);
        int pv = find_set(v);
        if (pu != pv)
            parent[pu] = pv;
    };

    std::vector<EdgeInfo> mstEdges;
    for (const auto& edge : adjustedEdges) {
        int u = edge.vertexA;
        int v = edge.vertexB;
        if (find_set(u) != find_set(v)) {
            unite_set(u, v);
            mstEdges.push_back(edge);
        }
    }

    // Build adjacency list from MST edges
    std::unordered_map<int, std::vector<int>> graph;
    for (const auto& edge : mstEdges) {
        graph[edge.vertexA].push_back(edge.vertexB);
        graph[edge.vertexB].push_back(edge.vertexA);
    }

    // Prune non-terminal leaves
    bool removed = true;
    while (removed) {
        removed = false;
        std::vector<int> leaves;
        for (const auto& node : graph) {
            if (graph[node.first].size() == 1 && terminalVertices.find(node.first) == terminalVertices.end()) {
                leaves.push_back(node.first);
            }
        }
        for (int leaf : leaves) {
            if (graph.find(leaf) != graph.end()) {
                int neighbor = graph[leaf][0];
                graph[neighbor].erase(std::remove(graph[neighbor].begin(), graph[neighbor].end(), leaf), graph[neighbor].end());
                graph.erase(leaf);
                removed = true;
            }
        }
    }

    // Extract path using DFS
    std::vector<int> path;
    std::unordered_set<int> visited;
    std::stack<int> stackDFS;
    int start = *terminalVertices.begin();
    stackDFS.push(start);
    visited.insert(start);
    while (!stackDFS.empty()) {
        int current = stackDFS.top();
        stackDFS.pop();
        path.push_back(current);
        for (int neighbor : graph[current]) {
            if (visited.find(neighbor) == visited.end()) {
                visited.insert(neighbor);
                stackDFS.push(neighbor);
            }
        }
    }

    // Compute total distance
    double totalDistance = computeTotalDistance(path);

    return {path, totalDistance};
}

/////////////////////////////////////////////////////////////////
// Method 2: Shortest Path Heuristic (SPH) Algorithm
/////////////////////////////////////////////////////////////////

SteinerTree WOCAlgorithm::generateSPHSteinerTree() {
    // Adjust edge weights based on frequencies
    std::vector<EdgeInfo> adjustedEdges = edges;
    for (auto& edge : adjustedEdges) {
        double freq = edge.frequency > 0 ? static_cast<double>(edge.frequency) : 1.0;
        edge.weight = edge.weight / freq;
    }

    // Build adjacency list with adjusted weights
    std::unordered_map<int, std::vector<std::pair<int, double>>> graph;
    for (const auto& edge : adjustedEdges) {
        graph[edge.vertexA].emplace_back(edge.vertexB, edge.weight);
        graph[edge.vertexB].emplace_back(edge.vertexA, edge.weight);
    }

    // Compute shortest paths between all pairs of terminals using Dijkstra's algorithm
    std::unordered_map<int, std::unordered_map<int, std::vector<int>>> shortestPaths;
    for (int terminal : terminalVertices) {
        // Dijkstra's algorithm from terminal
        std::unordered_map<int, double> dist;
        std::unordered_map<int, int> prev;
        for (int i = 0; i < n; ++i)
            dist[i] = std::numeric_limits<double>::infinity();
        dist[terminal] = 0.0;

        // Priority queue: (distance, node)
        using pii = std::pair<double, int>;
        std::priority_queue<pii, std::vector<pii>, std::greater<pii>> pq;
        pq.emplace(0.0, terminal);

        while (!pq.empty()) {
            auto [currentDist, u] = pq.top();
            pq.pop();

            if (currentDist > dist[u])
                continue;

            for (const auto& [v, weight] : graph[u]) {
                if (dist[u] + weight < dist[v]) {
                    dist[v] = dist[u] + weight;
                    prev[v] = u;
                    pq.emplace(dist[v], v);
                }
            }
        }

        // Store paths to other terminals
        for (int t : terminalVertices) {
            if (t == terminal) continue;
            if (dist[t] == std::numeric_limits<double>::infinity()) continue; // No path
            std::vector<int> path;
            int current = t;
            while (current != terminal) {
                path.push_back(current);
                current = prev[current];
            }
            path.push_back(terminal);
            std::reverse(path.begin(), path.end());
            shortestPaths[terminal][t] = path;
        }
    }

    // Collect all terminal edges with their paths
    struct TerminalEdge {
        int u, v;
        double weight;
        std::vector<int> path;
    };
    std::vector<TerminalEdge> terminalEdges;
    for (int u : terminalVertices) {
        for (int v : terminalVertices) {
            if (u >= v) continue;
            if (shortestPaths[u].find(v) != shortestPaths[u].end()) {
                double weight = computeTotalDistance(shortestPaths[u][v]);
                terminalEdges.push_back({u, v, weight, shortestPaths[u][v]});
            }
        }
    }

    // Kruskal's algorithm on terminal edges
    std::sort(terminalEdges.begin(), terminalEdges.end(), [](const TerminalEdge& a, const TerminalEdge& b) {
        return a.weight < b.weight;
    });

    std::unordered_map<int, int> parent;
    for (int t : terminalVertices)
        parent[t] = t;

    auto find_set_sph = [&](int u) -> int {
        while (parent[u] != u)
            u = parent[u];
        return u;
    };

    auto unite_set_sph = [&](int u, int v) {
        parent[find_set_sph(u)] = find_set_sph(v);
    };

    std::vector<int> steinerPath;
    for (const auto& edge : terminalEdges) {
        if (find_set_sph(edge.u) != find_set_sph(edge.v)) {
            unite_set_sph(edge.u, edge.v);
            steinerPath.insert(steinerPath.end(), edge.path.begin(), edge.path.end());
        }
    }

    // Remove duplicates while preserving order
    std::vector<int> uniquePath;
    std::unordered_set<int> seen;
    for (int node : steinerPath) {
        if (seen.find(node) == seen.end()) {
            uniquePath.push_back(node);
            seen.insert(node);
        }
    }

    // Compute total distance
    double totalDistance = computeTotalDistance(uniquePath);

    return {uniquePath, totalDistance};
}

/////////////////////////////////////////////////////////////////
// Method 3: Probabilistic Edge Selection Algorithm
/////////////////////////////////////////////////////////////////

SteinerTree WOCAlgorithm::generateProbabilisticSteinerTree() {
    // Compute total frequency
    double totalFrequency = 0.0;
    for (const auto& edge : edges)
        totalFrequency += static_cast<double>(edge.frequency);

    if (totalFrequency == 0.0) {
        std::cerr << "All edge frequencies are zero. Cannot proceed with Probabilistic Algorithm." << std::endl;
        return {};
    }

    // Create cumulative distribution
    std::vector<double> cumulativeFreq;
    cumulativeFreq.push_back(0.0);
    for (const auto& edge : edges)
        cumulativeFreq.push_back(cumulativeFreq.back() + static_cast<double>(edge.frequency) / totalFrequency);

    // Initialize Union-Find structure
    std::vector<int> parentUF(n);
    std::iota(parentUF.begin(), parentUF.end(), 0);

    auto find_set_p = [&](int u) -> int {
        while (parentUF[u] != u)
            u = parentUF[u];
        return u;
    };

    auto unite_set_p = [&](int u, int v) {
        parentUF[find_set_p(u)] = find_set_p(v);
    };

    // Initialize graph
    std::unordered_map<int, std::vector<int>> graph;

    // Random number generator
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);

    // Select edges probabilistically
    while (true) {
        double r = dis(gen);
        auto it = std::lower_bound(cumulativeFreq.begin(), cumulativeFreq.end(), r);
        int index = std::distance(cumulativeFreq.begin(), it) - 1;
        if (index < 0 || index >= static_cast<int>(edges.size())) continue;

        const auto& edge = edges[index];
        int u = edge.vertexA;
        int v = edge.vertexB;

        if (find_set_p(u) != find_set_p(v)) {
            unite_set_p(u, v);
            graph[u].push_back(v);
            graph[v].push_back(u);
        }

        // Check if all terminals are connected
        int parentTerminal = find_set_p(*terminalVertices.begin());
        bool allConnected = std::all_of(terminalVertices.begin(), terminalVertices.end(), [&](int t) {
            return find_set_p(t) == parentTerminal;
        });
        if (allConnected) break;
    }

    // Extract path using DFS
    std::vector<int> path;
    std::unordered_set<int> visited;
    std::stack<int> stackDFS;
    int start = *terminalVertices.begin();
    stackDFS.push(start);
    visited.insert(start);
    while (!stackDFS.empty()) {
        int current = stackDFS.top();
        stackDFS.pop();
        path.push_back(current);
        for (int neighbor : graph[current]) {
            if (visited.find(neighbor) == visited.end()) {
                visited.insert(neighbor);
                stackDFS.push(neighbor);
            }
        }
    }

    // Compute total distance
    double totalDistance = computeTotalDistance(path);

    return {path, totalDistance};
}

/////////////////////////////////////////////////////////////////
// Method 4: Iterative Edge Addition Algorithm
/////////////////////////////////////////////////////////////////

SteinerTree WOCAlgorithm::generateIterativeAdditionSteinerTree() {
    // Sort edges by frequency descending
    std::vector<EdgeInfo> sortedEdges = edges;
    std::sort(sortedEdges.begin(), sortedEdges.end(), [](const EdgeInfo& a, const EdgeInfo& b) {
        return a.frequency > b.frequency;
    });

    // Initialize Union-Find structure
    std::vector<int> parentUF(n);
    std::iota(parentUF.begin(), parentUF.end(), 0);

    auto find_set_i = [&](int u) -> int {
        while (parentUF[u] != u)
            u = parentUF[u];
        return u;
    };

    auto unite_set_i = [&](int u, int v) {
        parentUF[find_set_i(u)] = find_set_i(v);
    };

    // Initialize graph
    std::unordered_map<int, std::vector<int>> graph;

    // Add edges iteratively
    for (const auto& edge : sortedEdges) {
        int u = edge.vertexA;
        int v = edge.vertexB;
        if (find_set_i(u) != find_set_i(v)) {
            unite_set_i(u, v);
            graph[u].push_back(v);
            graph[v].push_back(u);
        }

        // Check if all terminals are connected
        int parentTerminal = find_set_i(*terminalVertices.begin());
        bool allConnected = std::all_of(terminalVertices.begin(), terminalVertices.end(), [&](int t) {
            return find_set_i(t) == parentTerminal;
        });
        if (allConnected) break;
    }

    // Extract path using DFS
    std::vector<int> path;
    std::unordered_set<int> visited;
    std::stack<int> stackDFS;
    int start = *terminalVertices.begin();
    stackDFS.push(start);
    visited.insert(start);
    while (!stackDFS.empty()) {
        int current = stackDFS.top();
        stackDFS.pop();
        path.push_back(current);
        for (int neighbor : graph[current]) {
            if (visited.find(neighbor) == visited.end()) {
                visited.insert(neighbor);
                stackDFS.push(neighbor);
            }
        }
    }

    // Compute total distance
    double totalDistance = computeTotalDistance(path);

    return {path, totalDistance};
}

/////////////////////////////////////////////////
// Method 5: Edge Frequency Thresholding Algorithm
/////////////////////////////////////////////////

SteinerTree WOCAlgorithm::generateThresholdSteinerTree() {
    // Calculate mean frequency
    double totalFrequency = 0.0;
    for (const auto& edge : edges)
        totalFrequency += static_cast<double>(edge.frequency);
    double meanFrequency = edges.empty() ? 0.0 : totalFrequency / static_cast<double>(edges.size());

    // Filter edges with frequency >= meanFrequency
    std::vector<EdgeInfo> highFreqEdges;
    for (const auto& edge : edges) {
        if (static_cast<double>(edge.frequency) >= meanFrequency)
            highFreqEdges.push_back(edge);
    }

    if (highFreqEdges.empty()) {
        std::cerr << "No edges above frequency threshold. Returning empty Steiner Tree." << std::endl;
        return {};
    }

    // Sort high-frequency edges by weight ascending
    std::sort(highFreqEdges.begin(), highFreqEdges.end(), [](const EdgeInfo& a, const EdgeInfo& b) {
        return a.weight < b.weight;
    });

    // Initialize Union-Find structure
    std::vector<int> parentUF(n);
    std::iota(parentUF.begin(), parentUF.end(), 0);

    auto find_set = [&](int u) -> int {
        while (parentUF[u] != u)
            u = parentUF[u];
        return u;
    };

    auto unite_set = [&](int u, int v) {
        parentUF[find_set(u)] = find_set(v);
    };

    // Initialize graph from selected edges
    std::unordered_map<int, std::vector<int>> selectedGraph;

    // Add high-frequency edges
    for (const auto& edge : highFreqEdges) {
        int u = edge.vertexA;
        int v = edge.vertexB;
        if (find_set(u) != find_set(v)) {
            unite_set(u, v);
            selectedGraph[u].push_back(v);
            selectedGraph[v].push_back(u);
        }

        // Check if all terminals are connected
        int parentTerminal = find_set(*terminalVertices.begin());
        bool allConnected = std::all_of(terminalVertices.begin(), terminalVertices.end(), [&](int t) {
            return find_set(t) == parentTerminal;
        });
        if (allConnected) break;
    }

    // If not all terminals are connected, add lower-frequency edges
    int parentTerminalCheck = find_set(*terminalVertices.begin());
    bool connected = std::all_of(terminalVertices.begin(), terminalVertices.end(), [&](int t) {
        return find_set(t) == parentTerminalCheck;
    });

    if (!connected) {
        // Collect lower-frequency edges
        std::vector<EdgeInfo> lowerFreqEdges;
        for (const auto& edge : edges) {
            if (static_cast<double>(edge.frequency) < meanFrequency)
                lowerFreqEdges.push_back(edge);
        }

        // Sort lower-frequency edges by weight ascending
        std::sort(lowerFreqEdges.begin(), lowerFreqEdges.end(), [](const EdgeInfo& a, const EdgeInfo& b) {
            return a.weight < b.weight;
        });

        // Add lower-frequency edges
        for (const auto& edge : lowerFreqEdges) {
            int u = edge.vertexA;
            int v = edge.vertexB;
            if (find_set(u) != find_set(v)) {
                unite_set(u, v);
                selectedGraph[u].push_back(v);
                selectedGraph[v].push_back(u);
            }

            // Check if all terminals are connected
            int parentTerm = find_set(*terminalVertices.begin());
            bool allConn = std::all_of(terminalVertices.begin(), terminalVertices.end(), [&](int t) {
                return find_set(t) == parentTerm;
            });
            if (allConn) break;
        }
    }

    // Final connectivity check
    parentTerminalCheck = find_set(*terminalVertices.begin());
    connected = std::all_of(terminalVertices.begin(), terminalVertices.end(), [&](int t) {
        return find_set(t) == parentTerminalCheck;
    });

    if (!connected) {
        std::cerr << "Unable to connect all terminal vertices. Returning empty Steiner Tree." << std::endl;
        return {};
    }

    // Extract a valid path using BFS to ensure all edges are valid
    int start = *terminalVertices.begin();
    std::vector<int> path = findValidPath(selectedGraph, start);

    if (path.empty()) {
        std::cerr << "Failed to extract a valid path from the Steiner Tree." << std::endl;
        return {};
    }

    // Compute total distance
    double totalDistance = 0.0;
    for (size_t i = 0; i < path.size() - 1; ++i) {
        int u = path[i];
        int v = path[i + 1];
        // Find the edge in the edge list to get its weight
        auto it = std::find_if(edges.begin(), edges.end(), [&](const EdgeInfo& e) {
            return (e.vertexA == u && e.vertexB == v) ||
                   (e.vertexA == v && e.vertexB == u);
        });
        if (it != edges.end()) {
            totalDistance += it->weight;
        } else {
            // This should not happen as we built the graph from existing edges
            std::cerr << "Edge (" << u << ", " << v << ") not found in edge list." << std::endl;
        }
    }

    return {path, totalDistance};
}

//////////////////////////////////////////////////////////
// Helper Function to Find a Valid Path Using BFS
//////////////////////////////////////////////////////////

std::vector<int> WOCAlgorithm::findValidPath(const std::unordered_map<int, std::vector<int>>& graph, int start) {
    std::vector<int> path;
    std::unordered_set<int> visited;
    std::queue<int> q;
    std::unordered_map<int, int> parentMap; // To reconstruct the path

    q.push(start);
    visited.insert(start);

    while (!q.empty()) {
        int current = q.front();
        q.pop();
        path.push_back(current);

        for (int neighbor : graph.at(current)) {
            if (visited.find(neighbor) == visited.end()) {
                visited.insert(neighbor);
                parentMap[neighbor] = current;
                q.push(neighbor);
            }
        }
    }

    // Reconstruct the path to include all terminals
    std::vector<int> finalPath;
    std::unordered_set<int> includedTerminals = terminalVertices;

    for (int terminal : includedTerminals) {
        if (terminal == start) continue;
        if (parentMap.find(terminal) == parentMap.end()) {
            std::cerr << "Terminal " << terminal << " is not reachable from " << start << "." << std::endl;
            continue;
        }

        // Reconstruct path from terminal to start
        std::vector<int> tempPath;
        int current = terminal;
        while (current != start) {
            tempPath.push_back(current);
            current = parentMap[current];
        }
        tempPath.push_back(start);
        std::reverse(tempPath.begin(), tempPath.end());

        // Append to finalPath without duplicating vertices
        for (int vertex : tempPath) {
            if (finalPath.empty() || finalPath.back() != vertex) {
                finalPath.push_back(vertex);
            }
        }
    }

    return finalPath;
}
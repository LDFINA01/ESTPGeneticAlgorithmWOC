// main.cpp

#include "gui/graphVisualizer.h"
#include "algorithms/geneticAlgorithm.h"
#include "algorithms/WOCAlgorithm.h"
#include "utilities/randomWeightedGraph.h"
#include "utilities/globalData.h"
#include "common.h"
#include <iostream>

int main() {
    Timer timer;

    // Parameters for the genetic algorithm
    int populationSize;
    int generations = 30;

    // Vector to store average distances for each node count
    std::vector<double> averageDistances;
    
    // Store best trees for visualization (only show GUI once)
    SteinerTree bestSteinerTree;
    SteinerTree bestWOCTree;
    bool showGUI = false;

    // Loop over node counts: 20, 40, 60, 80, 100
    for (int v = 20; v <= 40; v += 20) {
        double totalDistanceForVertex = 0.0;
        int totalRuns = 0;

        // Set the population size equal to the number of nodes
        populationSize = v;

        // Run 5 different random graphs for each node count
        for (int g = 0; g < 1; ++g) {
            // Set the global variable 'n' to the current number of nodes
            n = v;

            // Generate a random weighted graph
            randomWeightedGraph graph;
            graph.generateGraph();

            // Optional: Output the number of terminals if needed
            std::cout << "\nNumber of nodes: " << n << "\n";
            std::cout << "Number of terminal nodes: " << numTerminals << "\n";

            // Run the genetic algorithm 5 times on the same graph
            for (int ga = 0; ga < 1; ++ga) {
                timer.start();

                // Run the genetic algorithm
                std::vector<SteinerTree> best = runGeneticAlgorithm(populationSize, generations);

                timer.stop();

                // Assume best[0] is the best solution
                double distance = best[0].totalDistance;
                totalDistanceForVertex += distance;
                ++totalRuns;

                // Optional: Print the best distance for each run
                std::cout << "Run " << (ga + 1) << " on graph " << (g + 1)
                          << " for " << v << " nodes: Distance = " << distance << "\n";
                
                // Generate WOC tree using MST-based strategy (only for first run)
                if (!showGUI) {
                    WOCAlgorithm wocAlgorithm;
                    bestWOCTree = wocAlgorithm.generateSteinerTree(Strategy::MSTBased);
                    bestSteinerTree = best[0];
                    showGUI = true;
                }
            }
        }

        // Calculate the average distance for this node count
        double averageDistance = totalDistanceForVertex / totalRuns;
        averageDistances.push_back(averageDistance);

        // Output the average distance for this node count
        std::cout << "\nAverage distance for " << v << " nodes over " << totalRuns
                  << " runs: " << averageDistance << "\n\n";
    }

    // Optionally, output all average distances at the end
    std::cout << "Summary of average distances:\n";
    int nodeCount = 20;
    for (double avgDist : averageDistances) {
        std::cout << "Nodes: " << nodeCount << ", Average Distance: " << avgDist << "\n";
        nodeCount += 20;
    }

    // Display the GUI with the best trees
    if (showGUI) {
        std::cout << "\nOpening visualization window...\n";
        graphVisualizer visualizer;
        visualizer.displayPlot(bestSteinerTree, bestWOCTree);
    }

    return 0;
}

// main.cpp

#include "gui/graphVisualizer.h"
#include "algorithms/geneticAlgorithm.h"
#include "WOCAlgorithm.h"
#include "utilities/randomWeightedGraph.h"
#include "utilities/globalData.h"
#include "common.h"
#include <iostream>

int main() {
    Timer timer;

    // Parameters for the genetic algorithm
    int populationSize;
    int generations = 30;

    // Vector to store average distances for each vertex count
    std::vector<double> averageDistances;

    // Loop over vertex counts: 20, 40, 60, 80, 100
    for (int v = 20; v <= 200; v += 20) {
        double totalDistanceForVertex = 0.0;
        int totalRuns = 0;

        // Set the population size equal to the number of vertices
        populationSize = v;

        // Run 5 different random graphs for each vertex count
        for (int g = 0; g < 5; ++g) {
            // Set the global variable 'n' to the current number of vertices
            n = v;

            // Generate a random weighted graph
            randomWeightedGraph graph;
            graph.generateGraph();

            // Optional: Output the number of terminals if needed
            std::cout << "\nNumber of vertices: " << n << "\n";
            std::cout << "Number of terminal vertices: " << numTerminals << "\n";

            // Run the genetic algorithm 5 times on the same graph
            for (int ga = 0; ga < 3; ++ga) {
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
                          << " for " << v << " vertices: Distance = " << distance << "\n";
            }
        }

        // Calculate the average distance for this vertex count
        double averageDistance = totalDistanceForVertex / totalRuns;
        averageDistances.push_back(averageDistance);

        // Output the average distance for this vertex count
        std::cout << "\nAverage distance for " << v << " vertices over " << totalRuns
                  << " runs: " << averageDistance << "\n\n";
    }

    // Optionally, output all average distances at the end
    std::cout << "Summary of average distances:\n";
    int vertexCount = 20;
    for (double avgDist : averageDistances) {
        std::cout << "Vertices: " << vertexCount << ", Average Distance: " << avgDist << "\n";
        vertexCount += 20;
    }

    return 0;
}

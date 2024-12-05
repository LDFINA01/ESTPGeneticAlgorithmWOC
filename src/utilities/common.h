//
// Created by LuisF
//

#ifndef ESTPGENETICALGORITHMWOC_COMMON_H
#define ESTPGENETICALGORITHMWOC_COMMON_H

#include <chrono>
#include <vector>
#include "globalData.h"

class Timer {
public:
    Timer(); // Constructor to initialize the timer

    void start(); // Start the timer
    void stop();  // Stop the timer
    double elapsedMilliseconds() const; // Get elapsed time in milliseconds
    double elapsedSeconds() const;      // Get elapsed time in seconds

private:
    std::chrono::time_point<std::chrono::high_resolution_clock> startTime; // Start time point
    std::chrono::time_point<std::chrono::high_resolution_clock> endTime;   // End time point
    bool running; // Indicates if the timer is currently running
};

void printVector(const std::vector<int>& vec);
void printSteinerTree(const SteinerTree& tree);
void printEdges();

double computeTotalDistance(const std::vector<int>& path);
double computeTotalDistanceOfPopulation(const std::vector<SteinerTree>& population);
void incrementEdgeFrequencies(const std::vector<SteinerTree>& trees);
bool edgeExists(int u, int v);
std::vector<int> findPath(int start, int end);
int findClosestNode(int terminal, const std::vector<int>& path);
std::vector<int> repairPath(const std::vector<int>& path);


std::vector<SteinerTree> selectTopK(std::vector<SteinerTree> &population);


#endif //ESTPGENETICALGORITHMWOC_COMMON_H

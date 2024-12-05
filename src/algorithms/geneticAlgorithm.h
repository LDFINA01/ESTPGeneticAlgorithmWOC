// geneticAlgorithm.h
// Created by LuisF on 11/14/2024.

#ifndef ESTPGENETICALGORITHMWOC_GENETICALGORITHM_H
#define ESTPGENETICALGORITHMWOC_GENETICALGORITHM_H

#include <vector>
#include <unordered_map>
#include "globalData.h"


    SteinerTree randomSteinerTree();
    std::vector<int> removeConsecutiveDuplicates(const std::vector<int>& path);
    std::vector<SteinerTree> initialPopulation (int size);
    SteinerTree crossoverFunction(const SteinerTree &parentA, const SteinerTree &parentB);
    void randomMutation(std::vector<SteinerTree>& population, double mutationRate);
    std::vector<SteinerTree> runGeneticAlgorithm(int populationSize, int generations);



#endif //ESTPGENETICALGORITHMWOC_GENETICALGORITHM_H

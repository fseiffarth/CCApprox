//
// Created by Florian on 08.05.2021.
//

#ifndef CLOSURES_GENERATIONTEST_H
#define CLOSURES_GENERATIONTEST_H

#include "../Utils/Generators.h"

void generationTest(){
    std::mt19937_64 gen;
    std::vector<GraphData> graphs;
    Generators::generateTwoComponentGraphs("../out/graphs/", graphs, Properties(5, 1.3, 5, 1.3, 3, "Graph"), 10, gen);
}

#endif //CLOSURES_GENERATIONTEST_H

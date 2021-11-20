//
// Created by Florian on 09.05.2021.
//

#ifndef CLOSURES_EVALUATIONS_H
#define CLOSURES_EVALUATIONS_H


#include "../Data/GraphData.h"
#include "../Utils/graph_generators.h"

class Evaluations {
public:
    Evaluations(const Evaluations& evaluations);
    Evaluations(const GraphData &graphData, std::string closureOperator, std::string iuStrategy, double blockDensity, int blockConnection, int numGraphs, int trainingConfigurations, std::string patternType);
    void calculate_accuracies(GraphData& graphData, const Labels& prediction, std::vector<double>& acc);
    void evaluate(GraphData& graphData, const std::vector<Labels>& predictions, int eval_stepsize, std::vector<double>& runtimes);
    void save(const std::string& path);
    void copy(Evaluations &evaluation);

    std::vector<std::vector<double>> accuracies;
    std::vector<std::vector<double>> _runtimes;
    std::vector<int> iterations;
    int trainingSize = 0;
    int correctlyClassified = 0;
    int testSize = 0;
    int numGraphs = 0;
    int trainingConfigurations = 0;
    double blockDensity = 0;
    int blockConnections = 0;
    std::string patternType;
    GenerationType generationType;
    std::string _name;
    int _dataSize;
    int _dataEdges;
    std::string iuStrategy;
    std::string closureOperator;



};


#endif //CLOSURES_EVALUATIONS_H

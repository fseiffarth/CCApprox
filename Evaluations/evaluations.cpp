//
// Created by Florian on 09.05.2021.
//

#include "evaluations.h"

#include <utility>


Evaluations::Evaluations(const GraphData &graphData, std::string closureOperator, std::string iuStrategy, double blockDensity, int blockConnection, int numGraphs, const int trainingConfigurations, std::string patternType)
        : blockDensity(blockDensity), blockConnections(blockConnection), numGraphs(numGraphs), patternType(std::move(patternType)),
          trainingConfigurations(trainingConfigurations), _name(graphData.getName()), _dataSize(graphData.size()), _dataEdges(
                graphData.get_graph()->GetEdges()),
          iuStrategy(iuStrategy), closureOperator(closureOperator) {
}

Evaluations::Evaluations(const Evaluations& evaluations) : blockDensity(evaluations.blockDensity), blockConnections(evaluations.blockConnections), iterations(evaluations.iterations), numGraphs(evaluations.numGraphs), patternType(evaluations.patternType),
                                                        trainingConfigurations(evaluations.trainingConfigurations), _name(evaluations._name), _dataSize(evaluations._dataSize), _dataEdges(evaluations._dataEdges), accuracies(evaluations.accuracies), _runtimes(evaluations._runtimes), trainingSize(evaluations.trainingSize), correctlyClassified(evaluations.correctlyClassified), testSize(evaluations.testSize), iuStrategy(evaluations.iuStrategy), closureOperator(evaluations.closureOperator) {
}

void Evaluations::calculate_accuracies(GraphData& graphData, const Labels& prediction, std::vector<double>& acc) {
    this->correctlyClassified = 0;
    this->trainingSize = 0;
    for (const std::set<NodeId>& trainingSet: graphData.trainingSet) {
        this->trainingSize += static_cast<int>(trainingSet.size());
    }
    this->testSize = static_cast<int>(graphData.size()) - this->trainingSize;
    int Id = 0;
    for (Label label : graphData.labels()) {
        if (prediction[Id] == label){
            ++this->correctlyClassified;
        }
        ++Id;
    }
    acc.emplace_back((double) (this->correctlyClassified - this->trainingSize) / (double) this->testSize);
}

void Evaluations::evaluate(GraphData &graphData, const std::vector<Labels> &predictions, int eval_stepsize, std::vector<double>& runtimes) {
    int max_iterations = static_cast<int>(predictions.size());
    this->_runtimes.emplace_back(runtimes);
    int counter = 0;
    std::vector<double> acc;
    for (int iter = 1; iter <= max_iterations; iter += eval_stepsize) {
        if (this->iterations.size() < (int) max_iterations / eval_stepsize){
            this->iterations.push_back(iter);
        }
        std::vector<std::vector<int>> labelFrequencies = std::vector<std::vector<int>>(graphData.size(), std::vector<int>(graphData.class_num(), 0));
        std::vector<Label> prediction;
        for (const Labels& labels : std::vector<Labels>(predictions.begin(), predictions.begin() + iter)) {
            for (int i = 0; i < labels.size(); ++i) {
                if (labels[i] >=0) {
                    ++labelFrequencies[i][labels[i]];
                }
            }
        }
        for (auto & labelFrequency : labelFrequencies) {
            int maxLabel = 0;
            int maxFrequency = 0;
            for (int j = 0; j < labelFrequency.size(); ++j) {
                if (labelFrequency[j] > maxFrequency){
                    maxLabel = j;
                    maxFrequency = labelFrequency[j];
                }
            }
            prediction.emplace_back(maxLabel);
        }
        calculate_accuracies(graphData, prediction, acc);
    }
    this->accuracies.emplace_back(acc);
}


void Evaluations::save(const std::string& path) {
    if (static_cast<int>(iterations.size()) > 0) {
        bool newFile = std::filesystem::exists(path);
        std::ofstream fs;
        fs.open(path, std::ios_base::app);
        // write the file headers
        if (!newFile) {
            fs << ",Graphs" << "," << "GraphSize" << "," << "GraphEdgeNum" << "," << "BlockDensity" << "," << "BlockConnections" << "," << "PatternType" << ","
            << "TrainingSize" << "," << "TrainingConfigurations" << "," << "Accuracy" << "," << "Deviations" << "," << "NumGraphs" << ","
            << "Iterations" << "," << "Closure" << "," << "InitUpdateStrategy" << "," << "Generation" << "," << "Single Runtime Mean" << std::endl;
        }
        std::string EstimationString;
        std::string generationTypeString;
        switch (this->generationType) {
            case GenerationType::ONE_COMPONENT:
                generationTypeString = "OneComponent";
                break;
                case GenerationType::TWO_COMPONENTS:
                    generationTypeString = "TwoComponents";
                    break;
        }

        for (int i = 0; i < iterations.size(); ++i){
            std::vector<double> acc_per_iteration;
            std::vector<double> runtime_per_iteration;

            for (int num = 0; num < accuracies.size(); ++num) {
                acc_per_iteration.emplace_back(accuracies[num][i]);
                runtime_per_iteration.emplace_back(_runtimes[num][i]);
            }
            double acc_sum = std::accumulate(acc_per_iteration.begin(), acc_per_iteration.end(), 0.0);
            double acc_mean = acc_sum / static_cast<int>(acc_per_iteration.size());
            double acc_sq_sum = std::inner_product(acc_per_iteration.begin(), acc_per_iteration.end(), acc_per_iteration.begin(), 0.0);
            double acc_std = std::sqrt(acc_sq_sum / static_cast<int>(acc_per_iteration.size()) - acc_mean * acc_mean);

            double run_sum = std::accumulate(runtime_per_iteration.begin(), runtime_per_iteration.end(), 0.0);
            double run_mean = (run_sum / static_cast<int>(runtime_per_iteration.size())) / 1000000;

            fs << std::fixed << "," << this->_name << "," << this->_dataSize << ","
            << this->_dataEdges << "," << blockDensity << "," << blockConnections << "," << patternType << "," << trainingSize << ","
            << trainingConfigurations << "," << acc_mean << "," << acc_std << "," << numGraphs << "," << iterations[i] << "," << closureOperator << "," << iuStrategy << "," << generationTypeString << "," << run_mean << std::endl;
            fs << std::scientific;
        }
        fs.close();
    }
}



void Evaluations::copy(Evaluations &evaluation) {
    accuracies = evaluation.accuracies;
    _runtimes = evaluation._runtimes;
    iterations = evaluation.iterations;
    trainingSize = evaluation.trainingSize;
    correctlyClassified = evaluation.correctlyClassified;
    testSize = evaluation.testSize;
    numGraphs = evaluation.numGraphs;
    trainingConfigurations = evaluation.trainingConfigurations;
    blockDensity = evaluation.blockDensity;
    blockConnections = evaluation.blockConnections;
    patternType = evaluation.patternType;
    generationType = evaluation.generationType;

    _name = evaluation._name;
    _dataSize = evaluation._dataSize;
    _dataEdges = evaluation._dataEdges;
    iuStrategy = evaluation.iuStrategy;
    closureOperator = evaluation.closureOperator;

}



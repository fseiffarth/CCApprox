//
// Created by florian on 31.08.21.
//


#include "../../Data/GraphData.h"
#include "../../Experiments/Experiments.h"
#include "../../Utils/FileEvaluation.h"
#include "../../Experiments/DS2021/SetupExperiment.h"
#include "../../Utils/StaticFunctions.h"

void closure_runtime(int number) {
    int n = 1000;
    std::string out_name = "closure_runtime";

    std::vector<std::string> headers = {"Size",
                                        "Edges",
                                        "Density",
                                        "Input Size",
                                        "Outerplanar Generation Time",
                                        "Outerplanar Time",
                                        "Outerplanar New Time",
                                        "Tree Time",
                                        "Graph Time",
                                        "Samples"};

#pragma omp parallel for default(none) shared(n, headers, number, out_name)
    for (int j = 0; j < 8; ++j) {
        double p = 0.006 + 0.002*j;
        std::vector<std::pair<int, int>> graph_sizes;
        for (int i = 1; i <= 10; ++i) {
            graph_sizes.emplace_back(std::pair<int, int>{i * n, (int) (i * n * (i * n - 1) / 2 * p)});
        }
        for (auto const &size: graph_sizes) {
            ClosureParameters closureParameters;
            GraphClosureSP cl;
            std::vector<GraphData> outerplanar;
            std::vector<GraphData> graphs;
            std::vector<GraphData> trees;
            size_t outerplanar_time = 0;
            size_t outerplanar_generation_time = 0;
            size_t outerplanar_new_time = 0;
            size_t tree_time = 0;
            size_t graph_time = 0;
            int input_size = size.first / 100;
            FileEvaluation fileEvaluation = FileEvaluation("../out/ICDE_2022Exp/" + out_name);
            OuterplanarGraphStatistics statistics;
            std::vector<std::string> values;
            outerplanar.clear();
            graphs.clear();
            trees.clear();
            Experiments::GetOuterplanarSamples(size, number, outerplanar, graphs, trees);
            int edges = 0;
            for (int i = 0; i < outerplanar.size(); ++i) {
                statistics += OuterplanarGraphStatistics(outerplanar[i].get_graph());
                std::set<NodeId> input_set;
                StaticFunctions::generateInputSet(closureParameters.input_set, outerplanar[i], input_size, i);
                edges += outerplanar[i].edges();

                auto start = std::chrono::high_resolution_clock::now();
                cl.naive_closure(outerplanar[i], closureParameters);
                outerplanar_time += std::chrono::duration_cast<std::chrono::microseconds>(
                        std::chrono::high_resolution_clock::now() - start).count();

                start = std::chrono::high_resolution_clock::now();
                OuterplanarGraphData outerplanar_new = OuterplanarGraphData(outerplanar[i]);
                outerplanar_generation_time += std::chrono::duration_cast<std::chrono::microseconds>(
                        std::chrono::high_resolution_clock::now() - start).count();

                start = std::chrono::high_resolution_clock::now();
                cl.naive_closure(outerplanar_new, closureParameters);
                outerplanar_new_time += std::chrono::duration_cast<std::chrono::microseconds>(
                        std::chrono::high_resolution_clock::now() - start).count();

                start = std::chrono::high_resolution_clock::now();
                cl.naive_closure(graphs[i], closureParameters);
                graph_time += std::chrono::duration_cast<std::chrono::microseconds>(
                        std::chrono::high_resolution_clock::now() - start).count();

                start = std::chrono::high_resolution_clock::now();
                cl.naive_closure(trees[i], closureParameters);
                tree_time += std::chrono::duration_cast<std::chrono::microseconds>(
                        std::chrono::high_resolution_clock::now() - start).count();

                //std::cout << outerplanar[i].nodes() << " " << outerplanar[i].edges() << std::endl;
                //std::cout << graphs[i].nodes() << " " << graphs[i].edges() << std::endl;
                //std::cout << trees[i].nodes() << " " << trees[i].edges() << std::endl;
            }
            std::vector<std::string> stat_headers;
            std::vector<std::string> stat_values;
            statistics.evaluate(stat_headers, stat_values);
            edges /= (int) outerplanar.size();
            fileEvaluation.headerValueInsert(headers, {std::to_string(size.first),
                                                       std::to_string(edges),
                                                       std::to_string(size.second / ((double) size.first/2 * (size.first - 1))),
                                                       std::to_string(input_size),
                                                       std::to_string((double) outerplanar_generation_time/1000000.0),
                                                       std::to_string((double) outerplanar_time/1000000.0),
                                                       std::to_string((double) outerplanar_new_time/1000000.0),
                                                       std::to_string((double) tree_time/1000000.0),
                                                       std::to_string((double) graph_time/1000000.0),
                                                       std::to_string(outerplanar.size())});
            fileEvaluation.headerValueInsert(stat_headers, stat_values);
#pragma omp critical
            fileEvaluation.save();
        }
    }
}


int main(){
    int number = 100;
    closure_runtime(number);
}
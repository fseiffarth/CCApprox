#include <iostream>
#include "../Data/GraphData.h"
#include "../Utils/pattern_generator.h"
#include "../Utils/StaticFunctions.h"
#include "../Data/SimpleGraph.h"

//
// Created by florian on 07.09.21.
//
int main() {
    std::string out_path = "../out/Experiments/SamplingRuntime/sampling_runtime.csv";
    int iterations = 120;
    int seed = 0;
    int threads = 8;
    std::vector<PatternType> patternTypes = {PatternType::BFS_TREE, PatternType::OUTERPLANAR};
    std::vector<std::string> graph_paths = {"../../GraphData/roadNet-CA_seed_0_size_199300.edges"};
    std::vector<GraphData> graphs;
    if (threads == -1){
        threads = omp_get_max_threads();
    }
    omp_set_num_threads(threads);
    for (auto const &graph_path: graph_paths) {
        graphs.clear();
        graphs.emplace_back(GraphData(graph_path));
    }
    for (auto const& graph : graphs) {
        if (TSnap::IsConnected(graphs.back().get_graph())) {

            std::string type;
            //Get graph max degree
            int maxDegree = 0;
            for (auto nodeI = graph.get_graph()->BegNI(); nodeI != graph.get_graph()->EndNI(); nodeI++) {
                maxDegree = std::max(maxDegree, nodeI.GetDeg());
            }
            std::vector<NodeId> neighborIds(maxDegree);
            std::iota(neighborIds.begin(), neighborIds.end(), 0);

            std::cout << "Graph " << graphs.back().getName() << " is connected!" << std::endl;

            int size = graph.get_graph()->GetNodes();
            int num_edges = graph.get_graph()->GetEdges();
            OuterplanarSubgraphDFS outerPlanarSubgraphDfs = OuterplanarSubgraphDFS(graph.get_graph());
            for (auto const patternType: patternTypes) {
                double approximation_edges = 0;
                double generation_time = 0;
                OuterplanarGraphStatistics statistics = OuterplanarGraphStatistics();
#pragma omp parallel for default(none) shared(iterations, seed, patternType, type, graph, size) firstprivate(neighborIds, outerPlanarSubgraphDfs) private(statistics) reduction(+:generation_time, approximation_edges)
                    for (int i = 0; i < iterations; ++i) {
                        GraphData subgraphData = GraphData(new TUNGraph(), size);
                        std::mt19937_64 generator(i + iterations * seed);
                        auto start_generation = std::chrono::high_resolution_clock::now();
                        switch (patternType) {
                            case PatternType::BFS_TREE:
                                type = "Tree Approximation";
                                start_generation = std::chrono::high_resolution_clock::now();
                                GraphFunctions::bfsSubtree(graph.get_graph(), subgraphData.graph(), neighborIds,
                                                           generator);
                                subgraphData.graphType = GraphType::TREE;
                                generation_time += ((double) std::chrono::duration_cast<std::chrono::microseconds>(
                                        std::chrono::high_resolution_clock::now() - start_generation).count() /
                                                    1000000.0);
                                break;
                            case PatternType::OUTERPLANAR:
                                type = "Outerplanar Approximation";
                                start_generation = std::chrono::high_resolution_clock::now();
                                outerPlanarSubgraphDfs.generate(subgraphData.graph(), generator, false);
                                subgraphData.graphType = GraphType::OUTERPLANAR;
                                generation_time += ((double) std::chrono::duration_cast<std::chrono::microseconds>(
                                        std::chrono::high_resolution_clock::now() - start_generation).count() /
                                                    1000000.0);
                                statistics += OuterplanarGraphStatistics(subgraphData.get_graph());
                                break;
                        }
                        approximation_edges += subgraphData.get_graph()->GetEdges();
                    }

                statistics /= iterations;
                std::vector<std::string> stat_headers;
                std::vector<std::string> stat_values;
                statistics.evaluate(stat_headers, stat_values);
                std::vector<std::string> headers = {"Graph Name", "Graph Size",
                                                        "Num Edges",
                                                        "Approximation Type",
                                                        "Approximation Edges",
                                                        "Iterations",
                                                        "Runtime",
                                                        "Threads"};

                std::vector<std::string> values = {graph.getName(),
                                                       std::to_string(size),
                                                       std::to_string(num_edges),
                                                       type,
                                                       std::to_string((double) approximation_edges / iterations),
                                                       std::to_string(iterations),
                                                       std::to_string(generation_time / threads),
                                                       std::to_string(threads)};

                headers.insert(headers.end(), stat_headers.begin(), stat_headers.end());
                values.insert(values.end(), stat_values.begin(), stat_values.end());
                StaticFunctions::saveValuesToFile(out_path, headers, values);
            }
        }
    }
}
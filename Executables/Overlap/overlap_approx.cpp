//
// Created by florian on 30.09.21.
//

#include <iostream>
#include "../../Data/GraphData.h"
#include "../../ClosureOperators/GraphClosures.h"
#include "../../Utils/StaticFunctions.h"
#include "../../Experiments/DS2021/SetupExperiment.h"

// Comparison function to sort the vector elements
// by second element of tuples
bool sortbysecond(const std::tuple<std::string, int>& a,
               const std::tuple<std::string, int>& b)
{
    return (std::get<1>(a) < std::get<1>(b));
}

struct OverlapApproxParams{
    //Paths
    std::string in_path;
    std::string out_path;

    int thread_num = -1;

    int max_nodes = 100000;
    int max_edges = 200000;

    std::vector<int> generator_size = {10};
    int generator_seed = 0;
    int coreIterations = 3;
    int sample_number = 200;
    std::vector<double> threshold = {1.0 / (double) this->sample_number};
    double overall_threshold = 1.0 / (double) this->sample_number;
    bool overall = false;
    int sample_seed = 32487643;
    bool outerplanar_new = true;
    
    //Eval
    bool periphery = false;
    bool core_iteration = false;
    bool exact = true;
    bool approx_core = false;

    bool simple_approx = false;
    int runtime = 0;
    int simple_approx_iterations = 0;

    bool core_stats = false;
    bool outerplanar_statistics = false;
    bool threshold_percentage = false;
    bool tree_eval = true;
};


void get_core(GraphData& graph, int generator_size, int core_iterations, int seed, std::map<int, int>& degree_distribution, TIntV& coreNodes, double& runtime){
    GraphClosureSP gc;
    ClosureParameters closureParameters;
    std::vector<std::set<NodeId>> closures;
    std::cout << std::endl;
    std::set<NodeId> overlap;
    runtime = 0;
    std::chrono::time_point<std::chrono::system_clock> start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < core_iterations; ++i) {
        StaticFunctions::generateInputSet(closureParameters.input_set, graph, generator_size, i + core_iterations * seed);
        gc.naive_closure(graph, closureParameters);
        std::cout << "\tClosure Size: " << closureParameters.closed_set.size() << std::endl;
        closures.emplace_back(closureParameters.closed_set);
        if (i == 0) {
            overlap = closureParameters.closed_set;
        } else {
            std::vector<NodeId> v_intersection;
            std::set_intersection(overlap.begin(), overlap.end(),
                                  closures.back().begin(), closures.back().end(),
                                  std::back_inserter(v_intersection));
            overlap.clear();
            for (auto elem: v_intersection) {
                overlap.insert(elem);
            }
        }
        std::cout << "\tOverlap Size: " << overlap.size() << std::endl;
    }
    coreNodes.Clr();
    for (auto elem: overlap) {
        coreNodes.Add(elem);
    }
    //Measure runtime before getting statistics
    runtime = ((double) std::chrono::duration_cast<std::chrono::microseconds>(
            std::chrono::high_resolution_clock::now() - start).count() /
               1000000.0);

}

void simple_approx(GraphData &graph, int generator_size, int core_iterations, int seed, std::map<int, int>& degree_distribution,
                   TIntV& coreNodes, double& runtime, double time_constraint = -1, int iteration_constraint = 0) {
    GraphClosureSP gc;
    ClosureParameters closureParameters;
    std::vector<NodeId> v_intersection;
    std::cout << std::endl;
    std::set<NodeId> overlap;
    runtime = 0;
    std::chrono::time_point<std::chrono::system_clock> start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < core_iterations; ++i) {
        StaticFunctions::generateInputSet(closureParameters.input_set, graph, generator_size, i + core_iterations * seed);
        closureParameters.iteration_number = iteration_constraint;
        closureParameters.timeConstraint = time_constraint/core_iterations;
        gc.naive_closure(graph, closureParameters);
        std::cout << "\tClosure Size: " << closureParameters.closed_set.size() << std::endl;
        if (i == 0) {
            overlap = closureParameters.closed_set;
        } else {
            v_intersection.clear();
            std::set_intersection(overlap.begin(), overlap.end(),
                                  closureParameters.closed_set.begin(), closureParameters.closed_set.end(),
                                  std::back_inserter(v_intersection));
            overlap.clear();
            if (i == core_iterations -1) {
                closureParameters.closed_set.clear();
                for (auto elem: v_intersection) {
                    closureParameters.closed_set.insert(elem);
                }
            }
            else {
                for (auto elem: v_intersection) {
                    overlap.insert(elem);
                }
            }
        }
        std::cout << "\tOverlap Size: " << closureParameters.closed_set.size() << std::endl;
    }
    coreNodes.Clr();
    for (auto elem: closureParameters.closed_set) {
        coreNodes.Add(elem);
    }
    //Measure runtime before getting statistics
    runtime = ((double) std::chrono::duration_cast<std::chrono::microseconds>(
            std::chrono::high_resolution_clock::now() - start).count() /
               1000000.0);
}

void approximate_core(GraphData& graph, std::vector<GraphData>& samples, int generator_size, double threshold, OverlapApproxParams& params, std::map<int, int>& degree_distribution, TIntV& coreNodes, double& runtime){
    GraphClosureSP gc;
    ClosureParameters closureParameters;
    std::set<NodeId> overlap;
    std::vector<int> overall_approximation = std::vector<int>(graph.size(), 0);
    std::vector<NodeId> v_intersection;
    double closure_time = 0;
    runtime = 0;
    std::chrono::time_point<std::chrono::system_clock> start = std::chrono::high_resolution_clock::now();
    std::vector<int> approximation = std::vector<int>(graph.size(), 0);
    for (int i = 0; i < params.coreIterations; ++i) {
        StaticFunctions::generateInputSet(closureParameters.input_set, graph, generator_size, i + params.coreIterations * params.generator_seed);
        gc.approx_closure(samples, closureParameters, approximation, closure_time);
        for(int j = 0; j < approximation.size(); ++j){
            overall_approximation[j] += approximation[j];
        }
        if (!params.overall) {
            SetupExperiment::getClosedFromApproximation(approximation, closureParameters.closed_set,
                                                        (int) samples.size(), threshold);
            std::cout << "\tClosure Size: " << closureParameters.closed_set.size() << std::endl;
            if (i == 0) {
                overlap = closureParameters.closed_set;
            } else {
                v_intersection.clear();
                std::set_intersection(overlap.begin(), overlap.end(),
                                      closureParameters.closed_set.begin(), closureParameters.closed_set.end(),
                                      std::back_inserter(v_intersection));
                overlap.clear();
                if (i == params.coreIterations -1) {
                    closureParameters.closed_set.clear();
                    for (auto elem: v_intersection) {
                        closureParameters.closed_set.insert(elem);
                    }
                }
                else {
                    for (auto elem: v_intersection) {
                        overlap.insert(elem);
                    }
                }
            }
            std::cout << "\tOverlap Size: " << closureParameters.closed_set.size() << std::endl;
        }
    }
    if (params.overall){
        SetupExperiment::getClosedFromApproximation(overall_approximation, closureParameters.closed_set, (int) samples.size(), params.overall_threshold);
        std::cout << "\tOverlap Size: " << closureParameters.closed_set.size() << std::endl;
    }
    coreNodes.Clr();
    for (auto elem: closureParameters.closed_set) {
        coreNodes.Add(elem);
    }
    //Measure runtime before getting statistics
    runtime = ((double) std::chrono::duration_cast<std::chrono::microseconds>(
            std::chrono::high_resolution_clock::now() - start).count() /
               1000000.0);
}

void approximate_core(GraphData& graph, std::vector<OuterplanarGraphData>& samples, int generator_size, double threshold, OverlapApproxParams& params, std::map<int, int>& degree_distribution, TIntV& coreNodes, double& runtime){
    GraphClosureSP gc;
    ClosureParameters closureParameters;
    std::set<NodeId> overlap;
    std::vector<int> overall_approximation = std::vector<int>(graph.size(), 0);
    std::vector<NodeId> v_intersection;
    double closure_time = 0;
    runtime = 0;
    std::chrono::time_point<std::chrono::system_clock> start = std::chrono::high_resolution_clock::now();
    std::vector<int> approximation = std::vector<int>(graph.size(), 0);
    for (int i = 0; i < params.coreIterations; ++i) {
        StaticFunctions::generateInputSet(closureParameters.input_set, graph, generator_size, i + params.coreIterations * params.generator_seed);
        gc.approx_closure(samples, closureParameters, approximation, closure_time);
        for(int j = 0; j < approximation.size(); ++j){
            overall_approximation[j] += approximation[j];
        }
        if (!params.overall) {
            SetupExperiment::getClosedFromApproximation(approximation, closureParameters.closed_set,
                                                        (int) samples.size(), threshold);
            std::cout << "\tClosure Size: " << closureParameters.closed_set.size() << std::endl;
            if (i == 0) {
                overlap = closureParameters.closed_set;
            } else {
                v_intersection.clear();
                std::set_intersection(overlap.begin(), overlap.end(),
                                      closureParameters.closed_set.begin(), closureParameters.closed_set.end(),
                                      std::back_inserter(v_intersection));
                std::cout << "\tOverlap Size: " << v_intersection.size() << std::endl;
                overlap.clear();
                if (i == params.coreIterations -1) {
                    closureParameters.closed_set.clear();
                    for (auto elem: v_intersection) {
                        closureParameters.closed_set.insert(elem);
                    }
                }
                else {
                    for (auto elem: v_intersection) {
                        overlap.insert(elem);
                    }
                }
            }
        }
    }
    if (params.overall){
        SetupExperiment::getClosedFromApproximation(overall_approximation, closureParameters.closed_set, (int) samples.size(), params.overall_threshold);
        std::cout << "\tOverlap Size: " << closureParameters.closed_set.size() << std::endl;
    }
    coreNodes.Clr();
    for (auto elem: closureParameters.closed_set) {
        coreNodes.Add(elem);
    }
    //Measure runtime before getting statistics
    runtime = ((double) std::chrono::duration_cast<std::chrono::microseconds>(
            std::chrono::high_resolution_clock::now() - start).count() /
               1000000.0);
}



void overlap_eval(OverlapApproxParams& params) {
    std::vector<std::tuple<std::string, int>> paths_and_sizes;
    for (const auto &entry: std::filesystem::recursive_directory_iterator(params.in_path)) {
        std::string stripped_path =
                entry.path().parent_path().string() + "/" + entry.path().stem().string();
        if (entry.path().extension() == ".core" || (params.exact = false && entry.path().extension() == ".edges")){
            const GraphData& data = GraphData(stripped_path + ".edges");
            int nodes = data.nodes();
            int edges = data.edges();

            if (nodes < params.max_nodes && edges < params.max_edges) {
                paths_and_sizes.emplace_back(std::tuple<std::string, int>(stripped_path + ".edges", edges));
            }
            data.get_graph().Clr();
        }
    }
    std::sort(paths_and_sizes.begin(), paths_and_sizes.end(), sortbysecond);

    omp_set_num_threads(params.thread_num);
#pragma omp parallel for default(none) shared(paths_and_sizes, params)
    for (const auto &[path, size]: paths_and_sizes) {
        std::string stripped_path =
                std::filesystem::path(path).parent_path().string() + "/" + std::filesystem::path(path).stem().string();
        TIntV exact_core;
        StaticFunctions::load(stripped_path + ".core", exact_core);
        bool exact_computation;
        if (exact_core.Len() > 0) {
            exact_computation = true;
        } else {
            exact_computation = false;
        }
        FileEvaluation eval = FileEvaluation(params.out_path, (int) (params.generator_size.size()*params.threshold.size()), ".csv");

        std::vector<std::vector<std::string>> info_array;
        StaticFunctions::load_csv(stripped_path + ".core_info", info_array);
        double exact_runtime = std::stoi(info_array[1][9]);
        GraphData graph = GraphData(path), core_graph;
        //std::cout << "Graph Nodes: " << graph.nodes() << std::endl;
        GraphFunctions::GetLargestComponent(graph);
        //GraphFunctions::analyse_graph(graph.get_graph(), graph.getName(), false, &eval);
        eval.headerValueInsert({"Graph", "Nodes", "Edges", "Density", "Core Size", "Exact Runtime"},
                               {graph.getName(), std::to_string(graph.nodes()), std::to_string(graph.edges()),
                                std::to_string(graph.density()),
                                std::to_string(exact_core.Len()),
                                std::to_string(exact_runtime)},-1, true, true);

        int counter = 0;
        for (double threshold : params.threshold) {
            for (int generator_size: params.generator_size) {
                eval.headerValueInsert({"Threshold", "Generator Size"},
                                       {std::to_string(threshold),
                                        std::to_string(generator_size)}, counter, true, true);
                ++counter;
            }
        }
        std::map<int, int> core_degree_distribution;
        TIntV core_nodes;
        if (exact_computation) {
            GraphFunctions::GetCoreGraphStats(core_graph, graph, exact_core, exact_runtime, "Exact Core", &eval);
            //GraphFunctions::analyse_graph(core_graph.get_graph(), "Core");
//                if (params.periphery) {
//                    PUNGraph periphery_graph = GraphFunctions::RemoveSubgraph(graph.get_graph(),
//                                                                              core_graph.get_graph());
//                    GraphFunctions::analyse_graph(periphery_graph, "Periphery");
//                }
//                std::cout << "Exact Runtime: " << exact_runtime << "s" << std::endl;
//                if (params.core_iteration) {
//                    GraphData core2_graph;
//                    core_graph.set_graph(GraphFunctions::ResetGraphIds(core_graph.graph()));
//                    get_core(core_graph, params.generator_size, params.coreIterations, params.generator_seed,
//                             core_degree_distribution, exact_core, exact_runtime);
//                }
        }

        double tree_sampling_runtime = 0, outerplanar_sampling_runtime = 0;
        double tree_core_runtime = 0, outerplanar_core_runtime = 0;
        double conversion_runtime = 0;
        std::vector<std::string> out_stat_headers;
        std::vector<std::string> out_stat_values;
        if (params.approx_core) {
            std::vector<int> variation_check;
            std::vector<NodeId> neighborIds;
            GraphFunctions::generateNeighborVector(graph.get_graph(), neighborIds);
            OuterplanarGraphStatistics statistic;

            //Approximate cores
            double overall_tree = 0, overall_outerplanar = 0;
            for (auto type: {PatternType::BFS_TREE, PatternType::OUTERPLANAR}) {
                if (type == PatternType::BFS_TREE && params.tree_eval) {
                    //Tree samples
                    std::vector<GraphData> treeSamples;
                    GraphFunctions::GetSamples(graph, PatternType::BFS_TREE, treeSamples, nullptr, neighborIds,
                                               params.sample_number, params.sample_seed, tree_sampling_runtime);
                    //std::cout << "Tree Sampling Runtime: " << runtime << "s" << std::endl;

                    //std::cout << std::endl << "Approximate Tree Core:" << std::endl;
                    int counter = 0;
                    for (double threshold : params.threshold) {
                        for (int generator_size: params.generator_size) {
                            approximate_core(graph, treeSamples, generator_size, threshold, params, core_degree_distribution,
                                             core_nodes, tree_core_runtime);
                            if (params.core_stats) {
                                GraphFunctions::GetCoreGraphStats(core_graph, graph, core_nodes, tree_core_runtime,
                                                                  "Tree Core",
                                                                  &eval);
                            } else {
                                eval.headerValueInsert(
                                        {"Tree Core Nodes", "Tree Core Relative Nodes", "Tree Core Edges",
                                         "Tree Core Relative Edges", "Tree Core Out Edges"},
                                        {std::to_string(core_nodes.Len()),
                                         std::to_string((double)graph.nodes()/core_nodes.Len()),
                                         std::to_string(0),
                                         std::to_string(0),
                                         std::to_string(0)}, counter);
                            }
                            overall_tree = tree_sampling_runtime + tree_core_runtime;
                            //std::cout << "Tree Core Runtime: " << tree_core_runtime << "s" << std::endl;
                            //std::cout << "Tree Runtime: " << overall_tree << "s" << std::endl;
                            //GraphFunctions::analyse_graph(core_graph.get_graph(), "Approximation Core Tree");
                            eval.headerValueInsert({"Tree Core Runtime", "Tree Runtime"},
                                                   {std::to_string(tree_core_runtime),
                                                    std::to_string(overall_tree)}, counter, true,
                                                   true);
                            if (exact_computation) {
                                int intersection_length = core_nodes.IntrsLen(exact_core);
                                int union_length = core_nodes.UnionLen(exact_core);
                                double recall = (double) intersection_length / exact_core.Len();
                                double similarity = (double) intersection_length / union_length;
                                //std::cout << "Tree Recall: " << recall << std::endl;
                                //std::cout << "Tree Similarity: " << similarity << std::endl;
                                //Approx values
                                eval.headerValueInsert(
                                        {"Tree Recall", "Tree Similarity"},
                                        {
                                                std::to_string(recall),
                                                std::to_string(similarity)}, counter, true, true);
                            } else {
                                eval.headerValueInsert(
                                        {"Tree Recall", "Tree Similarity"},
                                        {
                                                std::to_string(0),
                                                std::to_string(0)}, counter, true, true);
                            }
                            ++counter;
                        }
                    }
                } else if (type == PatternType::OUTERPLANAR) {
                    std::vector<GraphData> outerplanarSamples;
                    std::vector<OuterplanarGraphData> outerplanarSamples_new;
                    //Outerplanar samples
                    OuterplanarSubgraphDFS outerPlanarSubgraphDfs = OuterplanarSubgraphDFS(graph.get_graph());
                    if (!params.outerplanar_new) {
                        GraphFunctions::GetSamples(graph, PatternType::OUTERPLANAR, outerplanarSamples,
                                                   &outerPlanarSubgraphDfs, neighborIds, params.sample_number,
                                                   params.sample_seed, outerplanar_sampling_runtime);
                        if (params.outerplanar_statistics) {
                            for (auto const &out_graph: outerplanarSamples) {
                                statistic += OuterplanarGraphStatistics(out_graph.get_graph());
                            }
                        }
                    } else {
                        GraphFunctions::GetOuterplanarSamples(graph, PatternType::OUTERPLANAR,
                                                              outerplanarSamples_new,
                                                              &outerPlanarSubgraphDfs, neighborIds,
                                                              params.sample_number,
                                                              params.sample_seed, outerplanar_sampling_runtime,
                                                              conversion_runtime, true);
                        if (params.outerplanar_statistics) {
                            for (auto const &out_graph: outerplanarSamples_new) {
                                statistic += OuterplanarGraphStatistics(out_graph.get_graph());
                            }
                        }
                    }

                    //std::cout << std::endl << "Approximate Outerplanar Core:" << std::endl;
                    int counter = 0;
                    for (double threshold: params.threshold) {
                        for (int generator_size: params.generator_size) {
                            if (!params.outerplanar_new) {
                                approximate_core(graph, outerplanarSamples, generator_size, threshold, params, core_degree_distribution,
                                                 core_nodes, outerplanar_core_runtime);
                            } else {
                                approximate_core(graph, outerplanarSamples_new, generator_size, threshold, params, core_degree_distribution,
                                                 core_nodes, outerplanar_core_runtime);
                            }
                            if (params.core_stats) {
                                GraphFunctions::GetCoreGraphStats(core_graph, graph, core_nodes,
                                                                  outerplanar_core_runtime,
                                                                  "Outerplanar Core", &eval);
                            } else {
                                eval.headerValueInsert({"Outerplanar Core Nodes", "Outerplanar Core Relative Nodes",
                                                        "Outerplanar Core Edges", "Outerplanar Core Relative Edges",
                                                        "Outerplanar Core Out Edges"},
                                                       {std::to_string(core_nodes.Len()),
                                                        std::to_string((double)graph.nodes()/core_nodes.Len()),
                                                        std::to_string(0),
                                                        std::to_string(0),
                                                        std::to_string(0)}, counter);
                            }
                            overall_outerplanar =
                                    outerplanar_sampling_runtime + conversion_runtime + outerplanar_core_runtime;
                            //std::cout << "Outerplanar Core Runtime: " << outerplanar_core_runtime << "s" << std::endl;
                            //std::cout << "Outerplanar Runtime: " << overall_outerplanar << "s" << std::endl;
                            //GraphFunctions::analyse_graph(core_graph.get_graph(), "Approximation Core Outerplanar");
                            eval.headerValueInsert({"Outerplanar Core Runtime", "Outerplanar Runtime"},
                                                   {std::to_string(outerplanar_core_runtime),
                                                    std::to_string(overall_outerplanar)}, counter, true, true);
                            if (exact_computation) {
                                int intersection_length = core_nodes.IntrsLen(exact_core);
                                int union_length = core_nodes.UnionLen(exact_core);
                                double recall = (double) intersection_length / exact_core.Len();
                                double similarity = (double) intersection_length / union_length;
                                //std::cout << "Outerplanar Recall: " << recall << std::endl;
                                //std::cout << "Outerplanar Similarity: " << similarity << std::endl;
                                //Approx values
                                eval.headerValueInsert(
                                        {"Outerplanar Recall",
                                         "Outerplanar Similarity"},
                                        {
                                                std::to_string(recall),
                                                std::to_string(similarity)}, counter,
                                        true, true);
                            } else {
                                eval.headerValueInsert(
                                        {"Outerplanar Recall",
                                         "Outerplanar Similarity"},
                                        {
                                                std::to_string(0),
                                                std::to_string(0)}, counter,
                                        true, true);
                            }
                            ++counter;
                        }
                    }
                }
            }
            if (params.outerplanar_statistics) {
                statistic.evaluate(out_stat_headers, out_stat_values);
            }
            //std::cout << "Outerplanar Sampling Runtime: " << runtime << "s" << std::endl;

//            int counter = 0;
//            for (double threshold: params.threshold) {
//                for (int generator_size: params.generator_size) {
//                    //Sampling param values
//                    eval.headerValueInsert(
//                            {"Generators", "Core Iterations", "Input Set Seed", "Threshold", "Overall Threshold",
//                             "Samples", "Sampling Seed"},
//                            {std::to_string(generator_size), std::to_string(params.coreIterations),
//                             std::to_string(params.generator_seed), std::to_string(threshold),
//                             std::to_string(params.overall_threshold),
//                             std::to_string(params.sample_number), std::to_string(params.sample_seed)}, counter);
//
//                    //Sampling param values
//                    eval.headerValueInsert(
//                            {"Tree Sampling Runtime", "Outerplanar Sampling Runtime", "Outerplanar Conversion Runtime"},
//                            {
//                                    std::to_string(tree_sampling_runtime),
//                                    std::to_string(outerplanar_sampling_runtime),
//                                    std::to_string(conversion_runtime)}, counter,true, true);
//                }
//            }
        }
        double simple_approx_runtime = 0;
        if (params.simple_approx) {
            int counter = 0;
            for (double threshold: params.threshold) {
                for (int generator_size: params.generator_size) {
                    //std::cout << std::endl << "Simple Approx Core:" << std::endl;
                    simple_approx(graph, generator_size, params.coreIterations, params.generator_seed,
                                  core_degree_distribution, core_nodes, simple_approx_runtime,
                                  outerplanar_core_runtime +
                                  (outerplanar_sampling_runtime + conversion_runtime) / params.sample_number);
                    GraphFunctions::GetCoreGraphStats(core_graph, graph, core_nodes, simple_approx_runtime,
                                                      "Simple Approx Core", &eval);
                    //std::cout << "Simple Approx Runtime: " << runtime << "s" << std::endl;
                    //GraphFunctions::analyse_graph(core_graph.get_graph(), "Simple approx graph");
                    eval.headerValueInsert({"Simple Approx Runtime"}, {std::to_string(simple_approx_runtime)}, counter, true,
                                           true);
                    if (exact_computation) {
                        int intersection_length = core_nodes.IntrsLen(exact_core);
                        int union_length = core_nodes.UnionLen(exact_core);
                        double recall = (double) intersection_length / exact_core.Len();
                        double similarity = (double) intersection_length / union_length;
                        //std::cout << "Simple Approx Recall: " << recall << std::endl;
                        //std::cout << "Simple Approx Similarity: " << similarity << std::endl;
                        //Approx values
                        eval.headerValueInsert(
                                {"Simple Approx Recall", "Simple Approx Similarity"},
                                {std::to_string(recall), std::to_string(similarity)}, counter,
                                true, true);
                    } else {
                        eval.headerValueInsert(
                                {"Simple Approx Recall", "Simple Approx Similarity"},
                                {std::to_string(0), std::to_string(0)}, counter,
                                true, true);
                    }
                    ++counter;
                }
            }
        }
        if (params.outerplanar_statistics) {
            eval.headerValueInsert(out_stat_headers, out_stat_values, -1, true, true);
        }


        graph.get_graph().Clr();
        core_graph.get_graph().Clr();
#pragma omp critical
        eval.save(true, true);
    }

}



int main(int argc, char *argv[]) {
    std::string out_path = "../out/ICDE_2022/";
    std::string in_path = "../../GraphData/";
    int max_nodes = std::numeric_limits<int>::max();
    int max_edges = std::numeric_limits<int>::max();

    std::vector<int> generators = {10};
    int generator_seed = 0;
    int core_iterations = 3;
    int samples = 100;
    std::vector<double> threshold = {1.0 / 100};
    double overall_threshold = 3.0 / 100;
    bool overall = false;
    int sample_seed = 1114352;
    bool outerplanar_new = true;
    int thread_num = omp_get_max_threads();
    bool set_out_path = false;
    bool outerplanar_statistics = false;
    bool core_stats = false;
    bool percentage = false;
    bool simple_approx_bool = false;

    int generator_begin = 5;
    int generator_next = 10;
    int generator_end = 500;
    int generator_step = 100;

    double threshold_begin = 0.01;
    double threshold_end = 0.04;
    double threshold_step = 0.01;
    bool tree_eval = true;


    for (int i = 0; i < argc; ++i) {
        if (std::strcmp(argv[i], "-i") == 0){
            in_path += std::string(argv[i+1]);
        }
        if (std::strcmp(argv[i], "-o") == 0){
            out_path += std::string(argv[i+1]);
            set_out_path = true;
        }
        if (std::strcmp(argv[i], "--generators") == 0){
            generator_next = atoi(argv[i+1]);
        }
        if (std::strcmp(argv[i], "--generators_begin") == 0){
            generator_begin = atoi(argv[i+1]);
        }
        if (std::strcmp(argv[i], "--threshold") == 0){
            threshold_begin = atof(argv[i+1]);
        }
        if (std::strcmp(argv[i], "--generators_end") == 0){
            generator_end = atoi(argv[i+1]);
        }
        if (std::strcmp(argv[i], "--threshold_end") == 0){
            threshold_end = atof(argv[i+1]);
        }
        if (std::strcmp(argv[i], "--generators_step") == 0){
            generator_step = atoi(argv[i+1]);
        }
        if (std::strcmp(argv[i], "--threshold_step") == 0){
            threshold_step = atof(argv[i+1]);
        }
        if (std::strcmp(argv[i], "--generators_seed") == 0){
            generator_seed = atoi(argv[i+1]);
        }
        if (std::strcmp(argv[i], "--core_iterations") == 0){
            core_iterations = atoi(argv[i+1]);
        }
        if (std::strcmp(argv[i], "--samples") == 0){
            samples = atoi(argv[i+1]);
        }

        if (std::strcmp(argv[i], "--overall_threshold") == 0){
            overall_threshold = atof(argv[i+1]);
        }
        if (std::strcmp(argv[i], "--max_nodes") == 0){
            max_nodes = atoi(argv[i+1]);
        }
        if (std::strcmp(argv[i], "--max_edges") == 0){
            max_edges = atoi(argv[i+1]);
        }
        if (std::strcmp(argv[i], "--overall") == 0){
            overall = true;
        }
        if (std::strcmp(argv[i], "--no_tree") == 0){
            tree_eval = false;
        }
        if (std::strcmp(argv[i], "--sample_seed") == 0){
            sample_seed = atoi(argv[i+1]);
        }
        if (std::strcmp(argv[i], "--naive_outerplanar") == 0){
            outerplanar_new = false;
        }
        if (std::strcmp(argv[i], "--simple_approx") == 0){
            simple_approx_bool = true;
        }
        if (std::strcmp(argv[i], "--outerplanar_stats") == 0){
            outerplanar_statistics = true;
        }
        if (std::strcmp(argv[i], "--core_stats") == 0){
            core_stats = true;
        }
        if (std::strcmp(argv[i], "--threads") == 0){
            thread_num = atoi(argv[i+1]);
        }
        if (std::strcmp(argv[i], "--percentage") == 0){
            percentage = atoi(argv[i+1]);
        }

    }
    if (!set_out_path){
        out_path += "core_approximation";
    }
    generators.clear();
    threshold.clear();

    generators.emplace_back(generator_begin);

    for (int i = generator_next; i < generator_end + generator_step; i += generator_step) {
        int j = i - i % generator_step;
        if (j == 0){
            j = i;
        }
        generators.emplace_back(j);
    }
    int steps = (int) ((threshold_end - threshold_begin)/threshold_step) + 1;
    for (int i = 0; i < steps; ++i) {
        threshold.emplace_back(threshold_begin + i*threshold_step);
    }


    OverlapApproxParams params = {in_path, out_path, thread_num, max_nodes, max_edges, generators, generator_seed, core_iterations, samples,
                                  threshold, overall_threshold, overall, sample_seed, outerplanar_new,false, false,
                                  true, true, simple_approx_bool, 0, 0,
                                  core_stats, outerplanar_statistics, percentage, tree_eval};
    overlap_eval(params);
}
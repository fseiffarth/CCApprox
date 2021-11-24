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

    int max_nodes = 100000;
    int max_edges = 200000;

    int generator_size = 10;
    int generator_seed = 0;
    int coreIterations = 3;
    int sample_number = 200;
    double threshold = 1.0 / (double) this->sample_number;
    double overall_threshold = 1.0 / (double) this->sample_number;
    bool overall = false;
    int sample_seed = 32487643;

    //Eval
    bool periphery = false;
    bool core_iteration = false;
    bool exact = true;
    bool approx_core = false;

    bool simple_approx = false;
};



void get_core(GraphData& graph, int generator_size, int core_iterations, int seed, std::vector<int>& intersection_loss, TIntV& coreNodes, double& runtime){
    GraphClosureSP gc;
    ClosureParameters closureParameters;
    std::vector<std::set<NodeId>> closures;
    std::cout << std::endl;
    std::set<NodeId> overlap;
    runtime = 0;
    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < core_iterations; ++i) {
        StaticFunctions::generateInputSet(closureParameters.input_set, graph, generator_size, i + core_iterations * seed);
        gc.naive_closure(graph, closureParameters);
        closures.emplace_back(closureParameters.closed_set);
        if (i == 0) {
            overlap = closureParameters.closed_set;
        } else {
            int overlap_size = (int) overlap.size();
            std::vector<NodeId> v_intersection;
            std::set_intersection(overlap.begin(), overlap.end(),
                                  closures.back().begin(), closures.back().end(),
                                  std::back_inserter(v_intersection));
            overlap.clear();
            for (auto elem: v_intersection) {
                overlap.insert(elem);
            }
            intersection_loss.emplace_back(overlap_size - overlap.size());
        }
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

void core_calc(OverlapApproxParams& params) {
    std::vector<std::tuple<std::string, int>> paths_and_sizes;
    for (const auto &entry: std::filesystem::recursive_directory_iterator(params.in_path)) {
        if (entry.path().extension() == ".edges") {
            const GraphData& data = GraphData(entry.path().string());
            int nodes = data.nodes();
            int edges = data.edges();
            if (nodes < params.max_nodes && edges < params.max_edges) {
                paths_and_sizes.emplace_back(std::tuple<std::string, int>(entry.path().string(), edges));
            }
        }
    }
    std::sort(paths_and_sizes.begin(), paths_and_sizes.end(), sortbysecond);

#pragma omp parallel for default(none) shared(paths_and_sizes, params)
    for (const auto &[path, size]: paths_and_sizes) {
        std::string stripped_path = std::filesystem::path(path).parent_path().string() + "/" + std::filesystem::path(path).stem().string();

//        TIntV load_core;
//        StaticFunctions::load(stripped_path + ".core", load_core);

        FileEvaluation eval = FileEvaluation(stripped_path, 1, ".core_info");
        GraphData graph = GraphData(path);
        GraphFunctions::GetLargestComponent(graph);
        eval.headerValueInsert({"Graph", "Nodes", "Edges", "Density"}, {graph.getName(), std::to_string(graph.nodes()), std::to_string(graph.edges()), std::to_string(graph.density())});

        std::vector<int> intersection_loss;
        TIntV exact_core;
        double exact_runtime = 0;
        get_core(graph, params.generator_size, params.coreIterations, params.generator_seed, intersection_loss, exact_core, exact_runtime);
        StaticFunctions::save(stripped_path, exact_core);

        eval.headerValueInsert( {"Core Size", "Generators", "CoreIterations", "Intersection Loss", "Exact Runtime"},
                                {std::to_string(exact_core.Len()),
                                 std::to_string(params.generator_size),
                                 std::to_string(params.coreIterations),
                                 StaticFunctions::print<std::vector<int>, int>(intersection_loss),
                                 std::to_string(exact_runtime)});
        //GraphFunctions::analyse_graph(core_graph.get_graph(), "Core");
        // std::cout << "Exact Runtime: " << exact_runtime << "s" << std::endl;

        eval.save();
    }
}



int main(int argc, char *argv[]) {
    std::string out_path = "../out/Experiments/Overlap/";
    std::string in_path = "../../GraphData/";
    int max_nodes = std::numeric_limits<int>::max();
    int max_edges = std::numeric_limits<int>::max();

    int generators = 10;
    int generator_seed = 0;
    int core_iterations = 3;
    bool overall = false;
    int sample_seed = 1114352;
    omp_set_num_threads(omp_get_max_threads());


    for (int i = 0; i < argc; ++i) {
        if (std::strcmp(argv[i], "-i") == 0){
            in_path += std::string(argv[i+1]);
        }
        if (std::strcmp(argv[i], "-o") == 0){
            out_path += std::string(argv[i+1]);
        }
        if (std::strcmp(argv[i], "-t") == 0){
            omp_set_num_threads(atoi(argv[i+1]));
        }
        if (std::strcmp(argv[i], "--generators") == 0){
            generators = atoi(argv[i+1]);
        }
        if (std::strcmp(argv[i], "--generators_seed") == 0){
            generator_seed = atoi(argv[i+1]);
        }
        if (std::strcmp(argv[i], "--core_iterations") == 0){
            core_iterations = atoi(argv[i+1]);
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
        if (std::strcmp(argv[i], "--sample_seed") == 0){
            sample_seed = atoi(argv[i+1]);
        }
    }

    OverlapApproxParams params = {in_path, out_path, max_nodes, max_edges, generators,generator_seed, core_iterations, 0, 0, 0, overall, sample_seed,
                                  false, false,
                                  true, true, true};
    core_calc(params);
}
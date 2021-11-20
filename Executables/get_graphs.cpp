#include "../Data/GraphData.h"
#include "../Utils/FileEvaluation.h"
#include "../Utils/GraphFunctions.h"

//
// Created by florian on 14.11.21.
//
int main() {
    std::string in_path = "../../GraphData/";
    for (const auto &entry: std::filesystem::recursive_directory_iterator(in_path)) {
        if (entry.path().extension() == ".edges") {
            FileEvaluation fileEvaluation = FileEvaluation("../out/ICDE_2022/graphs");
            std::string stripped_path =
                    entry.path().parent_path().string() + "/" + entry.path().stem().string();
            std::cout << "Consider " << stripped_path + ".edges" << std::endl;
            GraphData data = GraphData(stripped_path + ".edges");

            size_t nodes = data.nodes();
            size_t edges = data.edges();
            double density = (double) edges / ((double) (nodes * (nodes - 1)) / 2);
            double exact_runtime = 0;
            TIntV exact_core;
            StaticFunctions::load(stripped_path + ".core", exact_core);
            fileEvaluation.headerValueInsert({"Graph", "Size", "Edges", "Density"},
                                             {data.getName(),
                                              std::to_string(nodes),
                                              std::to_string(edges),
                                              std::to_string(density)
                                             });
            if (exact_core.Len() > 0) {
                std::vector<std::vector<std::string>> info_array;
                StaticFunctions::load_csv(stripped_path + ".core_info", info_array);
                exact_runtime = std::stoi(info_array[1][9]);
                GraphData core_graph = GraphData(TSnap::GetSubGraph(data.get_graph(), exact_core));
                GraphFunctions::GetCoreGraphStats(core_graph, data, exact_core, exact_runtime, "Exact Core", &fileEvaluation);
            }
            else{
                fileEvaluation.headerValueInsert(
                        {"Exact Core Nodes", "Exact Core Relative Nodes", "Exact Core Edges", "Exact Core Relative Edges", "Exact Core Out Edges", "Exact Core Runtime"},
                        {"0", "0", "0", "0", "0", "0"});
            }

            fileEvaluation.save();
            std::cout << "Finished " << stripped_path + ".edges" << std::endl;
        }
    }
}
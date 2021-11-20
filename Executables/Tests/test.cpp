//
// Created by florian on 29.09.21.
//

#include <iostream>
#include "../../Data/GraphData.h"
#include "../../Utils/StaticFunctions.h"
#include "../../Utils/GraphFunctions.h"
#include "../../ClosureOperators/GraphClosures.h"
#include "../../Experiments/DS2021/SetupExperiment.h"



void graph_eval(){
    std::vector<std::string> graph_paths;
    for (const auto &entry: std::filesystem::directory_iterator("../../GraphData/")) {
        if (entry.path().extension() == ".edges" && GraphData(entry.path().string()).get_data()->GetNodes() < 100000) {
            GraphData graph = GraphData(entry.path().string());
            GraphFunctions::GetLargestComponent(graph);
            std::map<int, int> degree_distribution;
            graph.getDegreeDistribution(degree_distribution);
            std::cout << graph.getName() << std::endl;
            std::cout << "Size: " << graph.size() << std::endl;
            std::cout << StaticFunctions::printMap(degree_distribution) << std::endl;
            std::cout << std::endl;
        }
    }
}

void createSynthetic(){
    int size = 5000;
    TRnd Rnd;
    for (int i = 1; i <= 100; ++i) {
        int edges = (int) ((double) (size*(size-1))/2 * ((double) i/2000));
        Rnd.PutSeed(i);
        GraphData graphData = GraphData(TSnap::GenRndGnm<PUNGraph>(size, edges, false, Rnd), "synthetic_nodes_" + std::to_string(size) + "_edges_" + std::to_string(edges));
        graphData.save_edges("../../GraphData/Synthetic/");
    }
}

void update(){
    for (const auto &entry: std::filesystem::recursive_directory_iterator("../../GraphData/")) {
        if (entry.path().extension() == ".edges" || entry.path().extension() == ".txt") {
            GraphData::Update(entry.path().string());
        }
    }
}

void component(){
    for (const auto &entry: std::filesystem::recursive_directory_iterator("../../GraphData/RealWorld/")) {
        if (entry.path().extension() == ".edges" || entry.path().extension() == ".txt") {
            GraphData graph = GraphData(entry.path().string());
            int size = graph.nodes();
            GraphFunctions::GetLargestComponent(graph);
            if (graph.nodes() < size) {
                std::string name = entry.path().stem().string();
                graph.setName(name + "_component");
                std::string parent_path = entry.path().parent_path().string();
                graph.save_edges(parent_path + "/");
            }
        }
    }
}

int main() {
    //update();
    component();
    //createSynthetic();
}
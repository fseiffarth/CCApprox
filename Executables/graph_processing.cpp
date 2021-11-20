//
// Created by florian on 31.08.21.
//

#include "../Data/GraphData.h"
#include "../Utils/GraphFunctions.h"

void take_graph_samples(const std::string &string, int size, int seed = 0);
void normalize_graphs();
void create_connected_components();

int main() {
    //normalize_graphs();
    //create_connected_components();
    take_graph_samples("../../../GraphData/", 200000, 0);
}

void take_graph_samples(const std::string &string, int size, int seed) {
    std::vector<GraphData> graphData;
    for (const auto & entry : std::filesystem::directory_iterator("../../GraphData/")) {
        if (entry.path().extension() == ".edges") {
            GraphData data = GraphData(entry.path().string());
            if (data.get_data()->GetNodes() > size) {
                if (!TSnap::IsConnected(data.get_data())) {
                    std::cout << "Graph " << data.getName() << " is not connected!" << std::endl;
                    std::cout << "Get biggest connected component!" << std::endl;
                    TCnComV ConComps;
                    TSnap::GetSccs(data.get_data(), ConComps);
                    TCnCom *biggestComponent = nullptr;
                    int Components = ConComps.Len();
                    for (int i = 0; i < Components; ++i) {
                        if (biggestComponent == nullptr || biggestComponent->NIdV.Len() < ConComps[i].NIdV.Len()) {
                            biggestComponent = &ConComps[i];
                        }
                    }
                    GraphData biggestComponentGraph = GraphData(TSnap::GetSubGraph(data.get_data(), biggestComponent->NIdV, true));
                    biggestComponentGraph.setName(data.getName());
                    if (biggestComponentGraph.get_data()->GetNodes() > size) {
                        graphData.emplace_back(biggestComponentGraph);
                    }
                }
                else{
                    graphData.emplace_back(data);
                }
            }
        }
    }
    std::mt19937_64 gen(seed);
    std::vector<GraphData> subGraphs;
    for (auto const& graph : graphData) {
        TIntV NodeIds = TIntV();
        GraphFunctions::getNodesInBall(graph.get_data(), NodeIds, gen,  size);
        subGraphs.emplace_back(GraphData(TSnap::GetSubGraph(graph.get_data(), NodeIds, true)));
        subGraphs.back().setName(graph.getName() + "_seed_" + std::to_string(seed) + "_size_" + std::to_string(subGraphs.back().size()));
    }
    for (auto const& graph : subGraphs) {
        graph.save_edges("../../GraphData/");
    }
}

void normalize_graphs() {
    auto path = "../../../GraphData/";
    for (const auto & entry : std::filesystem::directory_iterator(path)){
        if (entry.path().extension() == ".txt"){
            GraphData graphData = GraphData(entry.path().string());
            GraphData newGraph = GraphData(GraphFunctions::ResetGraphIds(graphData.get_data()));
            std::filesystem::path p = entry.path();
            newGraph.save_edges(p.replace_extension().string());
        }
    }

}


void create_connected_components() {
    auto path = "../../../GraphData/";
    std::vector<std::filesystem::directory_entry> paths;
    for (const auto & entry : std::filesystem::directory_iterator(path)){
        paths.emplace_back(entry);
    }
    for (const auto & entry : paths){
        if (entry.path().extension() == ".edges"){
            GraphData graphData = GraphData(entry.path().string());
            TCnComV ConComps;
            TSnap::Get1CnCom(graphData.get_data(), ConComps);
            if (ConComps.Len() > 1) {
                std::filesystem::path p = entry.path();
                for (int i = 0; i < ConComps.Len(); ++i) {
                    GraphData newGraph = GraphData(TSnap::GetSubGraph(graphData.get_data(), ConComps[i].NIdV));
                    std::string newPath = p.replace_extension().string() + "_" + std::to_string(i);
                    newGraph.save(newPath);
                }
            }
        }
    }
}



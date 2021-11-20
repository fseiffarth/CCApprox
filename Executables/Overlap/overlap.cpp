//
// Created by florian on 30.09.21.
//

#include <iostream>
#include "../../Data/GraphData.h"
#include "../../ClosureOperators/GraphClosures.h"
#include "../../ClosureOperators/BaseOperator.h"
#include "../../Utils/StaticFunctions.h"
#include "../../Experiments/DS2021/SetupExperiment.h"

void get_core(GraphData& graph, GraphData& core_graph, int generator_size, int iterations, std::map<int, int>& degree_distribution, TIntV& coreNodes){
    GraphClosureSP gc;
    ClosureParameters closureParameters;
    std::vector<std::set<NodeId>> closures;
    graph.getDegreeDistribution(degree_distribution);
    std::cout << graph.getName() << std::endl;
    std::cout << "Size: " << graph.size() << std::endl;
    std::cout << "Edges: " << graph.get_graph()->GetEdges() << std::endl;
    std::cout << StaticFunctions::printMap(degree_distribution) << std::endl;
    std::cout << std::endl;

    std::set<NodeId> overlap;
    for (int i = 0; i < iterations; ++i) {
        SetupExperiment::generateInputSet(closureParameters.input_set, graph, generator_size, i);
        gc.naive_closure(graph, closureParameters);
        //gc.fast_closure(graph, closureParameters);
        std::cout << "Closure Size: " << closureParameters.closed_set.size() << std::endl;
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
        std::cout << "Overlap Size: " << overlap.size() << std::endl;
    }
    TIntV nodesInOverlap;
    for (auto elem: overlap) {
        nodesInOverlap.Add(elem);
    }
    core_graph = GraphData(TSnap::GetSubGraph(graph.get_graph(), nodesInOverlap));

    std::cout << "Overlap Subgraph: " << std::endl;
    std::cout << "Size: " << core_graph.size() << std::endl;
    std::cout << "Edges: " << core_graph.get_graph()->GetEdges() << std::endl;

    int overlapOutEdges = 0;

    for(auto NodeIt = graph.get_graph()->BegNI(); NodeIt != graph.get_graph()->EndNI(); NodeIt++){
        NodeId Id = NodeIt.GetId();
        if(!core_graph.get_graph()->IsNode(Id)) {
            int degree = NodeIt.GetDeg();
            for (int n = 0; n < degree; ++n) {
                int NeighborId = NodeIt.GetNbrNId(n);
                if (core_graph.get_graph()->IsNode(NeighborId)){
                    ++overlapOutEdges;
                }
            }
        }
    }

    graph.getDegreeDistribution(degree_distribution, &overlap);
    std::cout << "Overlap degree distribution" << std::endl;
    std::cout << StaticFunctions::printMap(degree_distribution) << std::endl;
    std::cout << std::endl;
    std::cout << " Overlap out edges: " << std::endl;
    std::cout << overlapOutEdges << std::endl;
}


void overlap_eval(bool periphery = false, bool core_iteration = false) {
    std::vector<std::string> graph_paths;
    for (const auto &entry: std::filesystem::directory_iterator("../../GraphData/")) {
        if (entry.path().extension() == ".edges" && GraphData(entry.path().string()).get_data()->GetNodes() < 50000) {
            graph_paths.emplace_back(entry.path().string());
        }
    }

    //graph_paths = {"../../GraphData/CA-CondMat.edges"};
    for (const auto &path: graph_paths) {
        GraphData graph = GraphData(path);
        GraphFunctions::GetLargestComponent(graph);
        std::map<int, int> degree_distribution;

        GraphData core_graph;
        TIntV core_nodes;
        get_core(graph, core_graph, 10,3, degree_distribution, core_nodes);
        GraphFunctions::analyse_graph(core_graph.get_graph(), "Core");
        if(periphery) {
            PUNGraph periphery_graph = GraphFunctions::RemoveSubgraph(graph.get_graph(), core_graph.get_graph());
            GraphFunctions::analyse_graph(periphery_graph, "Periphery");
        }
        if(core_iteration) {
            GraphData core2_graph;
            core_graph.set_graph(GraphFunctions::ResetGraphIds(core_graph.graph()));
            get_core(core_graph, core2_graph, 10, 3, degree_distribution, core_nodes);
        }
    }
}

int main() {
    overlap_eval();
}
//
// Created by florian on 31.08.21.
//

#include <vector>
#include <iostream>
#include "../../Utils/typedefs.h"
#include "../../Data/GraphData.h"
#include "../../Utils/OuterplanarSubgraphDFS.h"
#include "../../Utils/GraphFunctions.h"

static bool bugFind(int& missingEdgeNum ,const std::pair<int, int>& size, int seed = 0, bool with_missing_edges = false){
    TRnd rand;
    rand.PutSeed(seed);
    auto graph = TSnap::GenRndGnm<PUNGraph>(size.first, size.second, false, rand);
    if (TSnap::IsConnected(graph)) {
        std::vector<NodePair> missingEdges;
        GraphData data = GraphData(graph, "TestGraph");
        OuterplanarSubgraphDFS outerplanarSubgraphDfs = OuterplanarSubgraphDFS(data.get_graph());
        GraphData subgraph = GraphData(new TUNGraph());
        std::mt19937_64 gen(seed);
        outerplanarSubgraphDfs.generate(subgraph, gen, true);
        bool out = !GraphFunctions::IsOuterPlanar(subgraph.get_graph());
        if (with_missing_edges) {
            GraphFunctions::IsMaximalOuterplanarSubgraph(data.get_graph(), subgraph.get_graph(), missingEdges);
            missingEdgeNum += (int) missingEdges.size();
        }
        std::cout << data.getName() << std::endl;
        data.print();
        std::cout << "Seed: " << seed << std::endl;
        std::cout << "Edges: " << subgraph.edges() << " / Missing: " << missingEdges.size();
        for (auto& edge : missingEdges) {
            std::cout << "(" << edge.first() << "," << edge.second() << ")";
        }
        std::cout << std::endl;
        subgraph.print();
        if(out){
            std::cout << "Incorrect!" << std::endl;
            data.save_dot("../Executables/Tests/");
        }
        else{
            std::cout << "Outerplanar!" << std::endl;
        }
        std::cout << std::endl;
        return out;

    }
    return false;
}

int main() {
    int missingEdges = 0;
    int runs = 1000;
    bool with_missing_edges = true;
    std::pair<int, int> size = {20, 50};
    bugFind(missingEdges, size, 62, with_missing_edges);
    for (int i = 0; i < runs; ++i) {
        if(bugFind(missingEdges, size, i, with_missing_edges)){
            break;
        }
    }
    if (with_missing_edges) {
        std::cout << "Missing Edges: " << missingEdges << std::endl;
        std::cout << "Relative Missing Edges: " << (double) missingEdges / (runs * size.second) << std::endl;
    }
//    closureTest();
//    closureTestB();
//    counterExample();
//    graphTest();
}
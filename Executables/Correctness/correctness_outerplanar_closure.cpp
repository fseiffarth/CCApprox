//
// Created by florian on 31.08.21.
//

#include <vector>
#include <iostream>
#include "../../Utils/typedefs.h"
#include "../../Data/GraphData.h"
#include "../../Utils/OuterplanarSubgraphDFS.h"
#include "../../Utils/GraphFunctions.h"
#include "../../ClosureOperators/GraphClosures.h"
#include "../../Utils/StaticFunctions.h"

struct Times{

    double gen_time;
    double gen_time_new;
    double outerplanar;
    double outerplanar_new;
    void print(){
        std::cout <<  "Subgraph Generation: " << (double) gen_time / 1000000 << std::endl;
        std::cout <<  "Subgraph To Outerplanar Graph: " << (double) gen_time_new / 1000000 << std::endl;
        std::cout <<  "Outerplanar Closure: " << (double) outerplanar / 1000000 << std::endl;
        std::cout <<  "Outerplanar New Closure: " << (double) outerplanar_new / 1000000 << std::endl;
    };
};

static bool bugFind(int& missingEdgeNum ,const std::pair<int, int>& size, Times& times, int seed = 0, bool with_missing_edges = false){
    TRnd rand;
    rand.PutSeed(seed);
    auto graph = TSnap::GenRndGnm<PUNGraph>(size.first, size.second, false, rand);
    if (TSnap::IsConnected(graph)) {
        std::vector<NodePair> missingEdges;
        GraphData data = GraphData(graph, "TestGraph");
        auto start = std::chrono::high_resolution_clock::now();
        OuterplanarSubgraphDFS outerplanarSubgraphDfs = OuterplanarSubgraphDFS(data.get_graph());
        GraphData subgraph = GraphData(new TUNGraph());
        std::mt19937_64 gen(seed);
        outerplanarSubgraphDfs.generate(subgraph, gen, true);
        subgraph.init();
        times.gen_time += std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start).count();
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

        start = std::chrono::high_resolution_clock::now();
        OuterplanarGraphData outSub = OuterplanarGraphData(subgraph);
        times.gen_time_new += std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start).count();

        GraphClosureSP cl = GraphClosureSP();
        ClosureParameters closureParameters1;
        ClosureParameters closureParameters2;
        StaticFunctions::generateInputSet(closureParameters1.input_set, subgraph, 10, seed);
        closureParameters2.input_set = closureParameters1.input_set;

        //subgraph.graphType = GraphType::GENERAL;
        start = std::chrono::high_resolution_clock::now();
        cl.naive_closure(subgraph, closureParameters1);
        times.outerplanar += std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start).count();

        start = std::chrono::high_resolution_clock::now();
        cl.naive_closure(outSub, closureParameters2);
        times.outerplanar_new += std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start).count();

        if (closureParameters1.closed_set != closureParameters2.closed_set){
            std::cout << "Not the same closure!" << std::endl;
            std::cout << StaticFunctions::print<std::set<NodeId>, NodeId>(closureParameters1.closed_set);
            std::cout << std::endl;
            std::cout << StaticFunctions::print<std::set<NodeId>, NodeId>(closureParameters2.closed_set);
            std::cout << std::endl;
            out = true;
        }

        return out;

    }
    return false;
}

int main() {
    int missingEdges = 0;
    int runs = 1000;
    bool with_missing_edges = false;
    std::pair<int, int> size = {1000, 10000};
    Times times{};
    bugFind(missingEdges, size, times, 12765, with_missing_edges);
    for (int i = 0; i < runs; ++i) {
        if(bugFind(missingEdges, size, times, i, with_missing_edges)){
            break;
        }
    }
    if (with_missing_edges) {
        std::cout << "Missing Edges: " << missingEdges << std::endl;
        std::cout << "Relative Missing Edges: " << (double) missingEdges / (runs * size.second) << std::endl;
    }
    times.print();
}
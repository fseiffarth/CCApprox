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


int main() {
    PUNGraph graph = new TUNGraph();
    for (int i = 0; i < 12; ++i) {
        graph->AddNode();
    }
    graph->AddEdge(10, 0);
    graph->AddEdge(0, 1);
    graph->AddEdge(1, 2);
    graph->AddEdge(1, 3);
    graph->AddEdge(3, 4);
    graph->AddEdge(3, 5);
    graph->AddEdge(4, 5);
    graph->AddEdge(5, 6);
    graph->AddEdge(5, 7);
    graph->AddEdge(6, 7);
    graph->AddEdge(5, 9);
    graph->AddEdge(5, 8);
    graph->AddEdge(8, 9);

    graph->AddEdge(4, 11);
    graph->AddEdge(8, 11);
    graph->AddEdge(5, 11);

//    for (int i = 0; i < 9; ++i) {
//        graph->AddEdge(i % 9, (i+1) % 9);
//    }
//    graph->AddEdge(1, 8);
//    graph->AddEdge(2, 8);
//    graph->AddEdge(3, 5);
//    graph->AddEdge(2, 6);
//    std::vector<PUNGraph> faces;
//
//    GraphFunctions::GetBiconnectedOuterplanarFaces(graph, faces);

    OuterplanarGraphData biconn = OuterplanarGraphData(graph);
    GraphClosureSP cl = GraphClosureSP();
    ClosureParameters closureParameters;
    closureParameters.input_set = {0, 11};
    cl.naive_closure(biconn, closureParameters);
    closureParameters.print();

}
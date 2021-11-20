//
// Created by Florian on 08.05.2021.
//

#ifndef CLOSURES_CORRECTNESS_CLOSURE_H
#define CLOSURES_CORRECTNESS_CLOSURE_H

#include <iostream>
#include "../../Utils/typedefs.h"
#include "../../Data/GraphData.h"
#include "../../Utils/Graphs.h"
#include "../../Utils/graphio.h"
#include "../../ClosureOperators/GraphClosures.h"
#include "../../Experiments/DS2021/SetupExperiment.h"
#include "../../Utils/StaticFunctions.h"

static bool comparison(GraphData& graph, ClosureParameters& output, std::vector<double>& runtimes){
    GraphClosureSP gc = GraphClosureSP("ShortestPath");
    std::vector<std::set<NodeId>> closedSets;
    double duration = 0;

    std::cout << "Naive Closure: " << std::endl;
    auto start = std::chrono::system_clock::now();
    gc.naive_closure(graph, output);
    closedSets.emplace_back(output.closed_set);
    std::cout << std::endl;
    duration = (double) std::chrono::duration_cast<std::chrono::microseconds>(
            std::chrono::system_clock::now() - start).count() / 1000000;
    std::cout << "Duration: " << duration << std::endl;
    output.print();
    std::cout << std::endl << std::endl;
    runtimes[0] += duration;

//    std::cout << "First Closure: " << std::endl;
//    start = std::chrono::system_clock::now();
//    gc.closure(graph, output);
//    closedSets.emplace_back(output.closed_set);
//    std::cout << std::endl;
//    std::cout << "Size: " << output.closed_set.size() << std::endl;
//    std::cout << "Duration: " << (double) std::chrono::duration_cast<std::chrono::microseconds>(
//            std::chrono::system_clock::now() - start).count() / 1000000 << std::endl;
//    output.print();
//    std::cout << std::endl << std::endl;


    std::cout << "Fast Closure: " << std::endl;
    start = std::chrono::system_clock::now();
    gc.fast_closure(graph, output);
    closedSets.emplace_back(output.closed_set);
    std::cout << std::endl;
    duration = (double) std::chrono::duration_cast<std::chrono::microseconds>(
            std::chrono::system_clock::now() - start).count() / 1000000;
    std::cout << "Duration: " << duration << std::endl;;
    output.print();
    std::cout << std::endl << std::endl;
    runtimes[1] += duration;

    std::set<NodeId> overlap = closedSets.front();
    for (int i = 1; i < closedSets.size(); ++i) {
        std::vector<NodeId> v_intersection;
        std::set_intersection(overlap.begin(), overlap.end(),
                              closedSets[i].begin(), closedSets[i].end(),
                              std::back_inserter(v_intersection));
        overlap.clear();
        for (auto elem: v_intersection) {
            overlap.insert(elem);
        }
    }
    bool allEqual = true;
    for (auto const & set : closedSets) {
        if (overlap.size() != set.size()){
            allEqual = false;
            break;
        }
    }
    if (allEqual){
        std::cout << "Closures coincide !!!" << std::endl;
    }
    else{
        std::cout << "There is some bug!!!" << std::endl;
    }
    return allEqual;
}

static bool bugFind(std::vector<double>& runtimes, int seed = 0){
    TRnd rand;
    rand.PutSeed(seed);
    auto graph = TSnap::GenRndGnm<PUNGraph>(14, 18, false, rand);
    if (TSnap::IsConnected(graph)) {
        GraphData data = GraphData(graph, "TestGraph");
        std::set<NodeId> input_set;
        StaticFunctions::generateInputSet(input_set, data, 3, seed);
        ClosureParameters output = ClosureParameters(input_set);

        bool out = !comparison(data, output, runtimes);
        data.print();
        std::cout << "Seed: " << seed << std::endl;
        if(out){
            data.save_dot("../Executables/Tests/");
        }
        return out;

    }
    return false;
}

static void gridTest(){
    PUNGraph graph = new TUNGraph();
    int n = 10;
    int m = 10;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {

        }
    }
}

static void closureTest(){
    PUNGraph graph = Graphs::completeBipartite(2, 3);
    GraphData data = GraphData(graph, "TestGraph");
    std::set<NodeId> input_set = std::set <NodeId> {2, 0};
    ClosureParameters output = ClosureParameters(input_set);
    std::vector<double> runtimes;
    comparison(data, output, runtimes);
}

static void graphTest(){
    GraphData data = GraphData("../../GraphData/CA-CondMat.edges");
    std::set<NodeId> input_set;
    StaticFunctions::generateInputSet(input_set, data, 10, 0);
    ClosureParameters output = ClosureParameters(input_set);
    std::vector<double> runtimes;
    comparison(data, output, runtimes);
}

static void counterExample() {
    PUNGraph g = new TUNGraph;
    int n = 1000;
    for (int i = 0; i < 2 * n + 1; ++i) {
        g->AddNode();
    }
    for (int i = 0; i < n; ++i) {
        int j = 2 * i;
        g->AddEdge(j, j + 2);
        g->AddEdge(j + 1, j + 2);
        if (i != n - 1) {
            g->AddEdge(j, j + 3);
            g->AddEdge(j + 1, j + 3);
        }
    }

    GraphData graph = GraphData(g);

    GraphClosureSP gc = GraphClosureSP("ShortestPath");
    std::set<NodeId> input_set = std::set<NodeId>{0, 1};
    ClosureParameters output = ClosureParameters(input_set);
    std::vector<double> runtimes;
    comparison(graph, output, runtimes);
}

static void closureTestB(){
    TRnd rnd;
    int n = 10;
    int e = 20;
    auto graph = TSnap::GenRndGnm<PUNGraph>(n, e, false, rnd);
    int i = 1;
    while (!TSnap::IsConnected(graph)){
        rnd.PutSeed(i);
        graph = TSnap::GenRndGnm<PUNGraph>(n, e, false, rnd);
        ++i;
    }
    GraphData data = GraphData(graph, "TestGraph");
    GraphClosureSP gc = GraphClosureSP("ShortestPath");
    std::set<NodeId> input_set = std::set <NodeId> {0, 1, 6};
    ClosureParameters output = ClosureParameters(input_set);
    std::vector<double> runtimes;
    comparison(data, output, runtimes);
}

static void extendedClosureTest() {
    GraphData data = GraphData("../out/graphs/OneComponent_A_5_dA_20_B_5_dB_20_C_0_0.possibleEdges");
    GraphClosureSP gc = GraphClosureSP("ShortestPath");
    std::set<NodeId> input_set = std::set <NodeId> {0,1, 2, 3, 4};
    ClosureParameters output = ClosureParameters(input_set);
    gc.closure(data, output);
    output.print();
}

static void checkGeneratedGraphs() {
    std::vector<GraphData> graphs;
    std::string graphPath = "../out/graphs/";
    LoadGraphData(graphPath, graphs, LoadProperties({{"TwoComponents", -1}}, -1));
    GraphClosureSP gc = GraphClosureSP("ShortestPath");

    for (auto& graph : graphs) {
        int size = graph.size();
        for (auto edge = graph.get_data()->BegEI(); edge != graph.get_data()->EndEI(); edge++) {
            if (edge.GetSrcNId() < size/2 && edge.GetDstNId() >= size/2){
                std::cout << "(" << edge.GetSrcNId() << "," << edge.GetDstNId() << ")" << " ";
            }
        }
        std::cout << std::endl;
        std::cout << graph.getName() << std::endl;
        std::set<NodeId> input_setA = std::set <NodeId>();
        std::set<NodeId> input_setB = std::set <NodeId>();

        for (int i = 0; i<size/2; ++i) {
            input_setA.insert(i);
            input_setB.insert(i+size/2);
        }
        ClosureParameters outputA = ClosureParameters(input_setA);
        gc.closure(graph, outputA);
        if(outputA.closed_set == input_setA){
            std::cout << "Generated Graph is _valid!" << std::endl;
        }
        else{
            std::cout << "Generated Graph is invalid!"<< std::endl;
        }
        outputA.print();
        ClosureParameters outputB = ClosureParameters(input_setB);
        gc.closure(graph, outputB);
        if(outputB.closed_set == input_setB){
            std::cout << "Generated Graph is _valid!" << std::endl;
        }
        else{
            std::cout << "Generated Graph is invalid!" << std::endl;
        }
        outputB.print();
    }
}

#endif //CLOSURES_CORRECTNESS_CLOSURE_H

//
// Created by Florian on 10.05.2021.
//

#ifndef CLOSURES_TESTS_H
#define CLOSURES_TESTS_H



//static void OuterplanarTest(){
//    bool outerplanar;
//    std::mt19937_64 gen(0);
//    PUNGraph graph = GraphFunctions::maximalOuterplanarSubgraphLinear(
//            GraphData("../out/tests/graphs/Graph_A_5_dA_20_B_5_dB_20_C_3_0.possibleEdges").get_data(), gen);
//    outerplanar = GraphFunctions::IsOuterPlanar(graph);
//    outerplanar = GraphFunctions::IsOuterPlanar(Graphs::triangle());
//    outerplanar = GraphFunctions::IsOuterPlanar(Graphs::completeBipartite(2, 3));
//    outerplanar = GraphFunctions::IsOuterPlanar(Graphs::componentGraph({Graphs::circle(4), Graphs::path(3), Graphs::circle(3)}, {{0, 0}, {2, 0}}));
//    outerplanar = GraphFunctions::IsOuterPlanar(Graphs::componentGraph({Graphs::circle(4), Graphs::path(3), Graphs::completeBipartite(2, 3)}, {{0, 0}, {2, 0}}));
//    outerplanar = GraphFunctions::IsOuterPlanar(
//            GraphData("../out/tests/outerplanar/Graph_A_5_dA_20_B_5_dB_20_C_3_0_Pattern0.possibleEdges").get_data());
//    outerplanar = GraphFunctions::IsOuterPlanar(
//            GraphData("../out/tests/graphs/Graph_A_5_dA_20_B_5_dB_20_C_3_0.possibleEdges").get_data());
//
//}
//
//static void OuterplanarDFSTest(){
//    int nonOuterPlanarNumber = 0;
//    std::vector<int> nonOuterplanarSeeds;
//    int nonMaximalNumber = 0;
//    std::vector<int> nonMaximalSeeds;
//    PUNGraph graph1 = Graphs::circle(5);
//    graph1->AddEdge(0, 2);
//    graph1->AddEdge(1, 3);
//    PUNGraph graph2 = Graphs::circle(5);
//    graph2->AddEdge(0, 3);
//    graph2->AddEdge(2, 4);
//
//    PUNGraph graph = Graphs::connectGraphs({graph1, graph2}, {{0, 2}, {3, 0}});
//    GraphFunctions::print(graph, "Graph");
//    int seedMin = 0;
//    int seedMax = 9999;
//    for(int seed = seedMin; seed < seedMax + 1; ++seed) {
//        std::mt19937_64 gen(seed);
//        std::cout << "Seed: " << seed << " //////////////////////////////////////" << std::endl;
//        GraphData outerplanarSubgraph = GraphData(new TUNGraph());
//        OuterplanarSubgraphDFS subgraphGeneration = OuterplanarSubgraphDFS(graph);
//        subgraphGeneration.generate(outerplanarSubgraph, gen, true);
//        std::vector<int> nonMaximalEdgeNum;
//        std::vector<int> outerplanarEdges;
//        std::vector<int> outerplanarMissingEdges;
//        checkingOuterplanarity(graph, outerplanarSubgraph.get_graph(), nonOuterPlanarNumber, nonMaximalNumber, nonOuterplanarSeeds, nonMaximalSeeds, nonMaximalEdgeNum, outerplanarMissingEdges, seed);
//    }
//    std::cout << "Graphs: " << seedMax - seedMin + 1 << std::endl;
//    std::cout << " Not Outerplanar: " << nonOuterPlanarNumber << " " << StaticFunctions::print<std::vector<int>, int>(nonOuterplanarSeeds) << std::endl;
//    std::cout << " Not Maximal: " << nonMaximalNumber << " " << StaticFunctions::print<std::vector<int>, int>(nonMaximalSeeds) << std::endl;
//}
//
//static void OuterplanarMitchellTest(){
//    int nonOuterPlanarNumber = 0;
//    std::vector<int> nonOuterplanarSeeds;
//    int nonMaximalNumber = 0;
//    std::vector<int> nonMaximalSeeds;
//    PUNGraph graph1 = Graphs::circle(5);
//    graph1->AddEdge(0, 2);
//    graph1->AddEdge(1, 3);
//    PUNGraph graph2 = Graphs::circle(5);
//    graph2->AddEdge(0, 3);
//    graph2->AddEdge(2, 4);
//
//    const PUNGraph graph = Graphs::connectGraphs({graph1, graph2}, {{0, 2}, {3, 0}});
//    GraphFunctions::print(graph, "Graph:");
//    int seedMin = 0;
//    int seedMax = 10000;
//    for(int seed = seedMin; seed < seedMax + 1; ++seed) {
//        std::mt19937_64 gen(seed);
//        std::cout << "Seed: " << seed << " //////////////////////////////////////" << std::endl;
//        GraphData outerplanarSubgraph = GraphData(new TUNGraph());
//        OuterPlanarSubgraphMitchell subgraphGeneration = OuterPlanarSubgraphMitchell(graph);
//        subgraphGeneration.generate(outerplanarSubgraph, gen, true);
//        GraphFunctions::print(outerplanarSubgraph.get_graph(), "Outerplanar Subgraph:");
//        std::vector<int> nonMaximalEdgeNum;
//        std::vector<int> outerplanarEdges;
//        std::vector<int> outerplanarMissingEdges;
//        checkingOuterplanarity(graph, outerplanarSubgraph.get_graph(), nonOuterPlanarNumber, nonMaximalNumber, nonOuterplanarSeeds, nonMaximalSeeds, nonMaximalEdgeNum, outerplanarMissingEdges, seed);
//    }
//    std::cout << "Graphs: " << seedMax - seedMin + 1 << std::endl;
//    std::cout << " Not Outerplanar: " << nonOuterPlanarNumber << " " << StaticFunctions::print<std::vector<int>, int>(nonOuterplanarSeeds) << std::endl;
//    std::cout << " Not Maximal: " << nonMaximalNumber << " " << StaticFunctions::print<std::vector<int>, int>(nonMaximalSeeds) << std::endl;
//}
//
//static void OuterPlanarRandomTest(const std::string& out_path, const std::string &algorithm,  TestParameter& testParameter) {
//    std::vector<int> neighbors = std::vector<int>();
//    int threads = testParameter.threads;
//    if (threads == -1){
//        threads = omp_get_max_threads();
//    }
//    omp_set_num_threads(threads);
//    for (auto const &graph_info: testParameter.graph_sizes) {
//        int nonOuterPlanarNumber = 0;
//        std::vector<int> nonOuterPlanarSeeds;
//        int nonMaximalNumber = 0;
//        std::vector<int> nonMaximalSeeds;
//        TRnd rnd = TInt::Rnd;
//
//        std::vector<int> maximalEdges;
//        double nonMaximalEdgeNumPerCent = 0.0;
//        int outerPlanarEdgeNum = 0;
//        int maximalEdgesNum = 0;
//        std::vector<int> algorithmMissingEdges;
//        int missingEdgesNum = 0;
//        size_t outerplanar_runtime = 0;
//        size_t tree_runtime = 0;
//        OuterplanarGraphStatistics privateStatistics = OuterplanarGraphStatistics();
//        OuterplanarGraphStatistics statistics = OuterplanarGraphStatistics();
//        int iterations = testParameter.seedMax + 1 - testParameter.seedMin;
//        rnd.PutSeed(0);
//        for (int seed = testParameter.seedMin; seed < testParameter.seedMax + 1; ++seed) {
//            auto graph = TSnap::GenRndGnm<PUNGraph>(graph_info.first, graph_info.second, false, rnd);
//            if (TSnap::IsConnected(graph)) {
//                std::cout << "Seed: " << seed << std::endl;
//                if (seed >= testParameter.seedMin) {
//                    std::mt19937_64 gen(seed);
//                    for (auto type: {PatternType::BFS_TREE, PatternType::OUTERPLANAR}) {
//                        GraphData subgraph = GraphData(new TUNGraph(), graph->GetNodes());
//                        if (type == PatternType::BFS_TREE) {
//                            GraphFunctions::generateNeighborVector(graph, neighbors);
//                            auto start = std::chrono::high_resolution_clock::now();
//                            GraphFunctions::bfsSubtree(graph, subgraph, neighbors, gen);
//                            tree_runtime += std::chrono::duration_cast<std::chrono::microseconds>(
//                                    std::chrono::high_resolution_clock::now() - start).count();
//                            if (testParameter.printLevel <= PrintLevel::ALL_GRAPHS) {
//                                GraphFunctions::print(graph, "Graph");
//                                //std::cout << "Seed: " << generator_seed << " //////////////////////////////////////" << std::endl;
//                            }
//                        } else{
//                            auto start = std::chrono::high_resolution_clock::now();
//                            OuterplanarSubgraphDFS subgraphGeneration = OuterplanarSubgraphDFS(graph);
//                            subgraphGeneration.generate(subgraph, gen, false);
//                            outerplanar_runtime += std::chrono::duration_cast<std::chrono::microseconds>(
//                                    std::chrono::high_resolution_clock::now() - start).count();
//                            privateStatistics += OuterplanarGraphStatistics(subgraph.get_graph());
//                            outerPlanarEdgeNum += subgraph.get_graph()->GetEdges();
//                            if (testParameter.outerplanarity_check) {
//                                checkingOuterplanarity(graph, subgraph.get_graph(), nonOuterPlanarNumber, nonMaximalNumber,
//                                                       nonOuterPlanarSeeds,
//                                                       nonMaximalSeeds, algorithmMissingEdges, maximalEdges,
//                                                       seed, testParameter.printLevel);
//                                maximalEdgesNum += (subgraph.get_graph()->GetEdges() + algorithmMissingEdges.back());
//                                nonMaximalEdgeNumPerCent += (
//                                        (double) (subgraph.get_graph()->GetEdges() - algorithmMissingEdges.back()) /
//                                        (double) subgraph.get_graph()->GetEdges());
//                                missingEdgesNum += algorithmMissingEdges.back();
//                            }
//                            if (testParameter.printLevel <= PrintLevel::ALL_GRAPHS) {
//                                GraphFunctions::print(subgraph.get_graph(), "Out Graph");
//                            }
//                        }
//                    }
//                }
//            }
//            else {
//                std::cout << "Not connected!" << std::endl;
//                --seed;
//            }
//            graph.Clr();
//        }
//        privateStatistics /= iterations;
//        statistics += privateStatistics;
//
//        std::vector<std::string> stat_headers;
//        std::vector<std::string> stat_values;
//        statistics.evaluate(stat_headers, stat_values);
//
//        std::vector<std::string> headers = {"Size",
//                                            "Edges",
//                                            "Density",
//                                            "Algorithm",
//                                            "Outerplanar Time",
//                                            "Tree Time",
//                                            "Factor",
//                                            "Samples",
//                                            "Not Outerplanar Samplings",
//                                            "Not Maximal Samplings",
//                                            "Average Maximal Edges",
//                                            "Average Edges",
//                                            "Missing Edges per Graph",
//                                            "% of Maximal outerplanar graph"};
//
//        std::vector<std::string> values = {std::to_string(graph_info.first),
//                                           std::to_string(graph_info.second),
//                                           std::to_string(graph_info.second / ((double) graph_info.first/2 * (graph_info.first - 1))),
//                                           algorithm,
//                                           std::to_string((double) outerplanar_runtime / (iterations * 1000000 * threads)),
//                                           std::to_string((double) tree_runtime / (iterations * 1000000 * threads)),
//                                           std::to_string((double) outerplanar_runtime / (double) tree_runtime),
//                                           std::to_string(iterations),
//                                           std::to_string(nonOuterPlanarNumber),
//                                           std::to_string(nonMaximalNumber),
//                                           std::to_string(maximalEdgesNum / (double) iterations),
//                                           std::to_string(outerPlanarEdgeNum / (double) iterations),
//                                           std::to_string(missingEdgesNum / (double) iterations),
//                                           std::to_string(nonMaximalEdgeNumPerCent / iterations * 100.0)};
//
//        headers.insert(headers.end(), stat_headers.begin(), stat_headers.end());
//        values.insert(values.end(), stat_values.begin(), stat_values.end());
//
//        StaticFunctions::saveValuesToFile(out_path, headers,values);
//        std::cout << std::endl << "******************************* FileEvaluation " << algorithm
//                  << " *********************************" << std::endl;
//        std::cout << "Runtime: " << (double) outerplanar_runtime << " microseconds" << ", Factor: "
//                  << (double) outerplanar_runtime / (double) tree_runtime << std::endl;
//        std::cout << "BFSTime: " << (double) tree_runtime << " microseconds" << std::endl;
//        std::cout << "Nodes: " << graph_info.first << std::endl;
//        std::cout << "Edges: " << graph_info.second << std::endl;
//        std::cout << "Density: " << graph_info.second / ((double) graph_info.first/2 * (graph_info.first - 1))<< std::endl;
//        std::cout << "Samples: " << testParameter.seedMax - testParameter.seedMin + 1 << std::endl;
//        std::cout << "\tNot Outerplanar: " << nonOuterPlanarNumber << " "
//                  << StaticFunctions::print<std::vector<int>, int>(nonOuterPlanarSeeds) << std::endl;
//        std::cout << "\tNot Maximal: " << nonMaximalNumber << " "
//                  << StaticFunctions::print<std::vector<int>, int>(nonMaximalSeeds) << std::endl;
//        std::cout << "Average Edges: " << outerPlanarEdgeNum / (double) iterations << std::endl;
//        std::cout << "Average Maximal Edges: "
//                  << std::accumulate(algorithmMissingEdges.begin(), algorithmMissingEdges.end(), 0.0) /
//                     (double) algorithmMissingEdges.size() << std::endl;
//        std::cout << "Missing Edges per Graph: "
//                  << std::accumulate(maximalEdges.begin(), maximalEdges.end(), 0.0) /
//                     (testParameter.seedMax - testParameter.seedMin + 1) << std::endl;
//        std::cout << "Edges per Graph in %: "
//                  << nonMaximalEdgeNumPerCent /(double) iterations * 100.0 << std::endl;
//    }
//}
//
//static void OuterplanarFaceTest() {
//    ExperimentalSetup setup("new_results", "../out/results/", "../out/faceTest/graphs/", "../out/faceTest/trees/",
//                            "../out/faceTest/outerplanar/", 1, {50},
//                            {2}, {10}, {1},
//                            5, 1, 1);
//    GraphClosureSP gc("shortestPath");
//    //Graph and pattern generation
//    setup.generateGraphs(10, GenerationType::TWO_COMPONENTS);
//    std::cout << "Finished graph two component generation!" << std::endl;
//    setup.generatePatterns(gc, 10, 5, GenerationType::TWO_COMPONENTS, PatternType::OUTERPLANAR);
//    std::cout << "Finished pattern two component generation!" << std::endl;
//    setup.analyzeFaces();
//}
//
//static void test_run(){
//    GraphClosureSP gc("shortestPath");
//    InitUpdateGreedy<GraphData> initUpdateGreedy("greedy");
//    InitUpdateGreedyRandom<GraphData> initUpdateGreedyRandom("greedyRandom");
//
//    ExperimentalSetup setup("test_results", "../out/tests/results/", "../out/tests/graphs/", "../out/tests/trees/",
//                            "../out/tests/outerplanar/", 1, {5},
//                            {2}, {}, {1, 2, 4},
//                            5, 10, 1);
//    //Graph and pattern generation
//    setup.generateGraphs(10, GenerationType::TWO_COMPONENTS);
//    std::cout << "Finished graph two component generation!" << std::endl;
//    setup.generatePatterns(gc, 10, 5, GenerationType::TWO_COMPONENTS, PatternType::BFS_TREE);
//    std::cout << "Finished pattern two component generation!" << std::endl;
//    setup.generatePatterns(gc, 10, 5, GenerationType::TWO_COMPONENTS, PatternType::OUTERPLANAR);
//    std::cout << "Finished pattern two component generation!" << std::endl;

    //Graph and pattern generation
//    setup.generateGraphs(10, GenerationType::ONE_COMPONENT);
//    std::cout << "Finished graph two component generation!" << std::endl;
//    setup.generatePatterns(gc, 10, 5, GenerationType::ONE_COMPONENT, PatternType::BFS_TREE);
//    std::cout << "Finished pattern two component generation!" << std::endl;
//    setup.generatePatterns(gc, 10, 5, GenerationType::ONE_COMPONENT, PatternType::OUTERPLANAR);
//    std::cout << "Finished pattern two component generation!" << std::endl;

//    //Algorithm run
//    for (int numPatterns : {5}) {
//        setup.run(gc, initUpdateGreedyRandom, 5, numPatterns, 1, GenerationType::TWO_COMPONENTS, PatternType::BFS_TREE);
//        setup.run(gc, initUpdateGreedy, 5, numPatterns, 1, GenerationType::ONE_COMPONENT, PatternType::BFS_TREE);
//    }
//    for (int numPatterns : {5}) {
//        setup.run(gc, initUpdateGreedyRandom, 5, numPatterns, 1, GenerationType::TWO_COMPONENTS, PatternType::OUTERPLANAR);
//        setup.run(gc, initUpdateGreedy, 5, numPatterns, 1, GenerationType::ONE_COMPONENT, PatternType::OUTERPLANAR);
//    }
//}

#endif //CLOSURES_TESTS_H

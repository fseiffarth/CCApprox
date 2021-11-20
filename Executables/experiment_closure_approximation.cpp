//
// Created by florian on 31.08.21.
//


#include "../Runs/ApproximationsRun.h"

int main() {
    run_approximations_test("../../GraphData/", "../out/Experiments/Approximations/test_pc.csv", 0, 3, 1);
    //run_approximations_test("", "../out/Experiments/Approximations/synthetic.csv", 0, 3, 1, true);
    //run_approximations_real_graphs_small("../out/Experiments/Approximations/real_graphs_small.csv");
    //run_approximations_road_subgraphs("../out/Experiments/Approximations/road_subgraphs.csv");
    //run_approximations_synthetic_graphs("../out/Experiments/Approximations/synthetic_graphs.csv");
}
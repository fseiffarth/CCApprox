//
// Created by florian on 31.08.21.
//

#include "../Experiments/Experiments.h"

void quality_experiment(){
    int n = 100;
    omp_set_num_threads(1);
#pragma omp parallel for default(none) shared(n)
    for (int j = 0; j < 10; ++j) {
        if (j % 4 == 3) {
            double p = 0.05 + 0.01 * j;
            std::vector<std::pair<int, int>> graph_sizes;
            for (int i = 1; i <= 5; ++i) {
                graph_sizes.emplace_back(std::pair<int, int>{i * n, (int) (((i * n) * ((i * n) - 1)) / 2) * p});
            }
            auto params = Experiments::TestParameter{graph_sizes, 0, 99, true, 1};
            std::string out_name = "sampling_quality";
            Experiments::OuterplanarSampling("../out/ICDE_2022/" + out_name, params);
        }
    }
}

int main(){
    quality_experiment();
}
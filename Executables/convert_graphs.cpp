//
// Created by Florian on 24.11.2021.
//

#include <string>
#include <filesystem>
#include <iostream>
#include "../Utils/FileEvaluation.h"
#include "../Data/GraphData.h"
#include "../Utils/GraphFunctions.h"

int main(int argc, char *argv[]) {
    std::string in_path = "../../GraphData/";
    for (int i = 0; i < argc; ++i) {
        if (std::strcmp(argv[i], "-i") == 0) {
            in_path += std::string(argv[i + 1]);
        }
    }

    for (const auto &entry: std::filesystem::recursive_directory_iterator(in_path)) {
        std::string extension;
        if (entry.path().extension() == ".edges"){
            extension = ".edges";
        }
        else if(entry.path().extension() == ".txt"){
            extension = ".txt";
        }
        if (!extension.empty()) {
            FileEvaluation fileEvaluation = FileEvaluation("../out/ICDE_2022/Other/converted_graphs");
            std::string stripped_path =
                    entry.path().parent_path().string() + "/" + entry.path().stem().string();
            std::cout << "Consider " << stripped_path + extension << std::endl;
            GraphData::Update(stripped_path + ".txt");
            GraphData data = GraphData(stripped_path + ".edges");
            GraphFunctions::GetLargestComponent(data);
            data.setName(entry.path().stem().string());
            data.save_edges(entry.path().parent_path().string() + "/");

            size_t nodes = data.nodes();
            size_t edges = data.edges();
            double density = (double) edges / ((double) (nodes * (nodes - 1)) / 2);
            fileEvaluation.headerValueInsert({"Graph", "Size", "Edges", "Density"},
                                             {data.getName(),
                                              std::to_string(nodes),
                                              std::to_string(edges),
                                              std::to_string(density)
                                             });
            fileEvaluation.save();
            std::cout << "Finished conversion of " << stripped_path + extension << std::endl;
        }
    }
}
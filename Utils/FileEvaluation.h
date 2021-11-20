//
// Created by florian on 21.10.21.
//

#ifndef CLOSURES_FILEEVALUATION_H
#define CLOSURES_FILEEVALUATION_H


#include <string>
#include <vector>
#include <iomanip>

class FileEvaluation {
    std::string out_path;
    std::vector<std::string> headers_all;
    std::vector<std::vector<std::string>> values_all;
    std::vector<std::string> headers_summary;
    std::vector<std::vector<std::string>> values_summary;
    std::string extension = ".csv";
public:
    explicit FileEvaluation(const std::string& out_path, int number = 1, const std::string& extension = ".csv") : out_path(out_path), extension(extension){
        for (int i = 0; i < number; ++i) {
            values_all.emplace_back(std::vector<std::string>());
            values_summary.emplace_back(std::vector<std::string>());
        }
    };
    void save(bool summary = false, bool both = true, std::_Ios_Openmode mode = std::ios_base::app);
    void headerValueInsert(const std::vector<std::string>& new_header, const std::vector<std::string>& new_values, int pos=0, bool summary = false, bool both = false);
    void clear();

};


#endif //CLOSURES_FILEEVALUATION_H

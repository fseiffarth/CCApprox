//
// Created by florian on 21.10.21.
//

#include <filesystem>
#include <fstream>
#include <numeric>
#include "FileEvaluation.h"
void FileEvaluation::save(bool summary, bool both, std::_Ios_Openmode mode) {
    if (summary && both){
        save(false, false, mode);
    }
    std::string appendix = extension;
    if (summary){
        appendix = "_summary" + extension;
    }
    bool newFile = std::filesystem::exists(this->out_path + appendix);
    std::ofstream fs;
    fs.open(this->out_path + appendix, mode);

    // write the file headers
    std::vector<std::string>* header = &headers_all;

    if (summary){
        header = &headers_summary;
    }
    if (!newFile) {
        for (const auto& string : *header) {
            fs << "," << string;
        }
        fs << std::endl;
    }
    for (int i = 0; i < values_all.size(); ++i) {
        std::vector<std::string> *values = &values_all[i];
        if (summary) {
            values = &values_summary[i];
        }
        fs << std::fixed;
        for (const auto &string: *values) {
            fs << "," << string;
        }
        fs << std::endl;
        fs << std::scientific;
    }

    fs.close();
}

void FileEvaluation::headerValueInsert(const std::vector<std::string> &new_header, const std::vector<std::string> &new_values,int pos, bool summary, bool both) {
    if (new_header.size() != new_values.size()){
        throw std::length_error("Header and Value lengths do not fit!");
    }
    std::vector<int> positions;
    if (pos == -1){
        positions = std::vector<int>(values_all.size());
        std::iota(positions.begin(),  positions.end(), 0);
    } else{
        positions = {pos};
    }


    int counter = 0;
    for (int x : positions) {
        std::vector<std::string> *header = nullptr;
        if (counter == 0 && pos < 1) {
            header = &headers_all;
        }
        std::vector<std::string> *values = &values_all[x];
        if (both) {
            if (counter == 0 && pos < 1) {
                header->insert(header->end(), new_header.begin(), new_header.end());
            }
            values->insert(values->end(), new_values.begin(), new_values.end());
            if (summary) {
                if (counter == 0 && pos < 1) {
                    header = &headers_summary;
                }
                values = &values_summary[x];
            }
            if (counter == 0 && pos < 1) {
                header->insert(header->end(), new_header.begin(), new_header.end());
            }
            values->insert(values->end(), new_values.begin(), new_values.end());
        } else {
            if (summary) {
                if (counter == 0 && pos < 1) {
                    header = &headers_summary;
                }
                values = &values_summary[x];
            }
            if (counter == 0 && pos < 1) {
                header->insert(header->end(), new_header.begin(), new_header.end());
            }
            values->insert(values->end(), new_values.begin(), new_values.end());
        }
        ++counter;
    }
}

void FileEvaluation::clear() {
    headers_all.clear();
    headers_summary.clear();
    values_all.clear();
    values_summary.clear();
}

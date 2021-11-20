#include "correctness_closure.h"



//
// Created by florian on 28.09.21.
//
int main() {
    std::vector<double> runtimes = {0, 0};
    bugFind(runtimes, 154);
//    for (int i = 0; i < 100000; ++i) {
//        if(bugFind(runtimes, i)){
//            break;
//        }
//    }
    std::cout << "Naive closure time: " << runtimes[0] << std::endl;
    std::cout << "Fast closure time: " << runtimes[1] << std::endl;
//    closureTest();
//    closureTestB();
//    counterExample();
//    graphTest();
}
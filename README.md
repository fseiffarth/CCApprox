# CCApprox

**Results of Experiments:**

-see *out/ICDE_2022/*

**Usage of Code:**

1. Clone this repository
2. Create folders *GraphData* and *ExternalLibraries* in the same location
3. Download the SNAP 6.0 library from http://snap.stanford.edu/snap/download.html and unpack in *ExternalLibraries* folder
4. Download the graphs from http://snap.stanford.edu/data/index.html and save unpacked *.txt* in *GraphData* folder or sub-folder
5. Create build folder and compile code with:
   ```mkdir build``` ```cd build``` ```cmake ..``` ```make -o2 -j 12```
6. Convert the graphs: ```./ConvertGraphs```
7. Compile and run the executables:
   1. ```./ExpSamplingRuntime``` for sampling runtime experiment
   2. ```./ExpSamplingQuality``` for sampling quality experiment
   3. ```./ExpClosureRuntime``` for closure runtime experiment
   4. ```./ExpExactCore``` for exact core computation
    
       | Optional Arguments | -i  | -o  | -t  | --generators | --generator_seed | --core_iterations | --max_nodes | --max_edges | --overall | --sample_seed |
       | :---:   | :-: | :-: | :-: | :----------: | :--------------: | :---------------: | :---------: | :---------: | :-------: | :-----------: |
       | Seconds | 301 | 283 | 301 | 283          | 301              | 283               | 301         | 283         | 301       | 283           |
   5. ```./ExpApproxCore``` for approximate core computation
   
      | Attempt | #1  | #2  |
      | :---:   | :-: | :-: |
      | Seconds | 301 | 283 |


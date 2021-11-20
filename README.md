# CCApprox

**Results of Experiments:**

-see *out/ICDE_2022/*

**Usage of Code:**

1. Download the graph data from http://snap.stanford.edu/data/index.html
2. Download the SNAP 6.0 library from http://snap.stanford.edu/snap/download.html
3. Convert the graphs
4. Compile and run the executables:

    4.1. ./ExpExactCore for exact core computation
    
    <| Attempt | #1  | #2  |>
    <| :---:   | :-: | :-: |>
    <| Seconds | 301 | 283 |>

    
    4.2. ./ExpApproxCore for approximate core computation
    
    <| Attempt | #1  | #2  |>
    <| :---:   | :-: | :-: |>
    <| Seconds | 301 | 283 |>

    
    4.3. ./ExpSamplingRuntime for sampling runtime experiment
    
    4.4. ./ExpSamplingQuality for sampling quality experiment
    
    4.5. ./ExpClosureRuntime for closure runtime experiment

# Installation

Please install the official Julia Binary from here https://julialang.org/downloads/. PLease download the Long-term support (LTS) release: v1.6.7 (July 19, 2022).

Now, open Julia REPL and input th following commands

``` 
import Pkg;
Pkg.add(url="https://github.com/pallavip11/SpinTrack.jl")
Pkg.add("Plots")

```
# Running Cases

To run the most basic example:

``` 
using SpinTrack;
using Plots;

m = FNAL_muon_ring()
sol = get_solution(u1_long(1e-6), m);
plot(sol)

```
